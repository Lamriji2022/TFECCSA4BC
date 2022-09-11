import sys
import random
import hashlib
import time
import collections


EllipticCurve = collections.namedtuple('EllipticCurve', 'name p a b g n h')

curve = EllipticCurve(
    'secp256k1',
    # Caractéristique du champ.
    p=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f,
    # Coefficients de courbe.
    a=0,
    b=7,
    # Point generateur.
    g=(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
       0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8),
    # Ordre des sous-groupes.
    n=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141,
    # Sous-groupe cofacteur.
    h=1,
)

# Inverse modulaire par la methode d'Euclid
def inverse_mod_Euclid(k, p):
    if k == 0:
        raise ZeroDivisionError('division par zero')

    if k < 0:
        return p - inverse_mod_Euclid(-k, p)

    # Algorithme euclidien étendu.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p

# Inverse modulaire par la methode d'Euler
def inverse_mod_Euler(k, p):
    if k == 0:
        raise ZeroDivisionError('division par zero')
    res = 1
    y = p-2
    k = k%p
    if k == 0:
        return 0
    while y > 0:
        if y%2 == 1:
            res = (res*k)%p
        y = y//2
        k = (k*k)%p
    return res

#fonction pour tester l'appartenance d'un point dans le curve
def is_on_curve(point):
    if point is None:
        # Aucun représente le point à l'infini.
        return True

    x, y = point

    return (y * y - x * x * x - curve.a * x - curve.b) % curve.p == 0

# la somme de deux points en utilisant 'l'inverse d'Euclid
def point_add_ECDSA1(point1, point2):
    if point1 is None:
        return point2
    if point2 is None:
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        return None

    if x1 == x2:
        m = (3 * x1 * x1 + curve.a) * inverse_mod_Euclid(2 * y1, curve.p)
    else:
        m = (y1 - y2) * inverse_mod_Euclid(x1 - x2, curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % curve.p, -y3 % curve.p)

    return result

# la somme de deux points en utilisant 'l'inverse d'Euler
def point_add_ECDSA2(point1, point2):
    if point1 is None:
        return point2
    if point2 is None:
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        return None

    if x1 == x2:
        m = (3 * x1 * x1 + curve.a) * inverse_mod_Euler(2 * y1, curve.p)
    else:
        m = (y1 - y2) * inverse_mod_Euler(x1 - x2, curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % curve.p, -y3 % curve.p)

    return result

# calculer -P
def point_neg(point):
    if point is None:
        return None

    x, y = point
    result = (x, -y % curve.p)

    return result

# le doubelement d'un point en utilisant 'l'inverse d'Euclid
def double_p_ECDSA1(point1):
    if point1 is None:
        return None

    x1, y1 = point1
    m = (3*x1*x1 + curve.a)*inverse_mod_Euclid(2*y1, curve.p)

    x3 = m * m - 2*x1
    y3 = m * (x1 - x3) - y1
    result = (x3 % curve.p, y3 % curve.p)

    return result

# le doubelement d'un point en utilisant 'l'inverse d'Euler
def double_p_ECDSA2(point1):
    if point1 is None:
        return None

    x1, y1 = point1
    m = (3*x1*x1 + curve.a)*inverse_mod_Euler(2*y1, curve.p)

    x3 = m * m - 2*x1
    y3 = m * (x1 - x3) - y1
    result = (x3 % curve.p, y3 % curve.p)

    return result

# Calcule de SM avec la methode doublement et addition en utilisant 'l'inverse d'Euclid
def mult_double_and_add_ECDSA1(k, point):
    Q = None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        k = list(bin(k)[2:])
        for i in range(len(k)):
            Q = double_p_ECDSA1(Q)
            if int(k[i]) == 1:
                Q = point_add_ECDSA1(Q,  point)
        return Q

# Calcule de SM avec la methode doublement et addition en utilisant 'l'inverse d'Euler
def mult_double_and_add_ECDSA2(k, point):
    Q = None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        k = list(bin(k)[2:])
        for i in range(len(k)):
            Q = double_p_ECDSA2(Q)
            if int(k[i]) == 1:
                Q = point_add_ECDSA2(Q,  point)
        return Q


# Calcule de SM avec la methode window en utilisant 'l'inverse d'Euclid
def mult_window_ECDSA1(k,w, point):
    assert is_on_curve(point)
    if k % curve.n == 0 or point is None:
        return None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        d=[]
        i=1
        
        while k>0:
            r = k%(2**(i*w))
            di = r//(2**((i-1)*w))
            d.append(di)
            k = k-r
            i = i+1
        
        Q=None
        R = mult_double_and_add_ECDSA1(pow(2,w),point)
        for j in range(len(d)):
            if j == 0:
                Q = mult_double_and_add_ECDSA1(d[j],point)
            else:
                if d[j] > 0:
                    Q = point_add_ECDSA1(Q, mult_double_and_add_ECDSA1(d[j],R)) 
                R = mult_double_and_add_ECDSA1(pow(2,w),R)
        return Q

# Calcule de SM avec la methode window en utilisant 'l'inverse d'Euler
def mult_window_ECDSA2(k,w, point):
    assert is_on_curve(point)
    if k % curve.n == 0 or point is None:
        return None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        d=[]
        i=1
        
        while k>0:
            r = k%(2**(i*w))
            di = r//(2**((i-1)*w))
            d.append(di)
            k = k-r
            i = i+1
        
        Q=None
        R = mult_double_and_add_ECDSA2(pow(2,w),point)
        for j in range(len(d)):
            if j == 0:
                Q = mult_double_and_add_ECDSA2(d[j],point)
            else:
                if d[j] > 0:
                    Q = point_add_ECDSA2(Q, mult_double_and_add_ECDSA2(d[j],R)) 
                R = mult_double_and_add_ECDSA2(pow(2,w),R)
        return Q

#Representer k en NAF
def naf(k):
    Naf = []
    i=0
    while k>0:
        if k%2 == 1:
            Naf.append(2 - k%4)
            k = k - Naf[i]
        else:
            Naf.append(0)
        k = k//2
        i = i+1
    return Naf

# Calcule de SM avec la methode NAF (Non Adjacent Form) en utilisant 'l'inverse d'Euclid
def mult_NAF_ECDSA1(k, point):
    assert is_on_curve(point)
    if k % curve.n == 0 or point is None:
        return None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        Naf=naf(k)
        
        Q = point
        for j in range(len(Naf)-2,-1,-1):
            Q = double_p_ECDSA1(Q)
            if Naf[j] == 1:
                Q = point_add_ECDSA1(Q,point)
            if Naf[j] == -1:
                Q = point_add_ECDSA1(Q,point_neg(point))
            
        return Q

# Calcule de SM avec la methode NAF en utilisant 'l'inverse d'Euler
def mult_NAF_ECDSA2(k, point):
    assert is_on_curve(point)
    if k % curve.n == 0 or point is None:
        return None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        Naf=naf(k)
        
        Q = point
        for j in range(len(Naf)-2,-1,-1):
            Q = double_p_ECDSA2(Q)
            if Naf[j] == 1:
                Q = point_add_ECDSA2(Q,point)
            if Naf[j] == -1:
                Q = point_add_ECDSA2(Q,point_neg(point))
            
        return Q

def mods(d,w):
    t=pow(2,w)
    if (d % t) >= t - 1:
        return (d % t) - t
    else:
        return d % t

#Representer k en W-NAF    
def w_naf(k,w):
    wNaf=[]
    i=0
    while k>0:
        if k%2 == 1:
            wNaf.append(mods(k,w))
            k = k - wNaf[i]
        else:
            wNaf.append(0)
        k = k//2
        i = i+1
    return wNaf

# Calcule de SM avec la methode Window NAF en utilisant 'l'inverse d'Euclid
def mult_w_Naf_ECDSA1(k,w,point):
    assert is_on_curve(point)
    if k % curve.n == 0 or point is None:
        return None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        wNaf = w_naf(k,w)
        
        R={}
        puis = pow(2,w)
        for l in range(1,puis):
            R[l]=mult_double_and_add_ECDSA1(l,point)
        Q = None
        for j in range(len(wNaf)-1,-1,-1):
            Q = double_p_ECDSA1(Q)
            if wNaf[j] > 0:
                Q = point_add_ECDSA1(Q,R[wNaf[j]])
            if wNaf[j] < 0:
                Q = point_add_ECDSA1(Q,point_neg(R[-wNaf[j]]))
            
        return Q

# Calcule de SM avec la methode Window NAF en utilisant 'l'inverse d'Euler
def mult_w_Naf_ECDSA2(k,w,point):
    assert is_on_curve(point)
    if k % curve.n == 0 or point is None:
        return None
    if k == 0:
        return None
    elif k == 1:
        return point
    else:
        wNaf = w_naf(k,w)
        
        R={}
        puis = pow(2,w)
        for l in range(1,puis):
            R[l]=mult_double_and_add_ECDSA2(l,point)
        if wNaf[len(wNaf)-1] > 0:
            Q = R[wNaf[len(wNaf)-1]]
        if wNaf[len(wNaf)-1] < 0:
            Q = point_neg(R[-wNaf[len(wNaf)-1]])
        for j in range(len(wNaf)-2,-1,-1):
            Q = double_p_ECDSA2(Q)
            if wNaf[j] > 0:
                Q = point_add_ECDSA2(Q,R[wNaf[j]])
            if wNaf[j] < 0:
                Q = point_add_ECDSA2(Q,point_neg(R[-wNaf[j]]))
            
        return Q


print("generation keys")
#Generer le paire de clés (dA,QA)
dA = random.randint(0, curve.n-1)
#  generer la paire de cles par la methode DBL et ADD avec inverseEuclid
iter = 1
st1 = time.time()
for i in range(iter):
    QA111 = mult_double_and_add_ECDSA1(dA,curve.g)
end1 = time.time()
print('DBL_ADD avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode DBL et ADD avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA11 = mult_double_and_add_ECDSA2(dA,curve.g)
end1 = time.time()
print('DBL_ADD avec InverseEuler  :',(end1 - st1)/iter)


#  generer la paire de cles par la methode window avec inverseEuclid
st1 = time.time()
for i in range(iter):
    QA2 = mult_window_ECDSA1(dA,4,curve.g)
end1 = time.time()
print('window avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode window avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA2 = mult_window_ECDSA2(dA,4,curve.g)
end1 = time.time()
print('window avec InverseEuler  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    QA3 = mult_NAF_ECDSA1(dA,curve.g)
end1 = time.time()
print('NAF avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode NAF avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA3 = mult_NAF_ECDSA2(dA,curve.g)
end1 = time.time()
print('NAF avec InverseEuler  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode W-NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    QA4 = mult_w_Naf_ECDSA1(dA,3,curve.g)
end1 = time.time()
print('W-NAF avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode W-NAF  avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA = mult_w_Naf_ECDSA2(dA,3,curve.g)
end1 = time.time()
print('W-NAF avec InverseEuler  :',(end1 - st1)/iter)


print("generation signature")

msg="Hello"
msg1="Hel lo"
# Signer le message


st2 = time.time()
h=int(hashlib.sha256(msg.encode()).hexdigest(),16)
k = random.randint(0, curve.n-1)
end2 = time.time()
t=end2 -st2
#  SM par la methode DBL et ADD avec inverseEuclid
st1 = time.time()
for i in range(iter):
    rpoint = mult_double_and_add_ECDSA1(k,curve.g)
end1 = time.time()
t1 = iter*t + (end1 - st1)

#  SM par la methode DBL et ADD avec inverseEuler
st1 = time.time()
for i in range(iter):
    rpoint = mult_double_and_add_ECDSA2(k,curve.g)
end1 = time.time()
t2 = iter*t + (end1 - st1)

#  SM par la methode window avec inverseEuclid
st1 = time.time()
for i in range(iter):
    rpoint = mult_window_ECDSA1(k,4,curve.g)
end1 = time.time()
t3 = iter*t + (end1 - st1)

#  SM par la methode window avec inverseEuler
st1 = time.time()
for i in range(iter):
    rpoint = mult_window_ECDSA2(k,4,curve.g)
end1 = time.time()
t4 = iter*t + (end1 - st1)

#  SM par la methode NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    rpoint = mult_NAF_ECDSA1(k,curve.g)
end1 = time.time()
t5 = iter*t + (end1 - st1)

#  SM par la methode NAF avec inverseEuler
st1 = time.time()
for i in range(iter):
    rpoint = mult_NAF_ECDSA2(k,curve.g)
end1 = time.time()
t6 = iter*t + (end1 - st1)

#  SM par la methode W-NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    rpoint = mult_w_Naf_ECDSA1(k,3,curve.g)
end1 = time.time()
t7 = iter*t + (end1 - st1)

#  SM par la methode W-NAF  avec inverseEuler
st1 = time.time()
for i in range(iter):
    rpoint = mult_w_Naf_ECDSA2(k,3,curve.g)
end1 = time.time()
t8 = iter*t + (end1 - st1)


st1 = time.time()
for i in range(iter):
    inv_k = inverse_mod_Euclid(k,curve.n)
end1 = time.time()
tEuc = end1 - st1

st1 = time.time()
for i in range(iter):
    inv_k = inverse_mod_Euler(k,curve.n)
end1 = time.time()
tEul = end1 - st1

st1 = time.time()
for i in range(iter):
    r = rpoint[0] % curve.n
    s = (inv_k*(h+r*dA)) % curve.n
end1 = time.time()
t = end1 - st1

print('DBL_ADD avec InverseEuclid  :',(t1+t+tEuc)/iter)
print('DBL_ADD avec InverseEuler  :',(t2+t+tEul)/iter)
print('window avec InverseEuclid  :',(t3+t+tEuc)/iter)
print('window avec InverseEuler  :',(t4+t+tEul)/iter)
print('NAF avec InverseEuclid  :',(t5+t+tEuc)/iter)
print('NAF avec InverseEuler  :',(t6+t+tEul)/iter)
print('W-NAF avec InverseEuclid  :',(t7+t+tEuc)/iter)
print('W-NAF avec InverseEuler  :',(t8+t+tEul)/iter)

print("verification sign")

# Verification de signature


st1 = time.time()
for i in range(iter):
    inv_s = inverse_mod_Euclid(s,curve.n)
end1 = time.time()
tEuc = end1 - st1

st1 = time.time()
for i in range(iter):
    inv_s = inverse_mod_Euler(s,curve.n)
end1 = time.time()
tEul = end1 - st1



st1 = time.time()
for i in range(iter):
    h=int(hashlib.sha256(msg.encode()).hexdigest(),16)
    u1=(h*inv_s) % curve.n
    u2=(r*inv_s) % curve.n
end1 = time.time()
t= end1 - st1

#  SM par la methode DBL et ADD avec inverseEuclid
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA1(mult_double_and_add_ECDSA1(u1,curve.g), mult_double_and_add_ECDSA1(u2,QA))
end1 = time.time()
t1= end1 - st1
#  SM par la methode DBL et ADD avec inverseEuler
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA2(mult_double_and_add_ECDSA2(u1,curve.g), mult_double_and_add_ECDSA2(u2,QA))
end1 = time.time()
t2= end1 - st1
#  SM par la methode window avec inverseEuclid
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA1(mult_window_ECDSA1(u1,4,curve.g), mult_window_ECDSA1(u2,4,QA))
end1 = time.time()
t3 =end1 - st1
#  SM par la methode window avec inverseEuler
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA2(mult_window_ECDSA2(u1,4,curve.g), mult_window_ECDSA2(u2,4,QA))
end1 = time.time()
t4 =end1 - st1
#  SM par la methode NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA1(mult_NAF_ECDSA1(u1,curve.g), mult_NAF_ECDSA1(u2,QA))
end1 = time.time()
t5 = end1 - st1

#  SM par la methode NAF avec inverseEuler
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA2(mult_NAF_ECDSA2(u1,curve.g), mult_NAF_ECDSA2(u2,QA))
end1 = time.time()
t6 = end1 - st1
#  SM par la methode W-NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA1(mult_w_Naf_ECDSA1(u1,3,curve.g), mult_w_Naf_ECDSA1(u2,3,QA))
end1 = time.time()
t7 = end1 - st1

#  SM par la methode W-NAF  avec inverseEuler
st1 = time.time()
for i in range(iter):
    P = point_add_ECDSA2(mult_w_Naf_ECDSA2(u1,3,curve.g), mult_w_Naf_ECDSA2(u2,3,QA))
end1 = time.time()
t8 = end1 - st1


st1 = time.time()
for i in range(iter):
    res = P[0] % curve.n
    #print (f"\nResult r={res}")
    if (res==r):
        sign = "valid"
        #print("Signature valid")
    else:
        sign = "invalid"
        #print("Signature invalid")
end1 = time.time()
tt = end1 - st1

print('DBL_ADD avec InverseEuclid  :',(t1+t+tt+tEuc)/iter)
print('DBL_ADD avec InverseEuler  :',(t2+t+tt+tEul)/iter)
print('window avec InverseEuclid  :',(t3+t+tt+tEuc)/iter)
print('window avec InverseEuler  :',(t4+t+tt+tEul)/iter)
print('NAF avec InverseEuclid  :',(t5+t+tt+tEuc)/iter)
print('NAF avec InverseEuler  :',(t6+t+tt+tEul)/iter)
print('W-NAF avec InverseEuclid  :',(t7+t+tt+tEuc)/iter)
print('W-NAF avec InverseEuler  :',(t8+t+tt+tEul)/iter)
