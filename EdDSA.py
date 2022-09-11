import sys
import random
import hashlib
import libnum
import time
import collections


EllipticCurve = collections.namedtuple('EllipticCurve', 'name p b a d g n')

curve = EllipticCurve(
    'ed25519',
    # Caractéristique du champ.
    p=(2**255) - 19,
    b = 256,
    # Coefficients de courbe.
    a=-1,
    d=37095705934669439343138083508754565189542113879843219016388785533085940283555,
    # Point generateur.
    g=(15112221349535400772501151409588531511454012693041857206046113283949847762202,46316835694926478169428394003475163141307993866256225615783033603165251855960),
    # Ordre des sous-groupes.
    n=(2**252) + 27742317777372353535851937790883648493,
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

    return (y*y + (curve.a*x*x) - 1 - (curve.d*x*x*y*y)) % curve.p == 0

# la somme de deux points en utilisant 'l'inverse d'Euclid
def point_add_EdDSA1(point1, point2):

    if point1 is None:
        return point2
    if point2 is None:
        return point1

    x1, y1 = point1
    x2, y2 = point2
    
    m1 = (x1*y2 + x2*y1)*inverse_mod_Euclid(1 + curve.d*x1*x2*y1*y2, curve.p)
    m2 = (y1*y2 - curve.a*x2*x1)*inverse_mod_Euclid(1 - curve.d*x1*x2*y1*y2, curve.p)
    
    result = (m1 % curve.p, m2 % curve.p)

    return result
# la somme de deux points en utilisant 'l'inverse d'Euler
def point_add_EdDSA2(point1, point2):
    if point1 is None:
        return point2
    if point2 is None:
        return point1

    x1, y1 = point1
    x2, y2 = point2

    m1 = (x1*y2 + x2*y1)*inverse_mod_Euler(1 + curve.d*x1*x2*y1*y2, curve.p)
    m2 = (y1*y2 - curve.a*x2*x1)*inverse_mod_Euler(1 - curve.d*x1*x2*y1*y2, curve.p)
    
    result = (m1 % curve.p, m2 % curve.p)
    return result

# calculer -P
def point_neg(point):

    if point is None:
        return None

    x, y = point
    result = (x, -y % curve.p)

    return result

# le doubelement d'un point en utilisant 'l'inverse d'Euclid
def double_p_EdDSA1(point1):
    x1, y1 = point1

    m1 = (2*x1*y1)*inverse_mod_Euclid(y1*y1 + curve.a*x1*x1, curve.p)
    m2 = (y1*y1 - curve.a*x1*x1)*inverse_mod_Euclid(2 - curve.a*x1*x1 - y1*y1, curve.p)
    
    result = (m1 % curve.p, m2 % curve.p)
    return result

# le doubelement d'un point en utilisant 'l'inverse d'Euler
def double_p_EdDSA2(point1):

    x1, y1 = point1

    m1 = (2*x1*y1)*inverse_mod_Euler(y1*y1 + curve.a*x1*x1, curve.p)
    m2 = (y1*y1 - curve.a*x1*x1)*inverse_mod_Euler(2 - curve.a*x1*x1 - y1*y1, curve.p)
    
    result = (m1 % curve.p, m2 % curve.p)
    return result

# Calcule de SM avec la methode doublement et addition en utilisant 'l'inverse d'Euclid
def mult_double_and_add_EdDSA1(k, point):
    if k == 0:
        return None
    elif k == 1:
        return point
    elif k%2 == 1:
        return point_add_EdDSA1(point, mult_double_and_add_EdDSA1(k - 1, point))
    else:
        return mult_double_and_add_EdDSA1(k//2, double_p_EdDSA1(point))

# Calcule de SM avec la methode doublement et addition en utilisant 'l'inverse d'Euler
def mult_double_and_add_EdDSA2(k, point):
    if k == 0:
        return None
    elif k == 1:
        return point
    elif k%2 == 1:
        return point_add_EdDSA2(point, mult_double_and_add_EdDSA2(k - 1, point))
    else:
        return mult_double_and_add_EdDSA2(k//2, double_p_EdDSA2(point))

# Calcule de SM avec la methode window en utilisant 'l'inverse d'Euclid    
def mult_window_EdDSA1(k,w, point):
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
        
        Q = (0,1)
        for j in range(len(d)-1,-1,-1):
            Q = mult_double_and_add_EdDSA1(2**w,Q)
            if d[j] > 0:
                Q = point_add_EdDSA1(Q, mult_double_and_add_EdDSA1(d[j],point)) 
        return Q
        
        
# Calcule de SM avec la methode window en utilisant 'l'inverse d'Euler
def mult_window_EdDSA2(k,w, point):
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
        
        Q = (0,1)
        for j in range(len(d)-1,-1,-1):
            Q = mult_double_and_add_EdDSA2(2**w,Q)
            if d[j] > 0:
                Q = point_add_EdDSA2(Q, mult_double_and_add_EdDSA2(d[j],point)) 
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
def mult_NAF_EdDSA1(k, point):
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
            Q = double_p_EdDSA1(Q)
            if Naf[j] == 1:
                Q = point_add_EdDSA1(Q,point)
            if Naf[j] == -1:
                Q = point_add_EdDSA1(Q,point_neg(point))
            
        return Q

# Calcule de SM avec la methode NAF en utilisant 'l'inverse d'Euler
def mult_NAF_EdDSA2(k, point):
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
            Q = double_p_EdDSA2(Q)
            if Naf[j] == 1:
                Q = point_add_EdDSA2(Q,point)
            if Naf[j] == -1:
                Q = point_add_EdDSA2(Q,point_neg(point))
            
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
def mult_w_Naf_EdDSA1(k,w,point):
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
            R[l]=mult_double_and_add_EdDSA1(l,point)
        if wNaf[len(wNaf)-1] > 0:
            Q = R[wNaf[len(wNaf)-1]]
        if wNaf[len(wNaf)-1] < 0:
            Q = point_neg(R[-wNaf[len(wNaf)-1]])
        for j in range(len(wNaf)-2,-1,-1):
            Q = double_p_EdDSA1(Q)
            if wNaf[j] > 0:
                Q = point_add_EdDSA1(Q,R[wNaf[j]])
            if wNaf[j] < 0:
                Q = point_add_EdDSA1(Q,point_neg(R[-wNaf[j]]))
            
        return Q


# Calcule de SM avec la methode Window NAF en utilisant 'l'inverse d'Euler
def mult_w_Naf_EdDSA2(k,w,point):
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
            R[l]=mult_double_and_add_EdDSA2(l,point)
        if wNaf[len(wNaf)-1] > 0:
            Q = R[wNaf[len(wNaf)-1]]
        if wNaf[len(wNaf)-1] < 0:
            Q = point_neg(R[-wNaf[len(wNaf)-1]])
        for j in range(len(wNaf)-2,-1,-1):
            Q = double_p_EdDSA2(Q)
            if wNaf[j] > 0:
                Q = point_add_EdDSA2(Q,R[wNaf[j]])
            if wNaf[j] < 0:
                Q = point_add_EdDSA2(Q,point_neg(R[-wNaf[j]]))
            
        return Q



print("generation keys")
#Generer le paire de clés (dA,QA)

dA = random.randint(1, curve.n-1)

#  generer la paire de cles par la methode DBL et ADD avec inverseEuclid
iter = 1
st1 = time.time()
for i in range(iter):
    QA = mult_double_and_add_EdDSA1(dA,curve.g)
end1 = time.time()
print('DBL_ADD avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode DBL et ADD avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA2 = mult_double_and_add_EdDSA2(dA,curve.g)
end1 = time.time()
print('DBL_ADD avec InverseEuler  :',(end1 - st1)/iter)


#  generer la paire de cles par la methode window avec inverseEuclid
st1 = time.time()
for i in range(iter):
    QA3 = mult_window_EdDSA1(dA,4,curve.g)
end1 = time.time()
print('window avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode window avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA4 = mult_window_EdDSA2(dA,4,curve.g)
end1 = time.time()
print('window avec InverseEuler  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    QA5 = mult_NAF_EdDSA1(dA,curve.g)
end1 = time.time()
print('NAF avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode NAF avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA6 = mult_NAF_EdDSA2(dA,curve.g)
end1 = time.time()
print('NAF avec InverseEuler  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode W-NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    QA7 = mult_w_Naf_EdDSA1(dA,3,curve.g)
end1 = time.time()
print('W-NAF avec InverseEuclid  :',(end1 - st1)/iter)
#  generer la paire de cles par la methode W-NAF  avec inverseEuler
st1 = time.time()
for i in range(iter):
    QA8 = mult_w_Naf_EdDSA2(dA,3,curve.g)
end1 = time.time()
print('W-NAF avec InverseEuler  :',(end1 - st1)/iter)



print("generation signature")

msg="Hello"
msg1="Hel lo"
# Signer le message





#  SM par la methode DBL et ADD avec inverseEuclid
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_double_and_add_EdDSA1(r,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t1 = (end1 - st1)

#  SM par la methode DBL et ADD avec inverseEuler
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_double_and_add_EdDSA2(r,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t2 = (end1 - st1)

#  SM par la methode window avec inverseEuclid
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_window_EdDSA1(r,4,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t3 = (end1 - st1)

#  SM par la methode window avec inverseEuler
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_window_EdDSA2(r,4,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t4 = (end1 - st1)

#  SM par la methode NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_NAF_EdDSA1(r,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t5 = (end1 - st1)

#  SM par la methode NAF avec inverseEuler
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_NAF_EdDSA2(r,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t6 = (end1 - st1)

#  SM par la methode W-NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_w_Naf_EdDSA1(r,3,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t7 = (end1 - st1)

#  SM par la methode W-NAF  avec inverseEuler
st1 = time.time()
for i in range(iter):
    M = str(dA)+msg
    h = int(hashlib.sha256(M.encode()).hexdigest(),16)
    r = h%curve.n
    rpoint = mult_w_Naf_EdDSA2(r,3,curve.g)
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    s = (r+dA*h)%curve.n
end1 = time.time()
t8 = (end1 - st1)





print('DBL_ADD avec InverseEuclid  :',(t1)/iter)
print('DBL_ADD avec InverseEuler  :',(t2)/iter)
print('window avec InverseEuclid  :',(t3)/iter)
print('window avec InverseEuler  :',(t4)/iter)
print('NAF avec InverseEuclid  :',(t5)/iter)
print('NAF avec InverseEuler  :',(t6)/iter)
print('W-NAF avec InverseEuclid  :',(t7)/iter)
print('W-NAF avec InverseEuler  :',(t8)/iter)

print("verification sign")

# Verification de signature



#  SM par la methode DBL et ADD avec inverseEuclid
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_double_and_add_EdDSA1(s,curve.g)
    P2 = point_add_EdDSA1(rpoint, mult_double_and_add_EdDSA1(h,QA))
end1 = time.time()
t1= end1 - st1
#  SM par la methode DBL et ADD avec inverseEuler
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_double_and_add_EdDSA2(s,curve.g)
    P2 = point_add_EdDSA2(rpoint, mult_double_and_add_EdDSA2(h,QA))
end1 = time.time()
t2= end1 - st1
#  SM par la methode window avec inverseEuclid
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_window_EdDSA1(s,4,curve.g)
    P2 = point_add_EdDSA1(rpoint, mult_window_EdDSA1(h,4,QA))
end1 = time.time()
t3 =end1 - st1
#  SM par la methode window avec inverseEuler
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_window_EdDSA2(s,4,curve.g)
    P2 = point_add_EdDSA2(rpoint, mult_window_EdDSA2(h,4,QA))
end1 = time.time()
t4 =end1 - st1
#  SM par la methode NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_NAF_EdDSA1(s,curve.g)
    P2 = point_add_EdDSA1(rpoint, mult_NAF_EdDSA1(h,QA))
end1 = time.time()
t5 = end1 - st1

#  SM par la methode NAF avec inverseEuler
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_NAF_EdDSA2(s,curve.g)
    P2 = point_add_EdDSA2(rpoint, mult_NAF_EdDSA2(h,QA))
end1 = time.time()
t6 = end1 - st1
#  SM par la methode W-NAF avec inverseEuclid
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_w_Naf_EdDSA1(s,3,curve.g)
    P2 = point_add_EdDSA1(rpoint, mult_w_Naf_EdDSA1(h,3,QA))
end1 = time.time()
t7 = end1 - st1

#  SM par la methode W-NAF  avec inverseEuler
st1 = time.time()
for i in range(iter):
    h = int(hashlib.sha256((str(rpoint)+str(QA)+msg).encode()).hexdigest(),16)
    h = h%curve.n
    P1 = mult_w_Naf_EdDSA2(s,3,curve.g)
    P2 = point_add_EdDSA2(rpoint, mult_w_Naf_EdDSA2(h,3,QA))
end1 = time.time()
t8 = end1 - st1


st1 = time.time()
for i in range(iter):
    if (P1==P2):
        sign = "valid"
    else:
        sign = "invalid"
end1 = time.time()
tt = end1 - st1

print('DBL_ADD avec InverseEuclid  :',(t1+tt)/iter)
print('DBL_ADD avec InverseEuler  :',(t2+tt)/iter)
print('window avec InverseEuclid  :',(t3+tt)/iter)
print('window avec InverseEuler  :',(t4+tt)/iter)
print('NAF avec InverseEuclid  :',(t5+tt)/iter)
print('NAF avec InverseEuler  :',(t6+tt)/iter)
print('W-NAF avec InverseEuclid  :',(t7+tt)/iter)
print('W-NAF avec InverseEuler  :',(t8+tt)/iter)
