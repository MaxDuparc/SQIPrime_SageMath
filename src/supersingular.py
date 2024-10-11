"""
Helper functions for the supersingular elliptic curve computations in FESTA
"""

# Sage Imports
from sage.all import *

# Local imports
from utilities.fast_sqrt import sqrt_Fp4
from utilities.discrete_log import weil_pairing_pari


#Compute a random point in E(Fp^4)[l^e].
def point_ord_4(E, l, e):
    #assert is_prime(l)
    p = E.base_ring().characteristic()
    assert (p^2 - 1) % l**e == 0
    
    P = (p^2 - 1)//(l**e)*E.random_point()
    while (l**(e-1)*P).is_zero():
        P = (p^2 - 1)//(l**e)*E.random_point()
    assert (l**e*P).is_zero()
    if P[1][0] >= (p^2 - 1)//2 or (P[1][0] == 0 and P[1][1] >= (p^2 - 1)//2):   # even if the seed of random_point is fixed, the sign of point is random.
        P = -P
    return P

#Compute a basis of E(Fp^4)[l^e].
def basis_4(E, l, e):
    P = point_ord_4(E, l, e)
    Q = point_ord_4(E, l, e)
    if l == 2:
        while (l**(e-1))*(P-Q) == E(0):
            Q = point_ord_4(E, l, e)
    else:
        while weil_pairing_pari((l**(e-1)*P), l**(e-1)*Q, l) == 1:
            Q = point_ord_4(E, l, e)
    return P, Q




#sample a random points in  E(Fp^4)[l^e].
def sample_ran_element(E, Fp4, Fp2, zeta, twist):
    """
    sample a random points in either E or its twist E'.
    Note that theses poitns are represented in Fp4.
    """
    w = Fp2.random_element()
    w = Fp4([w[0],0,w[1],0])
    A = E.a_invariants()[1]
    yw = w**3 + A*w**2+w 
    t = sqrt_Fp4(yw, zeta)
    if twist:
        while t[1] == 0:
            w = Fp2.random_element()
            w = Fp4([w[0],0,w[1],0])
            yw = w**3 + A*w**2+w 
            t = sqrt_Fp4(yw, zeta)
        return E([w,t])
    else:
        while t[1] != 0:
            w = Fp2.random_element()
            w = Fp4([w[0],0,w[1],0])
            yw = w**3 + A*w**2+w 
            t = sqrt_Fp4(yw, zeta)
        return E([w,t])



def Basis_comput(E,Fp4, Fp2, zeta, twist,l, e):
    """
    Compute a basis of E(Fp^4)[l^e]. Note that l^e must be such that 
    l^e| p+1 or p-1.  
    """
    p = Fp4.characteristic()
    P = sample_ran_element(E, Fp4, Fp2, zeta, twist)
    Q = sample_ran_element(E, Fp4, Fp2, zeta, twist)
    if twist:
        P = ((p-1)/(l**e))*P        
        Q = ((p-1)/(l**e))*Q
    else:
        P = ((p+1)/(l**e))*P        
        Q = ((p+1)/(l**e))*Q
    if l == 2:
        A = (l**(e-1))*P
        while A == E(0):
            P = sample_ran_element(E, Fp4, Fp2, zeta, twist)
            P = ((p+1)/(l**e))*P
            A = (l**(e-1))*P
        B = (l**(e-1))*Q
        while B == E(0) or A == B:
            Q = sample_ran_element(E, Fp4, Fp2, zeta, twist)
            Q = ((p+1)/(l**e))*Q
            B = (l**(e-1))*Q
    else:
        while weil_pairing_pari((l**(e-1)*P), l**(e-1)*Q, l) == 1:
            Q = sample_ran_element(E, Fp4, Fp2, zeta, twist)
            if twist:
                Q = ((p-1)/(l**e))*Q
            else:
                Q = ((p+1)/(l**e))*Q

    return P, Q

#Quick function that returns i s.t. order(P) = 2**i. 
def two_order(E,P,e):
    for i in range(e):
        if P == E(0):
            return i
        P = 2*P
    return e