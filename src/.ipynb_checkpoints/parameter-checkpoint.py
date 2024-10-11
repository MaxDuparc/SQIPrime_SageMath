# SageMath imports
from sage.all import *


import src.quaternion as quat
import random as rand
import src.endomorphism as endo

def calcFields(p):
    """
    Compute fields of Fp4, Fp2 and zeta s.t. zeta^2 =-1 all with characteristic p.
    """
    assert p % 4 == 3
    R = PolynomialRing(GF(p), name="x")
    x = R.gens()[0]
    Fp2 = GF(p**2, modulus=x**2+1, name="i")
    z = Fp2.random_element()
    while z.is_square():
        z = Fp2.random_element()
    t = ZZ(z + z**p)
    n = ZZ(z**(p+1))
    R = PolynomialRing(ZZ, name="x")
    x = R.gens()[0]
    Fp4 = GF(p**4, modulus=x**4 - t*x**2 + n, name="z")
    Fp2 = Fp4.subfield(2)
    i = Fp2(-1).sqrt(all=True)[0]
    return Fp4, Fp2, i



#Algorithm 3 of https://eprint.iacr.org/2024/773
def FindPrecomputedBasis(E,Fp4, basis, Qalgebra, Basis0, image_action, zeta, d):
    """
    Algorithm 3 of https://eprint.iacr.org/2024/773
    """
    assert E == EllipticCurve(Fp4, [0, Fp4(0), 0, 1, 0])
    p = E.base_ring().characteristic()
    P,Q = basis[0],basis[1]
    Ideal_0 = Qalgebra.ideal(Basis0)
    assert d*P == E(0) and d*Q == E(0)
    c=1
    while gcd(c,d) != 1 or c * d < p:
        c = rand.randint(p>>64,p)
    c *=  d
    alpha = quat.FullRepresentInteger(c,p)
    a, b = rand.randint(1, d), rand.randint(1, d)
    R,S = E(0), E(0)
    while R == E(0):
        while( a % d == 0 or b % d == 0):
            a, b = rand.randint(1, d), rand.randint(1, d)
        vP = endo.action_image(alpha,image_action, basis)
        R = a*vP[0] + b*vP[1]
    alpha_endo = Qalgebra(0)
    for j in range(4):
        alpha_endo += Basis0[j]*alpha[j]
    Ideal_P = Ideal_0*(alpha_endo.conjugate())+ (Ideal_0*d)
    
    inter_image = endo.image_matrices((R,S), d, zeta, Fp4)
    n = rand.randint(p<<10, p<<32)
    while gcd(n,d) != 1:
        n = rand.randint(p<<10, p<<32)
    iota = quat.FullRepresentInteger(n,p)
    wP = endo.action_image(iota, inter_image ,(R,S))
    S = wP[0] 
    while R.weil_pairing(S, d) == 1:
        while gcd(n,d)  != 1 or n*d<p:
            n = rand.randint(4*p, p*d)
        iota = quat.FullRepresentInteger(n,p)
        wP = endo.action_image(iota, inter_image, vP)
        S  = wP[0] 

    return (R,S), iota, Ideal_P




def power_of_l_dividing(n,l):
    e=0
    while n % l == 0:
        e += 1
        n =n/l
    return l**e




def smooth_part(n,B):
    i = 2
    m = 1
    while i <= B:
        if n % i == 0:
            m = m*i
            n = n/i
        else:
            i = next_prime(i)
    return m



    
def Find_SQIPrime_prime(lamb):
    """
    Function that seach for good SQIprime prime that provide at least lamb bits of security.
    """
    print(f'p = 2*(2^2n*f)-1, q = (2^n - 1)/m1 or q = (2^n + 1)/m2;\n')
    r = lamb 
    for t in range(20):
       x = 2**(r+t)
       for f in range(1,1500,2):
           p = 2*(x*f)**2-1
           if is_prime(p):
               B = (x*f+1).nbits() - lamb
               alpha = power_of_l_dividing(p+1,2)
               m1 = smooth_part(x*f-1,2**B)
               q1 = Integer((x*f-1)/m1)
               if is_prime(q1) and alpha > q1*sqrt(p) and q1.nbits() > lamb:
                   print(f'log p: {p.nbits()},  alpha: {alpha.nbits()}, log q: {q1.nbits()}, n: {r+t}, f:{f}, m1: {m1}')
               m2 = smooth_part(x*f+1,2**B)
               q2 = Integer((x*f+1)/m2)
               if is_prime(q2) and alpha > q2*sqrt(p) and q2.nbits()> lamb:
                   print(f'log p: {p.nbits()},  alpha: {alpha.nbits()}, log q: {q2.nbits()}, n: {r+t}, f:{f}, m2: {m2}')
    
    












