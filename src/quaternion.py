# SageMath imports
from sage.all import (
    ZZ,
    Primes,
    is_pseudoprime,
    Mod,
    sqrt,
    floor,
    randint
)

# Return a, b such that a^2 + b^2 = n if such a, b exist.
def SumOf2Sq(n):
    if n < 0:
        return None
    if n == 0:
        return 0, 0
    if n == 1:
        return 1, 0
    a, b = 1, 0
    for q in Primes()[:100]:
        e = 0
        while n % q == 0:
            n //= q
            e += 1
        a *= q**(e//2)
        b *= q**(e//2)
        if e % 2 == 1:
            if q % 4 == 3:
                return None
            else:
                s, t = cornacchia_smith(q)
                at = a*s - b*t
                bt = a*t + b*s
                a, b = at, bt
    if n % 4 == 1 and is_pseudoprime(n):
        s, t = cornacchia_smith(n)
        at = a*s - b*t
        bt = a*t + b*s
        a, b = at, bt
    elif n > 1:
        return None
    return a, b

# Return a, b such that a^2 + b^2 = q, where q is a prime.
def cornacchia_smith(q):
    if q == 2:
        a, b = 1, 1
    else:
        x = Mod(-1, q).sqrt()
        a, b, c = q, ZZ(x), floor(sqrt(q))
        while b > c:
            a, b = b, a % b
        t = q - b**2
        a, b = b, sqrt(t)

    # randomize output
    r = randint(0, 1)
    if r == 0:
        a, b = b, a
    r1 = randint(0, 1)
    r2 = randint(0, 1)
    return (-1)**r1*a, (-1)**r2*b

# return a, b, c, d s.t. the norm of a + bi + c(i+j)/2 + d(1 + ij)/2 is N.
def FullRepresentInteger(N, p):
    upper_bound = floor(4*N/p)
    assert upper_bound > 0 # N > p
    while True:
        z = randint(0, floor(sqrt(upper_bound)))
        t = randint(0, floor(sqrt(upper_bound - z**2)))
        M = 4*N - p*(z**2 + t**2)
        tmp = SumOf2Sq(M)
        if tmp:
            x, y = tmp
            if (x - t) % 2 == 0 and (y - z) % 2 == 0:
                return [ZZ(v) for v in [(x-t)/2, (y-z)/2, z, t]]

# norm of a + bi + c(i+j)/2 + d(1 + ij)/2
def norm(alpha, p):
    a, b, c, d = alpha
    return (a + d/2)**2 + (b + c/2)**2 + p*(c**2 + d**2)/4

# a - b - c(i+j)/2 + d(1 - ij)/2
def involution(alpha):
    a, b, c, d = alpha
    return [a + d, -b, -c, -d] 
















