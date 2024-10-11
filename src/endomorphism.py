# SageMath imports
from sage.all import (
    EllipticCurve,
    PolynomialRing,
    discrete_log,
    ZZ,
    factor,
    matrix,
    identity_matrix,
    vector
)

#import elliptic_curve as ec
from  utilities.discrete_log import BiDLP 

# action of square root of -1
def i_action(P, zeta2):
    F = P.base_ring()
    E = P.curve()
    assert zeta2 in F
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    X, Y, Z = P
    return E([-X, zeta2*Y, Z])

# Frobenius endomorphism
def Frobenius(P):
    p = P.base_ring().characteristic()
    E = P.curve()
    X, Y, Z = P
    return E([X**p, Y**p, Z**p])

# return retP s.t. 2*retP = P. Note that retP is over an extension field.
def half_point(P, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    assert F.is_subring(F2)

    E2 = EllipticCurve(F2, [1, 0])
    R = PolynomialRing(F2, name="X")
    X = R.gens()[0]
    if P.is_zero():
        f = X**3 + X
    else:
        f = P[2]*(X**2 - 1)**2 - P[0]*4*(X**3 + X)
    xs = f.roots(multiplicities=False)
    assert len(xs) > 0
    x = xs[0]
    y = (x**3 + x).sqrt()
    retP = E2([x, y])
    if not E(2*retP) == P:
        retP = -retP
    assert E(2*retP) == P
    return retP

# the action of (i + j)/2
def i_j_2_action(P, zeta2, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(i_action(halfP, zeta2) + Frobenius(halfP))

# the action of (1 + ij)/2
def one_ij_2_action(P, zeta2, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(halfP + i_action(Frobenius(halfP), zeta2))

# the action of a + bi + c(i + j)/2 + d(1 + ij)/2
def action(alpha, P, zeta2, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    a, b, c, d = alpha
    ret = a*P
    ret += b*i_action(P, zeta2)
    ret += c*i_j_2_action(P, zeta2, F2)
    ret += d*one_ij_2_action(P, zeta2, F2)
    return ret

# matrix of the multiplication by alpha w.r.t. basis of N-torsion subgroup
def action_matrix(alpha, basis, N, zeta2, F2):
    P, Q = basis
    aP, aQ = [action(alpha, R, zeta2, F2) for R in [P, Q]]
    a, c = BiDLP(aP, P, Q, N)
    b, d = BiDLP(aQ, P, Q, N)
    assert aP == a*P + c*Q
    return matrix([[a, b], [c, d]])

# matrices of the multiplications by i, (i + j)/2, (1 + ij)/2
def action_matrices(basis, N, zeta2, F2):
    one = [1,0,0,0]
    Ms = []
    for alpha in [one[-i:]+one[:-i] for i in range(1, 4)]:
        Ms.append(action_matrix(alpha, basis, N, zeta2, F2))
    return Ms

# the action of a + bi + c(i + j)/2 + d(1 + ij)/2 w.r.t. vector representation by fixed basis
def action_by_matrices(alpha, v, action_matrices):
    Ms = [identity_matrix(2)] + action_matrices
    M = sum(alpha[i]*Ms[i] for i in range(4))
    return M * vector(v)

# The image of the basis through i, (i + j)/2, (1 + ij)/2 . (used when we do cannot perform discrete logs over our basis)
def image_matrices(basis, N, zeta2, F2):
    one = [1,0,0,0]
    Ms = []
    for alpha in [one[-i:]+one[:-i] for i in range(1, 4)]:
        Ms.append(image_matrix(alpha, basis, N, zeta2, F2))
    return Ms

# The image of the basis through i, (i + j)/2, (1 + ij)/2 . (used when we do cannot perform discrete logs over our basis)
def image_matrix(alpha, basis, N, zeta2, F2):
    P, Q = basis
    aP, aQ = [action(alpha, R, zeta2, F2) for R in [P, Q]]
    return (aP,aQ)

# the action of a + bi + c(i + j)/2 + d(1 + ij)/2 w.r.t. (used when we do cannot perform discrete logs over our basis)
def action_image(alpha, images_list, basis):
    retP, retQ = alpha[0]*basis[0], alpha[0]*basis[1]
    for i in range(1,4):
        retP += alpha[i]* images_list[i-1][0]        
        retQ += alpha[i]* images_list[i-1][1]
    return (retP, retQ)

# Helper function to construct a quaternion element in O0 from its different component, given by a list alpha
def make_endo(Quat_alg, Basis0, alpha):
    assert len(alpha) == 4
    alpha_endo = Quat_alg(0)
    for j in range(4):
        alpha_endo += Basis0[j]*alpha[j]
    return alpha_endo

# Given gamma a quaternion in O0, compute its respective component in the standard basis (1, i, (i + j)/2, (1 + ij)/2)
def component_endo(gamma):
    return [gamma[0]-gamma[3], gamma[1]-gamma[2],2*gamma[2],2*gamma[3]]

# Cut the components above q. Used to speed-up some computations. 
def endo_reduction(gamma, q):
    return [gamma[i]%q for i in range(len(gamma))]

