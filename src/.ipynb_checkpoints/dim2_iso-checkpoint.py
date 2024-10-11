from sage import *

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
from utilities.fast_sqrt import sqrt_Fp4


from utilities.discrete_log import weil_pairing_pari

import src.endomorphism as End
import src.quaternion as quat







#Function used to compute the image of a basis through an isogeny Phi.
def eval_basis_unsmooth(E0, E1, Phi, N, Basis, basis_order, position=True):
    """
    Compute the image of a Basis of either E0[basis_order] (position = True) or 
    E1[basis_order] (position = False)  through phi, where phi a dim 1 isogeny 
    of degree N embedded in Phi via Kani's Lemma.
    Phi was computed using theta structures described in https://eprint.iacr.org/2023/1747.
    """
    P = Basis[0]    
    Q = Basis[1]
    PQ = P+Q
    
    if position:
        X = Phi(CouplePoint(P, E1(0)))
        Y = Phi(CouplePoint(Q, E1(0)))
        XY = Phi(CouplePoint(PQ, E1(0)))
        if not (X[0] + Y[0] == XY[0] or X[0] + Y[0] == -XY[0]):
            Y[0] = -Y[0]
        Pd, Qd = X[0], Y[0]
        if not weil_pairing_pari(Pd,Qd, basis_order) == weil_pairing_pari(P,Q, basis_order)**N:
            Good_order = False
            if not (X[1] + Y[1] == XY[1] or X[1] + Y[1] == -XY[1]):
                Y[1] = -Y[1]
            Pd, Qd = X[1], Y[1]
        #assert weil_pairing_pari(Pd,Qd, basis_order) == weil_pairing_pari(P,Q, basis_order)**N
        
    if not position:
        X = Phi(CouplePoint(E0(0),P))
        Y = Phi(CouplePoint(E0(0),Q))
        XY = Phi(CouplePoint(E0(0),PQ))
        if not (X[0] + Y[0] == XY[0] or X[0] + Y[0] == -XY[0]):
            Y[0] = -Y[0]
        Pd, Qd = X[0], Y[0]
        if not weil_pairing_pari(Pd,Qd, basis_order) == weil_pairing_pari(P,Q, basis_order)**N:
            Good_order = False
            if not (X[1] + Y[1] == XY[1] or X[1] + Y[1] == -XY[1]):
                Y[1] = -Y[1]
            Pd, Qd = X[1], Y[1]
        #assert weil_pairing_pari(Pd,Qd, basis_order) == weil_pairing_pari(P,Q, b_order)**N

    return Pd,Qd


def compute_endo(E0, p , q , e , basis, action_matrices):
    """
    Compute the torsion points and prepare them for computation for Kani's Lemma 
    """
    P0,Q0= basis[0],basis[1]

    #Sample the Endomorphisms
    alpha = quat.FullRepresentInteger(q*(2**e-q),p)
    vP = End.action_by_matrices(alpha, [1, 0], action_matrices)
    vP = (vP[0] % 2**e, vP[1] % 2**e)
    alphaP = vP[0]*P0 + vP[1]*Q0
    vQ = End.action_by_matrices(alpha, [0, 1], action_matrices)
    vQ = (vQ[0] % 2**e, vQ[1] % 2**e)

    alphaQ = vQ[0]*P0 + vQ[1]*Q0

    #Compute points
    P1 = (2**e - q)*P0
    Q1 = (2**e - q)*Q0
    P2 = alphaP
    Q2 = alphaQ

    return P1, Q1, P2,Q2, alpha  



def Special_isogeny_1(E0, P1, Q1, P2, Q2, e, Rs):
    """
    Subfunction used to compute the (2,2) isogeny whose kernel is {(R, R) | R in E0[2]}  (E0xE0)-specific case
    """
    P1, P2 = P1 - P2, P1 + P2
    Q1, Q2 = Q1 - Q2, Q1 + Q2
    e = e-1
    for i in range(len(Rs)):
        Rs[i] = (Rs[i][0]-Rs[i][1],Rs[i][0] + Rs[i][1])
    return P1, Q1, P2, Q2, e, Rs
        
def iota_auto(R,E0, zeta):
    return E0([-R[0], zeta*R[1], R[2]])

def Special_isogeny_2(E0, P1, Q1, P2, Q2, e, zeta,  Rs):
    """
    Subfunction used to compute the (2,2) isogeny whose kernel is <((0,0),(0,0)), ((zeta2,0),(-zeta2,0))>. (E0xE0)-specific case
    """

    P1, P2 = (P1 + iota_auto(P2, E0, zeta)), (iota_auto(P1, E0, zeta) + P2)
    Q1, Q2 = (Q1 + iota_auto(Q2, E0, zeta)), (iota_auto(Q1, E0, zeta) + Q2)
    e = e-1
    for i in range(len(Rs)):
        Rs[i] = (Rs[i][0] + iota_auto(Rs[i][1], E0, zeta), iota_auto(Rs[i][0],E0 ,zeta) + Rs[i][1])
    return P1, Q1, P2, Q2, e, Rs


#Special function used to compute (2,2)-0isogeny with domain E0xE0.  
def eval_basis_unsmooth_SpecialDomain(E0, P1 ,Q1 , P2, Q2, e, N, Basis, basis_order ,zeta, position):    
    """
    Special function used to compute and evaluate points through a dim 2 isogeny of degree 2^e when its domain is of the form 
    E0xE0 (j(E0)=1728) (although this code works for any domain ExE except if j(E)==0 (automorphism reasons)).
    P1, Q1, P2, Q2 given by compute_endo. For the others imputs:
    - N is a list of the degree of the dim 1 isogenies that we are evaluating through Phi using Kani's Lemma
    - Basis is a list of basis of E0 that we are evaluating over. Basis[i] = E0[basis_order[i]]
    - basis_order is a list of the order of our respective matrices.
    - position is a list that tells us wether basis is to be taken as a basis of the first (top-left) or second (bottom-right) E0. (See diagramm bellow)
                  d
            E0 -------> E1
            |            ∧    
    2^e-d   |            | 2^e-d   
            |            |    
            ∨      d     |
            E2<---------E0 
    """
    Rs  = []
    Res = []
    E1 = E0
    E2 = E0
    for i in range(len(Basis)):
        P = Basis[i][0]    
        Q = Basis[i][1]
        PQ = P+Q
        if position[i]:
            Rs.append((P , E0(0)))
            Rs.append((Q , E0(0)))
            Rs.append((PQ , E0(0))) 
        else:
            Rs.append((E0(0) , P))
            Rs.append((E0(0) , Q))
            Rs.append((E0(0) , PQ)) 


    #Check bad cases in E0xE0
    while e >=0 :
        
        T1 = 2**(e-1)*P1
        T2 = 2**(e-1)*P2
        S1 = 2**(e-1)*Q1
        S2 = 2**(e-1)*Q2

        if (T1.is_zero() or T2.is_zero() or S1.is_zero() or S2.is_zero()):
            #(2,2) isogeny is diagonal 
            
            if T1.is_zero():
                T1, S1 = S1, T1
                T2, S2 = S2, T2
            assert not T1.is_zero()
            if not S1.is_zero():
                assert S1 == T1
                S2 -= T2
            assert T2.is_zero() or T2 == S2
            phi1 = E1.isogeny(T1)
            phi2 = E2.isogeny(S2)
            e = e-1
            P1, Q1 = phi1(P1), phi1(Q1)
            P2, Q2 = phi2(P2), phi2(Q2)
            E1 = P1.curve()
            E2 = P2.curve()
            for i in range(len(Rs)):
                Rs[i] = (phi1(Rs[i][0]),phi2(Rs[i][1]))
        else:
            if P1.curve() == P2.curve() == E0:
                if S1[0] == 0:
                    T1, S1 = S1, T1
                    T2, S2 = S2, T2
                elif not T1[0] == 0:
                    T1 += S1
                    T2 += S2
                assert T1[0] == 0
                
                if T2[0] == 0:
                    if S1 == S2:
                        #"Special Case n1"
                        P1, Q1, P2, Q2, e, Rs = Special_isogeny_1(E0, P1, Q1, P2, Q2, e, Rs)
                            
            
                    else:
                        #Special Case n2"
                        P1, Q1, P2, Q2, e, Rs = Special_isogeny_2(E0, P1, Q1, P2, Q2, e, zeta,  Rs)
                            
                else:
                    break
            else:
                    break

    
    Pa, Qa = CouplePoint(P1,P2), CouplePoint(Q1, Q2)
    kernel = (Pa, Qa)
    Phi = EllipticProductIsogenySqrt(kernel, e, None, zeta)
    Xs = [Phi(CouplePoint(R[0],R[1])) for R in Rs]
    
    for i in range(0,int(len(Xs)/3)):
         X = Xs[3*i]
         Y = Xs[3*i+1]
         XY = Xs[3*i+2]
         if not (X[0] + Y[0] == XY[0] or X[0] + Y[0] == -XY[0]):
             Y[0] = -Y[0]
         Pd, Qd = X[0], Y[0]
         if not weil_pairing_pari(Pd,Qd, basis_order[i]) == weil_pairing_pari(Basis[i][0],Basis[i][1], basis_order[i])**N[i]:
            if not (X[1] + Y[1] == XY[1] or X[1] + Y[1] == -XY[1]):
                Y[1] = -Y[1]
            Pd, Qd = X[1], Y[1]
            #assert Pd.weil_pairing(Qd, basis_order[i]) == Basis[i][0].weil_pairing(Basis[i][1], basis_order[i])**N[i]
         Res.append((Pd,Qd))
    return Res








