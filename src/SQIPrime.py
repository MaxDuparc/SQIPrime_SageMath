
from sage.all import *

from sage.modules.free_module_integer import IntegerLattice





import src.quaternion as quat
import src.endomorphism as End
import random as rand
import src.compression as comp
import sympy
import src.supersingular as supersingular
import src.dim2_iso as dim2_iso
import src.parameter as param

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
from utilities.discrete_log import BiDLP_power_two
from utilities.fast_sqrt import sqrt_Fp4


from Crypto.Hash import SHAKE256

def arity(n):
    return (ZZ(int(n) & int(~(int(n - 1)))).nbits())-1





#------------------------------------------------------------------------------
                                #SQIPrime
#------------------------------------------------------------------------------

class SQIPrime_Sigma :
    def __init__(self, p, q, e, window = [], pp = None):

        self.p = p
        self.q = q
        self.e = e
        self.window = window
        
        assert (p+1) % 2**e == 0
        assert (p-1) % q == 0

        if pp == None:
            #Compute the standard fields used in SQIPrime
            Fp4, Fp2, zeta = param.calcFields(p)
        
            self.Fp4 = Fp4
            self.Fp2 = Fp2
            self.zeta = zeta 
            self.Fq = GF(q) 

            self.E0 = EllipticCurve(Fp4, [0, 0, 0, 1, 0])

            self.Basis_E0_2 = supersingular.Basis_comput(self.E0, Fp4, Fp2, zeta, False,2, e)
            Basis_E0_q = supersingular.Basis_comput(self.E0, Fp4, Fp2, zeta, True,q, 1)
            self.action_matrices_2 = End.action_matrices(self.Basis_E0_2, 2**e, zeta, Fp4)
        
            image_action_q =  End.image_matrices(Basis_E0_q, q, zeta, Fp4)
            
            #Compute the Precomputed Basis.
            self.Quat_alg = QuaternionAlgebra(-1, -p)
            self.Order0_basis = (self.Quat_alg((1,0,0,0)),self.Quat_alg((0,1,0,0)),self.Quat_alg((0,1/2,1/2,0)),self.Quat_alg((1/2,0,0,1/2)))

            Basis_q, iota, Ideal_P = param.FindPrecomputedBasis(self.E0, Fp4, Basis_E0_q, self.Quat_alg, self.Order0_basis, image_action_q, zeta, q)
            self.Basis_E0_q = Basis_q
            self.Ideal_P = Ideal_P
            self.iota_end = End.make_endo(self.Quat_alg, self.Order0_basis, iota)

            #Check vailidy of the Precomputed basis
            image_action_q = End.image_matrices(Basis_q, q, zeta, Fp4)
            Im = End.action_image(iota, image_action_q, Basis_q)
        
            assert Im[0]==Basis_q[1]
            assert q*Basis_q[0] == self.E0(0) and q*Basis_q[0] == self.E0(0)

        else:
            #load Fields info
            t,n = pp[0]
            
            R = PolynomialRing(ZZ, name="x")
            x = R.gens()[0]
            self.Fp4 = GF(p**4, modulus=x**4 - t*x**2 + n, name="z")
            self.Fp2 = self.Fp4.subfield(2)
            squareroot = self.Fp2(-1).sqrt(all = True)
            if ZZ(squareroot[0][0]) %2 == 0:
                self.zeta = squareroot[1]
            else:
                self.zeta = squareroot[0]
            self.Fq = GF(q)

            
            self.E0 = EllipticCurve(self.Fp4, [0, 0, 0, 1, 0])

            #load Basis_E0_2

            X1= self.Fp4([pp[1][0][0][0],0, pp[1][0][0][1], 0])
            Y1 = self.Fp4([pp[1][0][1][0],0, pp[1][0][1][1], 0])
            X2= self.Fp4([pp[1][1][0][0],0, pp[1][1][0][1], 0])
            Y2 = self.Fp4([pp[1][1][1][0],0, pp[1][1][1][1], 0])
            
            self.Basis_E0_2 = self.E0([X1,Y1]), self.E0([X2,Y2])

            #load action_matrices_2 
            self.action_matrices_2 = [matrix(pp[2][0]),matrix(pp[2][1]),matrix(pp[2][2])] 
            
            self.Quat_alg = QuaternionAlgebra(-1, -p)
            self.Order0_basis = (self.Quat_alg((1,0,0,0)),self.Quat_alg((0,1,0,0)),self.Quat_alg((0,1/2,1/2,0)),self.Quat_alg((1/2,0,0,1/2)))

            #load Basis_E0_q

            X3= self.Fp4([pp[3][0][0][0],0, pp[3][0][0][1], 0])
            Y3 = self.Fp4([0,pp[3][0][1][0],0, pp[3][0][1][1]])
            X4= self.Fp4([pp[3][1][0][0],0, pp[3][1][0][1], 0])
            Y4 = self.Fp4([0,pp[3][1][1][0],0, pp[3][1][1][1]])

            self.Basis_E0_q = self.E0([X3,Y3]), self.E0([X4,Y4])

            #load Ideal_P 

            self.Ideal_P = self.Quat_alg.ideal(pp[4])

            #load iota 

            self.iota_end = End.make_endo(self.Quat_alg, self.Order0_basis, pp[5])
            
            #Verify the precomputed basis  
            image_action_q = End.image_matrices(self.Basis_E0_q, self.q, self.zeta, self.Fp4)
            Im = End.action_image(pp[5], image_action_q, self.Basis_E0_q)
        
            assert Im[0]==self.Basis_E0_q[1]
            assert q*self.Basis_E0_q[0] == self.E0(0) and q*self.Basis_E0_q[0] == self.E0(0)

    def __repr__(self):
        return f"SQIPrime instance over the prime p={self.p} and q={self.q}"
    
#------------------------------------------------------------------------------
                                #KEYGEN
#------------------------------------------------------------------------------


    def KeyGen(self):

        
        # Compute an endomorphism alpha of size q(2^q.nbits()-q)(2^e -q(2^q.nbits()-q))
        P1, Q1, P2,Q2, alpha = dim2_iso.compute_endo(self.E0, self.p , self.q*(2**self.q.nbits()-self.q) , self.e , self.Basis_E0_2, self.action_matrices_2)

        # Compute a basis of E0[2**q.nbits()]  
        P0_, Q0_ = 2**(self.e-self.q.nbits())*self.Basis_E0_2[0],2**(self.e-self.q.nbits())*self.Basis_E0_2[1]

        # Compute the ideal corresponding to rho. 
        Ideal_0 = self.Quat_alg.ideal(self.Order0_basis)
        alpha_endo = End.make_endo(self.Quat_alg, self.Order0_basis, alpha)
        I_rho=(Ideal_0*conjugate(alpha_endo)) + (Ideal_0)*((2**self.q.nbits()-self.q)*(2**self.e-self.q*(2**self.q.nbits()-self.q)))
    

        size1= self.q*(2**self.q.nbits() - self.q)
        size2 = 2**self.e-self.q*(2**self.q.nbits()-self.q)

        # Evaluate via the first Kani's lemma of ExtendedKaniDoublePath the following basis:
        #E0[2^qnbits()] through hat(rho_2) ◯ tau, 
        # E0[q] through rho_1
        # E0[2^e] through rho_1.
        
        Res = dim2_iso.eval_basis_unsmooth_SpecialDomain(self.E0, P1,Q1,P2,Q2, self.e, 
                                                [size1, size2,size2], 
                                                [(P0_,Q0_), self.Basis_E0_q,self.Basis_E0_2], [2**self.q.nbits(), self.q, 2**self.e],self.zeta, [True,False, False]
                                               )
        
        Image_Basis_rho1_q = Res[1]
        Image_Basis_rho1_2 = Res[2]
        P2_, Q2_ = Res[0][0],Res[0][1]
        
        P1_ = (-self.q)*P0_
        Q1_ = (-self.q)*Q0_
        
    
        #Construct the second Kani, that correspond to an isogeny of degree 2^q.nbits()
        Pa, Qa = CouplePoint(P1_,P2_), CouplePoint(Q1_, Q2_)
        kernel = (Pa, Qa)
        Phi_A2 = EllipticProductIsogenySqrt(kernel, self.q.nbits(), None, self.zeta)

        E_A1 = Phi_A2.domain()[1]

        #Evaluate through tau E0[2^e]
        Image_Basis_tau_E_key_2 = dim2_iso.eval_basis_unsmooth(self.E0, E_A1, Phi_A2, self.q, self.Basis_E0_2, 2**self.e, True)
        E_key = Image_Basis_tau_E_key_2[0].curve()
    
        #Evaluate through rho E0[2^e] and E0[q]
        Image_Basis_rho_E_key_q = dim2_iso.eval_basis_unsmooth(self.E0, E_A1, Phi_A2, 2**self.q.nbits()-self.q, Image_Basis_rho1_q,self.q, False)
        Image_Basis_rho_E_key_2 = dim2_iso.eval_basis_unsmooth(self.E0, E_A1, Phi_A2, 2**self.q.nbits()-self.q, Image_Basis_rho1_2,2**self.e, False)
    
        #assert Image_Basis_rho_E_key_q[0].curve() == E_key
        #assert Image_Basis_rho1_q[0].curve() == E_A1

        #Mask our basis of E_key[q]
        Mask_matrix = random_matrix(self.Fq,2,2) 
        R = (Mask_matrix[0][0]*Image_Basis_rho_E_key_q[0]+Mask_matrix[0][1]*Image_Basis_rho_E_key_q[1])
        S = (Mask_matrix[1][0]*Image_Basis_rho_E_key_q[0]+Mask_matrix[1][1]*Image_Basis_rho_E_key_q[1])

        # Ccnstruct our respective piblic_key, secret_key.
        public_key = E_key,(R,S)
        secret_key = I_rho, Image_Basis_tau_E_key_2, Image_Basis_rho_E_key_2, Mask_matrix
        return secret_key, public_key

#------------------------------------------------------------------------------
                                #COMMIT
#------------------------------------------------------------------------------

    def Commit(self):

        #Find a prime degree between [0,p^(2/3)].
        deg_1 = 1
        while not sympy.isprime(deg_1:=rand.randint(0,2**int(2*self.e/3))):
            pass

        # Compute an endomorphism alpha of size env. deg_1(2**e - deg_1)
 
        P1, Q1, P2,Q2, alpha = dim2_iso.compute_endo(self.E0, self.p , deg_1 , self.e , self.Basis_E0_2 , self.action_matrices_2)

        # Compute the ideal corresponding to psi. 
        Ideal_0 = self.Quat_alg.ideal(self.Order0_basis)
        alpha_endo = End.make_endo(self.Quat_alg, self.Order0_basis, alpha)
        I_psi=(Ideal_0*(alpha_endo)) + (Ideal_0)*(deg_1)

        # Evaluate E0[2**e] through psi.
        Image_Basis_E_com_2 =  dim2_iso.eval_basis_unsmooth_SpecialDomain(self.E0, P1 ,Q1 , P2, Q2,self.e, [deg_1], [self.Basis_E0_2], [2**self.e],self.zeta, [True])[0]  
        E_com = Image_Basis_E_com_2[0].curve()
    
            
        commit = E_com
        secret_commit = I_psi, Image_Basis_E_com_2, deg_1
        return secret_commit, commit

#------------------------------------------------------------------------------
                                #RESPONSE
#------------------------------------------------------------------------------

    def Response(self, secret_key, public_key, secret_commit, commit, chall):
        #Loading
        E_key,(R,S) = public_key 
        I_rho, Image_Basis_tau_E_key_2, Image_Basis_rho_E_key_2, Mask_matrix = secret_key
        E_com = commit
        I_psi, Image_Basis_E_com_2, deg_1 = secret_commit
        ker_chall = R + chall*S 

        
        deg_rho = (2**self.q.nbits()-self.q)*(2**self.e-self.q*(2**self.q.nbits()-self.q))
        
        #Compute the challenge ideal
        matrixbc = Mask_matrix.transpose() * Matrix([[self.Fq(1)],[self.Fq(chall)]])
        Precomputed_basis_endo = self.Quat_alg(matrixbc[0][0])+self.Quat_alg(matrixbc[1][0])*self.iota_end
    
        I_rhoI_varphi = (self.Ideal_P*(Precomputed_basis_endo**(-1))).intersection(I_rho)
    
        # Compute the answer ideal
        J = I_psi.conjugate()*I_rhoI_varphi
        norm_J = self.q*deg_1*deg_rho
        beta = J.minimal_element()
        d = beta.reduced_norm()/norm_J

        # Check beta is small-enough. 
        assert d*self.q < 2**self.e 

        # Check the parity of our number.
        e2 = arity(d)
        EVEN = (e2 == 0)
        e_auxiliary = self.e - e2
        d_auxiliary = d>>e2

        I_sigma = J * (beta.conjugate() / norm_J)
        I_psiI_sigma = I_psi*I_sigma 
    
        
        #Compute corresponding endomorphism
        gamma = (I_psiI_sigma.conjugate()*I_rhoI_varphi).minimal_element()    
        gamma_list = End.component_endo(gamma)


    
        # Compute images of kappa using EvalTorsion.
        vP = End.action_by_matrices(gamma_list, [1,0], self.action_matrices_2)
        vQ = End.action_by_matrices(gamma_list, [0,1], self.action_matrices_2)
    
        X3 = vP[0]*Image_Basis_E_com_2[0] + vP[1]*Image_Basis_E_com_2[1]
        Y3 = vQ[0]*Image_Basis_E_com_2[0] + vQ[1]*Image_Basis_E_com_2[1]
        
        Z2e = Integers(2**self.e)
        lambda_ = Z2e(deg_1)**(-1)
        
        X4, Y4 =  lambda_*X3, lambda_*Y3
    
        
    
        
        
    
        # Construct auxiliary endomorphism of degree d_auxiliary*(2**e_auxiliary-(self.q*d_auxiliary)
        
        sigma = quat.FullRepresentInteger(d_auxiliary*(2**e_auxiliary-(self.q*d_auxiliary)),self.p)
        vP = End.action_by_matrices(sigma, [1, 0], self.action_matrices_2)
        vP = ((vP[0])%2**e_auxiliary, (vP[1])%2**e_auxiliary)    
        
        vQ = End.action_by_matrices(sigma, [0, 1], self.action_matrices_2)
        vQ = ((vQ[0]) % 2**e_auxiliary, (vQ[1]) % 2**e_auxiliary)
        
        sigmaP = (2**e2*vP[0])*self.Basis_E0_2[0] + (2**e2*vP[1])*self.Basis_E0_2[1]
        sigmaQ = (2**e2*vQ[0])*self.Basis_E0_2[0] + (2**e2*vQ[1])*self.Basis_E0_2[1]
    
        
        P2 = -(self.q*(d_auxiliary**2)*2**e2)*Image_Basis_tau_E_key_2[0]
        Q2 = -(self.q*(d_auxiliary**2)*2**e2)*Image_Basis_tau_E_key_2[1]
    
        P1 = (d_auxiliary*self.q)*sigmaP
        Q1 = (d_auxiliary*self.q)*sigmaQ
    
        # Compute delta the auxiliary isogeny of degree 2^e_auxiliary -q*d_auxiliary.

        Pa, Qa = CouplePoint(P1,P2), CouplePoint(Q1, Q2)
        kernel = (Pa, Qa)
        Phi_aux = dim2_iso.EllipticProductIsogenySqrt(kernel, e_auxiliary, None, self.zeta)
    
        
         
        #Evaluate points
        Image_Basis_E_aux_2 = dim2_iso.eval_basis_unsmooth(self.E0, E_key, Phi_aux, 2**e_auxiliary-self.q*d_auxiliary, Image_Basis_rho_E_key_2 , 2**self.e, False)
        V = Phi_aux(CouplePoint(self.E0(0), ker_chall))
        E_aux = Image_Basis_E_aux_2[0].curve()
        
        if V[0].curve() == E_aux:
            V = V[0]
        else:
            V = V[1]
        
        
        #Sample new basis of E_aux[2**e]
        Basis_E_aux_2 = supersingular.Basis_comput(E_aux, self.Fp4, self.Fp2, self.zeta, False ,2, self.e)
    
    
        #Compute discrete log for the dual of delta. 
        disc_log_image=[]
        disc_log_image.append(BiDLP_power_two(Basis_E_aux_2[0],Image_Basis_E_aux_2[0],Image_Basis_E_aux_2[1], self.e, self.window))
        disc_log_image.append(BiDLP_power_two(Basis_E_aux_2[1],Image_Basis_E_aux_2[0],Image_Basis_E_aux_2[1], self.e, self.window))
    

        #Compute T, U = [(2^e_aux -qd)^(-1) mod 2^e_aux]kappa hat(delta)(Basis_E_aux_2)
        
        T = -(disc_log_image[0][0]*X4 +  disc_log_image[0][1]*Y4)
        U = -(disc_log_image[1][0]*X4 +  disc_log_image[1][1]*Y4)
        

        signature = E_com, E_aux, Basis_E_aux_2, T,U,V, e2
        return signature


#-------------------------------------------------------------------------------------------------------
                                #VERIFY
#-------------------------------------------------------------------------------------------------------

    
    def Verify(self, public_key, chall, signature):
        
        #Load 
        E_key ,(R,S) = public_key
        E_com, E_aux, Basis_E_aux_2, T,U,V, e2 = signature
        ker_chall = R + chall*S 
        pos = 0
        
        e_auxiliary = self.e - e2
        E_com_auxiliary = E_com
        
        #Check points are in the rightr curve
        if not (V.curve() ==E_aux or (U.curve() == T.curve() == E_com)):
            return False
    
        #Even Case
        if e2 !=0:
            
            T_aux = 2**e_auxiliary*T
            U_aux = 2**e_auxiliary*U
    
            if E_com(0) != 2**(e2-1)*T_aux:
                kernel_sigma2 = T_aux
            else:
                kernel_sigma2 = U_aux
    
            sigma2 = E_com.isogeny(kernel_sigma2)
            E_com_auxiliary = sigma2.codomain()
            T_inter, U_inter = sigma2(T), sigma2(U)
    
            Pa, Qa = CouplePoint((2**e2)*Basis_E_aux_2[0],T_inter), CouplePoint((2**e2)*Basis_E_aux_2[1], U_inter)
      
        else:
            #Odd case
            Pa, Qa = CouplePoint(Basis_E_aux_2[0],T), CouplePoint(Basis_E_aux_2[1], U)

        #Compute the dim-2 isogeny used for verification.
        kernel = (Pa, Qa)
        Phi_verif = EllipticProductIsogenySqrt(kernel, e_auxiliary, None, self.zeta)
    
        #Check codomain
        if Phi_verif.codomain()[0].j_invariant() == E_key.j_invariant():
            pos = 0
        elif Phi_verif.codomain()[1].j_invariant() == E_key.j_invariant():
            pos = 1
        else:
            return False
    
        
        #Compute point
        Inter_point = Phi_verif(CouplePoint(V,E_com_auxiliary(0)))
        bool1 = (Phi_verif.codomain()[pos]).isomorphism_to(E_key)(Inter_point[pos]).x() == ((2**e_auxiliary) * (ker_chall)).x()
        bool2 = Inter_point[1-pos] == Phi_verif.codomain()[1-pos](0)
    
        #Compute second point to counter falsification.
        W = supersingular.sample_ran_element(E_aux, self.Fp4, self.Fp2, self.zeta, True)
        W = ((self.p-1)//self.q)*W
        #W = point_ord_4(E_delta, q,1)
        Inter_point2 = Phi_verif(CouplePoint(W,E_com_auxiliary(0)))
        bool3 = (Inter_point2[1-pos] != Phi_verif.codomain()[1-pos](0))
        return bool1 and bool2 and bool3
        
        
    
    
    


class SQIPrime(SQIPrime_Sigma):
        def __init__(self, p, q, e, window=[], pp = None):
            super().__init__(p, q, e, window, pp)
            self.byte_size = (p.nbits()+7)//8
        
        def Hash(self, m):
            shake = SHAKE256.new(m)
            return shake.read(self.byte_size) 
        
        def KeyGen(self):
            sk, pk = super().KeyGen()
            #if byte_rep:
            #    pk_byte
            #    comp.compress_curve(pk[0], self.byte_size) 
            #    comp.compress_points(pk[0],self.zeta, pk[1],self.byte_size )
            return sk, pk
            
        def Sign(self, sk, pk, message):

            sec, com = self.Commit()
            chall = comp.bytes_to_integer(self.Hash(message)) % self.q #Imperfect but this is only a PoC.
            sign = self.Response(sk, pk, sec, com, chall)

            return sign

        def Verif(self, pk, message, sign):

            chall = comp.bytes_to_integer(self.Hash(message)) % self.q #Imperfect but this is only a PoC.
            
            return self.Verify(pk, chall,sign)
            












