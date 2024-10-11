"""
This file is written to hold the compression helper functions for the 
main (de)compression functions of FESTA(+).
"""

# Sage imports
from sage.all import EllipticCurve, inverse_mod, ZZ


from Crypto.Hash import SHAKE256
# Local imports
from utilities.fast_sqrt import sqrt_Fp4

# ======================================================== #
#     Helper Functions for compression and decompression   #
# ======================================================== #
def integer_to_bytes(n, byte_len=None):
    """
    Represent an integer n as bytes in big-endian
    """
    n = int(n)
    if byte_len is None:
        byte_len = (n.bit_length() + 7) // 8
    return n.to_bytes(byte_len, 'big')

def bytes_to_integer(b):
    """
    Represent bytes as an integer in big-endian
    """
    return ZZ(int.from_bytes(b, 'big'))

def compress_curve(E, Fp4, p_byte_len ):
    """
    Represent a Montgomery curve as bytes in big-endian
    """
    A = E.a_invariants()[1]
    return integer_to_bytes(A[0], p_byte_len)+ integer_to_bytes(A[2], p_byte_len)

def decompress_curve(Fp4,R,p_byte_len):
    """
    Transform a Montgomery curve given as bytes into the corresponding Elliptic curve sage object
    """

    a_bytes = R[:p_byte_len]
    b_bytes = R[p_byte_len:]
    a_int = bytes_to_integer(a_bytes)
    b_int = bytes_to_integer(b_bytes)
    A = Fp4([a_int,0, b_int,0])
    return EllipticCurve(Fp4, [0, A, 0, 1, 0])


def compress_points(E, zeta, Ps, p_byte_len):
    """
    Represent points Ps on the curve E into bytes in big-endian. 
    They are soly represented by their x-component + 1bit for the y-component.
    """
    
    A = E.a_invariants()[1]
    Rs=[]
    for P in Ps:
        w = P.x()
        x = A+w
        x = (w*x+1)*w
        bool = sqrt_Fp4(x, zeta) == P.y()        
        Rs.append((integer_to_bytes(w[0], p_byte_len)+ integer_to_bytes(w[2], p_byte_len), bool))
    return Rs

def decompress_points(E,Fp4,zeta,Rs,p_byte_len):
    """
    Transform a list of points Ps on the curve E in bytes into the sage object corresponding to 
    theses points.
    """

    A = E.a_invariants()[1]
    Ps=[]
    for R in Rs:
        a_bytes = R[0][:p_byte_len]
        b_bytes = R[0][p_byte_len:]
        a_int = bytes_to_integer(a_bytes)
        b_int = bytes_to_integer(b_bytes)
        X = Fp4([a_int,0, b_int,0])
        d = A+X
        d = (X*d+1)*X

        Y= sqrt_Fp4(d, zeta)
        if not R[1]:
            Y = -Y
        Ps.append(E([X,Y]))
    return Ps 

