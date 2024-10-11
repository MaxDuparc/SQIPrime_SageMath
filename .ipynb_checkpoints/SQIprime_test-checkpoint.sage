#!/usr/bin/env sage

# Python imports
import sys
import time
from random import randint 

# Local imports
from src.SQIPrime import SQIPrime
from utilities.utils import print_info, speed_up_sagemath
from src.compression import integer_to_bytes
from precomputed_values.precomputed_public_parameters import *

# Sage goes vroom!
# Sets up a monkey-patch to cache vector_space which
# dramatically helps with performance of the genus-two
# arithmetic and isogenies
speed_up_sagemath()



print_info("This is a proof of concept build of SQIPrime using SAGE-math. Note that it uses curves defined over Fp^4 instead of the desired Fp^2.")

#Setting the randomness


sqiprime117 = SQIPrime(p117, q117, e117, window117, pp117)
sqiprime130 = SQIPrime(p130, q130, e130, window130, pp130)
sqiprime186 = SQIPrime(p186, q186, e186, window186, pp186)
sqiprime240 = SQIPrime(p240, q240, e240, window240, pp240)



SECURITY = ["117","130","186","240"]
SQIPRIME = [sqiprime117, sqiprime130,sqiprime186, sqiprime240]

for i in range(len(SECURITY)):

    NAME = "SQIPrime_" + SECURITY[i]
    print_info(f"Running {NAME}")
    
    # ============================== #
    #              Keygen            #
    # ============================== #
    t0 = time.time()
    sk, pk = (SQIPRIME[i]).KeyGen()
    keygen_time = time.time() - t0
    print_info(f"Keygen took: {keygen_time:.3f} seconds")
    
    # ============================== #
    #           Sign                 #
    # ============================== #
    mes = integer_to_bytes(randint(0,1<<512), 64)
    t0 = time.time()
    sign = (SQIPRIME[i]).Sign(sk, pk, mes)
    sign_time = time.time() - t0
    print_info(f"Sign took: {sign_time:.3f} seconds")
    
    # ============================== #
    #           Verify               #
    # ============================== #
    t0 = time.time()
    bool = (SQIPRIME[i]).Verif(pk, mes, sign)
    Verify_time = time.time() - t0
    print_info(f"Verify took: {Verify_time:.3f} seconds")
    print_info(f"Is the signature valid: {bool}")



    print_info("\n\n\n")




