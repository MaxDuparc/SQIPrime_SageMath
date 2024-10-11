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


nmb_rep = 50

print_info(" This is the benchmarking of SQIPrime. Each function will be called " + str(nmb_rep) + " times.")



sqiprime117 = SQIPrime(p117, q117, e117, window117, pp117)
sqiprime130 = SQIPrime(p130, q130, e130, window130, pp130)
sqiprime186 = SQIPrime(p186, q186, e186, window186, pp186)
sqiprime240 = SQIPrime(p240, q240, e240, window240, pp240)



SECURITY = ["117","130","186","240"]
SQIPRIME = [sqiprime117, sqiprime130,sqiprime186, sqiprime240]
for i in range(len(SECURITY)):

    t_key = 0
    t_sign = 0
    t_verif = 0
    mes = integer_to_bytes(randint(0,1<<512), 64)

    NAME = "SQIPrime_" + SECURITY[i]
    print_info(f"Running {NAME}")


    for j in range(nmb_rep):
        #KeyGen
        t0 = time.time()
        sk, pk = (SQIPRIME[i]).KeyGen()
        t_key += time.time() - t0
        
        t0 = time.time()
        sign = (SQIPRIME[i]).Sign(sk, pk, mes)
        t_sign += time.time() - t0

        t0 = time.time()
        bool = (SQIPRIME[i]).Verif(pk, mes, sign)
        t_verif += time.time() - t0

        if not bool:
            print_info(f"ERROR: signature invalid")

    print_info(f"Keygen takes on average: {t_key/nmb_rep:.5f} seconds")
    print_info(f"Sign takes on average: {t_sign/nmb_rep:.5f} seconds")
    print_info(f"Verif takes on average: {t_verif/nmb_rep:.5f} seconds")

    print_info("\n\n")





