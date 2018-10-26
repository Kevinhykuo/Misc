# Algorithum used to calculate pi

import numpy as np
import math as m
from math import factorial as f

# =======
# My way:
# =======

T = 426880 * (10005)**0.5

n = 0
sum = 0
while n<=500:
    t = f(6*n) * (545140134*n+13591409)
    b = f(3*n)*(f(n)**3)*(-262537412640768000)**n
    B = t/b
    sum += B
    n += 1

pi = T/sum

print(str.format('{0:.50f}', pi))

# ==========
# Wikipedia:
# ==========
# from decimal import Decimal as Dec, getcontext as gc
#
# def PI(maxK=70, prec=1008, disp=1007): # parameter defaults chosen to gain 1000+ digits within a few seconds
#     gc().prec = prec
#     K, M, L, X, S = 6, 1, 13591409, 1, 13591409
#     for k in range(1, maxK+1):
#         M = (K**3 - (K<<4)) * M / k**3
#         L += 545140134
#         X *= -262537412640768000
#         S += Dec(M * L) / X
#         K += 12
#     pi = 426880 * Dec(10005).sqrt() / S
#     pi = Dec(str(pi)[:disp]) # drop few digits of precision for accuracy
#     print("PI(maxK=%d iterations, gc().prec=%d, disp=%d digits) =\n%s" % (maxK, prec, disp, pi))
#     return pi
#
# Pi = PI()
# print("\nFor greater precision and more digits (takes a few extra seconds) - Try")
# print("Pi = PI(317,4501,4500)")
# print("Pi = PI(353,5022,5020)")
