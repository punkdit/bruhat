#!/usr/bin/env python3

import sys, os
from math import floor

from mpmath import kleinj, mpc
import numpy
from PIL import Image

from bruhat.argv import argv

#z = mpc(1.0, 0.5)
#print(kleinj(z))

N = argv.get("N", 256)
tol = argv.get("tol", 0.5)

re0 = -1.
im0 = 1e-4

delta = 1./N


A = numpy.zeros((N, N, 3), dtype=numpy.uint8)

for i in range(N):
  for j in range(N):
    z = kleinj(mpc(re0 + 2*i*delta, im0 + 2*j*delta))
    x = abs(z.imag)
    if x > tol:
        continue
    #print(".", end="", flush=True)
    cl = int(floor(255 * (1. - x/tol)))
    if z.real < 0:
        cl = (cl, 0, 0)
    else:
        cl = (0, cl, 0)
    A[N-j-1, i] = cl
    
print()

#print(A)

im = Image.fromarray(A)
im.save("image.png")






