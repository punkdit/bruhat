#!/usr/bin/env python3

import sys, os
from math import floor

from mpmath import kleinj, mpc # sympy
import numpy
from PIL import Image

from bruhat.argv import argv

#z = mpc(1.0, 0.5)
#print(kleinj(z))

EPSILON = 1e-6

N = argv.get("N", 256)
tol = argv.get("tol", 0.5)

re0 = -1.
im0 = 1e-4

delta = 1./N

r = 3**0.5

A = numpy.zeros((N, N, 3), dtype=numpy.uint8)

for i in range(N):
  for j in range(N):
    #z = mpc(re0 + 2*i*delta, im0 + 2*j*delta)
    z = mpc(2*(1-EPSILON)*i/(N-1)-1., 2*(1-EPSILON)*j/(N-1)-1.)
    z = (z + 1.j)/(1.j*z + 1.)
    if z.imag < EPSILON:
        continue
    z0 = kleinj(z/r)
    z1 = kleinj(z*r)
    x0 = min(abs(z0.imag), tol)
    x1 = min(abs(z1.imag), tol)
    #print(".", end="", flush=True)
    cl0 = int(floor(255 * (1. - x0/tol)))
    cl1 = int(floor(255 * (1. - x1/tol)))
    #if z.real < 0:
    #    cl = (cl, 0, 0)
    #else:
    #    cl = (0, cl, 0)
    cl = (cl0, cl1, 0)
    A[i, j] = cl


    
print()

#print(A)

im = Image.fromarray(A)
im.save("image.png")






