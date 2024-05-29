#!/usr/bin/env python

from cmath import *

#from mpmath import kleinj, mpc # sympy
import numpy
from PIL import Image

from numba import jit

from bruhat.argv import argv

EPSILON = 1e-6


@jit
def weirstrass(z, tau=1j, N=5):
    z = complex(z)
    u = 1 / (z**2)
    for i in range(-N, N+1):
      for j in range(-N, N+1):
        if i==j==0:
            continue
        lmda = i + tau*j
        #if abs(z-lmda) < EPSILON:
        #    assert 0, cstr(z)
        #assert abs(lmda) > EPSILON, (i,j)
        u += 1 / ((z-lmda)**2) - 1/(lmda**2)
    return u


def cstr(z):
    s = "%.4f+%.4f"%(z.real, z.imag)
    s = s.replace("+-", "-")
    if s.startswith("-"):
        pass
    else:
        s = " "+s
    return s


def render(N, tau):
    A = numpy.zeros((N, N), dtype=complex)
    for i in range(1, N):
      for j in range(1, N):
        z = i + j*1j + (0.098 + 0.12j)
        z = 4*z/N
        u = weirstrass(z, tau, N=5)
        #u1 = weirstrass(z, N=10)
        #print(cstr(z), cstr(u), cstr(u1), cstr(weirstrass(z,20)))
        A[i, j] = u
      #print('.', end='', flush=True)
    print()

    A0 = A.real
    A0 -= A0.min() - 0.1
    #A0 /= A0.max()

    A1 = A.imag
    A1 -= A1.min() - 0.1
    #A1[A1>100] = 100
    #A1 /= A1.max()

    #print()
    #print(A0.min(), A0.max())
    #print(A1.min(), A1.max())

    A = numpy.zeros((N, N, 3), dtype=numpy.uint8)
    A[:,:,0] = A0
    A[:,:,1] = A1
    A[:,:,2] = 0.2*(A1 + A0)
    im = Image.fromarray(A)
    return im


def main():

    theta = 2*pi/3
    tau = cos(theta) + 1j*sin(theta)
    N = 128*8

    im = render(N, 1j)
    im.save("weierstrass_gauss.png")

    #im = render(N, tau)
    #im.save("weierstrass_eisenstein.png")




if __name__ == "__main__":

    from time import sleep, time
    start_time = time()
    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        main()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))



