#!/usr/bin/env python

"""
abelian galois extensions / cyclotomic _number fields 
"""

from time import time
start_time = time()
import string

from sage.all_cmdline import *
from sage import all_cmdline 

import numpy

from bruhat.argv import argv
from bruhat.smap import SMap
from huygens.namespace import *

yellow = RGB(1.0, 0.8, 0.2)
orange = RGB(252/255, 132/255, 3/255) * (1./0.8)


def main():
    K = CyclotomicField(12)
    R = K.ring_of_integers()

    z = R.gen(1)

    K_G = CyclotomicField(4)
    G = K_G.ring_of_integers()
    i = G.gen(1)
    assert i**2 == -1

#    #hom = Hom(G, R)
#    #f = hom([1],[1])
#    f = G.hom([1, z**3])
#    print(f)
#    print(f(i))

    N = 30
    ps = []
    for a in range(N):
     for b in range(N):
        n = a+i*b
        if n.is_prime():
            ps.append(n)
        
    def hom(n):
        n = str(n)
        n = n.replace("^", "**")
        n = eval(n, {"zeta4" : z**3})
        return n*R.one()

    def get_coord(n):
        n = eval(str(n), {"zeta4": 1j}) + 0j
        return n.real, n.imag

    cvs = Canvas()
    cvs.stroke(path.line(0, 0, N, 0), st_arrow)
    cvs.stroke(path.line(0, 0, 0, N), st_arrow)
    cvs.clip(path.rect(0, 0, N, N))

    for x in range(N):
     for y in range(N):
        cvs.fill(path.circle(x, y, 0.03))

    #A = numpy.array([[3,2],[-2,3]])
    A = numpy.array([[3,0],[0,3]])
    for i0 in range(-7, N//3+1):
      for i1 in range(-2, N//3+2):
        v0 = numpy.dot(A, [i0,i1])
        vx = numpy.dot(A, [i0+1,i1])
        vy = numpy.dot(A, [i0,i1+1])
        #cvs.fill(path.circle(*v0, 0.1), [orange])
        cvs.stroke(path.line(*v0, *vx), [orange])
        cvs.stroke(path.line(*v0, *vy), [orange])

    for src in ps:
        tgt = hom(src)
        #print(get_coord(src), tgt.is_prime())
        x, y = get_coord(src)
        p = path.circle(x, y, 0.1)
        if tgt.is_prime():
            cvs.fill(p)
        else:
            fs = (tgt.factor())
            #if len(fs)==1:
            #    print(fs, src, len(fs), )
            if src==2+3*i:
                print(src, "-->", tgt, "-->", fs)
            cvs.fill(p, [white])
            cvs.stroke(p, [red])
    #cvs.writePDFfile("images/gauss_split_3.pdf")
    cvs.writePDFfile("images/gauss_split_total.pdf")
    print("images/gauss_split_total.pdf")


def main_34():
    """
    we can think of one of the "combined" quantities as 
    "(a + bi) + (c+di)w", or equivalently as "a + bi + cw + diw"
    ....  and the eisenstein conjugate is "a + bi - cw -diw"
    .  i guess there's also the gaussian conjugate 
    "a - bi + cw - diw" but i don't think we're going to use that
    in this context ....
    """

    cvs = Canvas()
    R = 4.0
    cvs.stroke(path.circle(0, 0, R), [orange])
    cvs.stroke(path.line(-2*R, 0, +2*R, 0), st_arrow)
    cvs.stroke(path.line(0, -2*R, 0, +2*R), st_arrow)

    def roots(N, r=0.1, st=[], fill=False):
        for i in range(N):
            theta = 2*pi*i/N
            x, y = R*cos(theta), R*sin(theta)
            p = path.circle(x, y, r)
            if fill:
                cvs.fill(p, st)
            else:
                cvs.stroke(p, st)

    roots(12, 0.05, fill=True)
    roots(4)
    roots(3, 0.15, [red])

    cvs.writePDFfile("images/cyclotomic_12.pdf")

    


if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


