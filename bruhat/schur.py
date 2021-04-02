#!/usr/bin/env python3

"""
Earlier version: qupy.ldpc.cell
Used by: bruhat.morse
See also: bruhat.vec 
"""

from random import shuffle, seed
from math import sin, cos, pi

import numpy

from bruhat.argv import argv
from bruhat import element
from bruhat import elim
from bruhat.elim import (
    zeros, rand, dot, identity, eq, coequalizer, compose,
    rank, pseudo_inverse)
#from bruhat.solve import parse
from bruhat import solve

def parse(ring, decl):
    A = solve.parse(decl)
    B = zeros(ring, *A.shape)
    B = B+A
    return B


def shortstr(A):
    s = elim.shortstr(A)
    s = str(s)
    s = s.replace(" 0 ", " . ")
    return s


def kron(ring, A, B):
    if 0 in A.shape or 0 in B.shape:
        C = zeros(ring, A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])
    else:
        #print("kron", A.shape, B.shape)
        C = numpy.kron(A, B)
        #print("\t", C.shape)
    return C


def sum_swap(ring, a, b):
    "swap isomorphism a+b -> b+a"
    one = ring.one
    g = zeros(ring, a+b, a+b)
    for i in range(a):
        g[i, i+a] = one
    for i in range(b):
        g[i+b, i] = one
    return g

def direct_sum(ring, f, g):
    mf, nf = f.shape
    mg, ng = g.shape
    h = zeros(ring, mf+mg, nf+ng)
    h[:mf, :nf] = f
    h[mf:, nf:] = g
    return h


def tensor_swap(ring, m1, m2):
    A = zeros(ring, m2*m1, m1*m2)
    one = ring.one
    for i in range(m1):
      for j in range(m2):
        A[j*m1 + i , i*m2 + j] = one
    return A


def test():

    p = argv.get("p", 2)
    ring = element.FiniteField(p)

    s = tensor_swap(ring, 3, 4)
    si = tensor_swap(ring, 4, 3)
    #print(shortstr(s))
    assert eq(dot(ring, si, s), identity(ring, 3*4))
    assert eq(dot(ring, s, si), identity(ring, 4*3))

    m, n = 2, 3
    A1 = rand(ring, m, n, 1, 1)
    A2 = rand(ring, m, n, 1, 1)

    B = kron(ring, A1, A2)

    for m in range(1, 5):
        I = identity(ring, m*m)
        s = tensor_swap(ring, m, m)
        f = coequalizer(ring, I, s)
    
        assert eq(compose(ring, s, f), f)
        assert rank(ring, f) == [1, 3, 6, 10][m-1]

    # ---------------------------------

    m = argv.get("m", 3)
    n = argv.get("n", 4)
    H = rand(ring, m, n, 1, 1)
    #print("H:")
    #print(shortstr(H))

    if argv.toric:
        H = zeros(ring, m, m)
        for i in range(m):
            H[i, i] = ring.one
            H[i, (i+1)%m] = ring.one
        print("H:")
        print(shortstr(H))

    hypergraph_product(ring, H, H)


def schur(ring, m):
    I = identity(ring, m*m)
    s = tensor_swap(ring, m, m)
    f = coequalizer(ring, I, s)
    return f


def is_zero(ring, A):
    return eq(A, zeros(ring, *A.shape))


def hypergraph_product(ring, A, check=False):
    #print("hypergraph_product: A=%s, B=%s"%(A.shape, B.shape))

    n, m = A.shape

    In = identity(ring, n)
    Im = identity(ring, m)

    Hzs = kron(ring, Im, A), -kron(ring, A, Im)
    Hz = numpy.concatenate(Hzs, axis=0) # horizontal concatenate

    Hxs = kron(ring, A, In), kron(ring, In, A)
    #print("Hxs:", Hxs[0].shape, Hxs[1].shape)
    Hx = numpy.concatenate(Hxs, axis=1) # horizontal concatenate

    assert is_zero(ring, dot(ring, Hx, Hz))

    H1 = Hz
    H0 = Hx
    #print(shortstr(dot(ring, H0, H1)))
    
    assert H1.shape == (n*m+m*n, m*m)
    assert H0.shape == (n*n, n*m+m*n)

    f2 = schur(ring, m)
    f0 = schur(ring, n)

    a = direct_sum(ring,
        tensor_swap(ring, m, n),
        tensor_swap(ring, n, m))
    b = -sum_swap(ring, n*m, m*n)

    I = identity(ring, m*n+n*m)
    f1 = coequalizer(ring, I, compose(ring, a, b))

    #print(shortstr(f1))

    f1i = pseudo_inverse(ring, f1)
    assert(eq(compose(ring, f1i, f1), identity(ring, f1i.shape[1])))
    f2i = pseudo_inverse(ring, f2)
    assert(eq(compose(ring, f2i, f2), identity(ring, f2i.shape[1])))

    J1 = compose(ring, f2i, H1, f1)
    J0 = compose(ring, f1i, H0, f0)

    #print(shortstr(compose(ring, J1, J0)))
    assert is_zero(ring, compose(ring, J1, J0))

    print("H1:", H1.shape)
    print(shortstr(H1))
    print("H0:", H0.shape)
    print(shortstr(H0))

    print("J1:", J1.shape)
    print(shortstr(J1))
    print("J0:", J0.shape)
    print(shortstr(J0))

    assert eq(compose(ring, H1, f1), compose(ring, f2, J1))
    assert eq(compose(ring, H0, f0), compose(ring, f1, J0))
    


if __name__ == "__main__":


    fn = argv.next() or "test"
    fn = eval(fn)
    fn()

    print("OK")


