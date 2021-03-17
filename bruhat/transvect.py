#!/usr/bin/env python3


from random import shuffle
from functools import reduce
from operator import matmul, add, mul

import numpy

from bruhat import element
from bruhat import elim

#from bruhat.vec import Space, Hom, Map
#from bruhat.action import mulclose, mulclose_names, Perm, Group
from bruhat.argv import argv


p = argv.get("p", 2)
ring = element.FiniteField(p)
one = ring.one
zero = ring.zero

def shortstr(A):
    s = elim.shortstr(A)
    s = str(s)
    s = s.replace(" 0 ", " . ")
    return s

def array(A):
    A = numpy.array(A, dtype=object)
    for idx in numpy.ndindex(A.shape):
        A[idx] = ring.promote(A[idx])
    return A

def zeros(*shape):
    A = numpy.empty(shape, dtype=object)
    A[:] = zero
    return A

def identity(m):
    A = numpy.empty((2*m, 2*m), dtype=object)
    A[:] = zero
    for i in range(2*m):
        A[i,i] = one
    return A

_cache = {}
def omega(n):
    if n in _cache:
        return _cache[n]
    A = zeros(2*n, 2*n)
    for i in range(n):
        A[i, i+n] = one
        A[i+n, i] = -one
    _cache[n] = A # warning: not immutable!
    return A

def rowvec(A):
    if len(A.shape)==1:
        A = A.view()
        A.shape = (1,)+A.shape
    assert len(A.shape)==2
    return A

def sy(A, B):
    "inner product of row-vectors"
    A = rowvec(A)
    B = rowvec(B)
    n = A.shape[1]
    assert n == B.shape[1] 
    assert n%2==0
    O = omega(n//2)
    Bt = B.transpose()
    r = numpy.dot(A, numpy.dot(O, Bt))
    return r

def transvect(v):
    v = rowvec(v)
    m = v.shape[1]//2
    O = omega(m)
    I = identity(m)
    vtv = numpy.dot(v.transpose(), v)
    Tv = I + numpy.dot(O, vtv)
    return Tv
    

def toint(A):
    B = numpy.zeros(A.shape, dtype=int)
    for idx in numpy.ndindex(A.shape):
        B[idx] = A[idx].value # unwrap
    return B


def upper(A):
    A = A.copy()
    for (i, j) in numpy.ndindex(A.shape):
        if i>=j:
            A[i, j] = 0
    return A



def main():

    #F = array([[1,0], [0,1]])
    #print(shortstr(F))

    m = argv.get("m", 5)

    O = omega(m)

    V = array([
        [1,0,0,1,1, 1,0,0,1,0],
        [1,1,0,0,0, 1,1,1,1,0],
        [0,1,1,0,1, 1,0,0,1,0],
        [0,0,0,1,0, 1,1,1,1,0],
        [0,1,0,0,0, 1,0,0,0,1],
    ])

    print(shortstr(V))
    m = len(V)

    F0 = reduce(numpy.dot, [transvect(v) for v in V])
    print("F0=")
    print(shortstr(F0))

    A = zeros(m, m)
    A = sy(V, V)
    print(shortstr(A))

    A = toint(A)
    Au = upper(A)
    print(Au)

    B = numpy.zeros((m, m), dtype=int)
    P = numpy.identity(m, dtype=int)
    while numpy.any(P):
        B += P
        P = numpy.dot(Au, P)
    print(B)

    B = array(B)
    print(shortstr(B))

    Vt = V.transpose()
    F = identity(m) + numpy.dot(O, numpy.dot(Vt, numpy.dot(B, V)))

    print("F=")
    print(shortstr(F))

    #print(F==F0)

    assert numpy.all(F0==F)


if __name__ == "__main__":

    main()





