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
from bruhat.elim import eq
#from bruhat.elim import (
#    zeros, rand, dot, identity, coequalizer, compose,
#    rank, pseudo_inverse)
#from bruhat.solve import parse
from bruhat import solve


def shortstr(A):
    s = elim.shortstr(A)
    s = str(s)
    s = s.replace(" 0 ", " . ")
    return s


class Space(object):
    def __init__(self, ring):
        self.ring = ring

    def parse(self, decl):
        ring = self.ring
        A = solve.parse(decl)
        B = self.zeros(*A.shape)
        B = B+A
        return B

    def zeros(self, *args):
        return elim.zeros(self.ring, *args)
    
    def rand(self, *args):
        return elim.rand(self.ring, *args)
    
    def dot(self, *args):
        return elim.dot(self.ring, *args)
    
    def identity(self, *args):
        return elim.identity(self.ring, *args)
    
    def coequalizer(self, *args):
        return elim.coequalizer(self.ring, *args)
    
    def compose(self, *args):
        return elim.compose(self.ring, *args)
    
    def rank(self, *args):
        return elim.rank(self.ring, *args)
    
    def pseudo_inverse(self, *args):
        return elim.pseudo_inverse(self.ring, *args)
    
    def kron(self, A, B):
        ring = self.ring
        if 0 in A.shape or 0 in B.shape:
            C = self.zeros(A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])
        else:
            #print("kron", A.shape, B.shape)
            C = numpy.kron(A, B)
            #print("\t", C.shape)
        return C
    
    def sum_swap(self, a, b):
        "swap isomorphism a+b -> b+a"
        ring = self.ring
        one = ring.one
        g = self.zeros(a+b, a+b)
        for i in range(a):
            g[i, i+a] = one
        for i in range(b):
            g[i+b, i] = one
        return g
    
    def direct_sum(self, f, g):
        ring = self.ring
        mf, nf = f.shape
        mg, ng = g.shape
        h = self.zeros(mf+mg, nf+ng)
        h[:mf, :nf] = f
        h[mf:, nf:] = g
        return h
    
    def tensor_swap(self, m1, m2):
        ring = self.ring
        A = self.zeros(m2*m1, m1*m2)
        one = ring.one
        for i in range(m1):
          for j in range(m2):
            A[j*m1 + i , i*m2 + j] = one
        return A

    def schur(self, m):
        I = self.identity(m*m)
        s = self.tensor_swap(m, m)
        f = self.coequalizer(I, s)
        return f
    
    def is_zero(self, A):
        return eq(A, self.zeros(*A.shape))
    
    def is_identity(self, A):
        return eq(A, self.identity(A.shape[0]))


class Chain(object):
    def __init__(self, space, diffs={}):
        # map domain grade to array
        assert type(diffs) is dict
        if not diffs:
            diffs = { 0 : space.zeros(0, 0) } # ?
        grades = list(diffs.keys())
        grades.sort()
        dims = {}
        for grade in grades:
            A = diffs[grade]
            m, n = A.shape
            dim = dims.setdefault(grade, n)
            assert dim==n, "mismatch at grade %d"%(grade,)
            dim = dims.setdefault(grade-1, m)
            assert dim==m, "mismatch at grade %d"%(grade-1,)
    
        self.diffs = dict(diffs)
        self.grades = grades
        self.mingrade = min(grades)
        self.maxgrade = max(grades) + 1
        self.dims = dims

    def __getitem__(self, grade):
        dims = self.dims
        space = self.space
        A = self.diffs.get(grade)
        if A is None:
            n = dims.get(grade, 0)
            m = dims.get(grade-1, 0)
            A = space.zeros(m, n)
            self.diffs[grade] = A
        return A

    #def __setitem__(self, grade, A):
    #    self.diffs[grade] = A

    def dual(self):
        diffs = {}
        for grade in self.grades:
            A = self.diffs[grade]
            A = A.transpose()
            diffs[grade - 1] = A
        return Chain(diffs)

    def direct_sum(self, other):
        space = self.space
        diffs = {}
        grades = list(set(self.grades + other.grades))
        grades.sort()
        for grade in grades:
            A, B = self[grade], other[grade]
            diffs[grade] = space.direct_sum(A, B)
        return Chain(diffs)

    def tensor(self, other):
        space = self.space
        tensor = space.tensor
        diffs = {}
        g0 = min(self.mingrade, other.mingrade)
        g1 = max(self.maxgrade, other.maxgrade)
        for grade in range(g0, g1+1):
            # map grade -> grade-1
            rows, cols = [], []
            for a in range(g0, g1+1):
                b = grade - a
                A = self[a]
                B = other[b]
                 # ARGGHH


def test():

    p = argv.get("p", 2)
    ring = element.FiniteField(p)
    space = Space(ring)
    zeros = space.zeros
    rand = space.rand
    dot = space.dot
    kron = space.kron
    direct_sum = space.direct_sum
    identity = space.identity
    coequalizer = space.coequalizer
    compose = space.compose
    rank = space.rank
    pseudo_inverse = space.pseudo_inverse
    tensor_swap = space.tensor_swap
    sum_swap = space.sum_swap
    schur = space.schur
    is_zero = space.is_zero
    is_identity = space.is_identity

    s = tensor_swap(3, 4)
    si = tensor_swap(4, 3)
    #print(shortstr(s))
    assert eq(dot(si, s), identity(3*4))
    assert eq(dot(s, si), identity(4*3))

    m, n = 2, 3
    A1 = rand(m, n, 1, 1)
    A2 = rand(m, n, 1, 1)

    B = kron(A1, A2)

    for m in range(1, 5):
        I = identity(m*m)
        s = tensor_swap(m, m)
        f = coequalizer(I, s)
    
        assert eq(compose(s, f), f)
        assert rank(f) == [1, 3, 6, 10][m-1]

    # ---------------------------------

    m = argv.get("m", 3)
    n = argv.get("n", 4)

    if argv.toric:
        A = zeros(m, m)
        for i in range(m):
            A[i, i] = ring.one
            A[i, (i+1)%m] = -ring.one
    elif argv.surface:
        A = zeros(m-1, m)
        for i in range(m-1):
            A[i, i] = ring.one
            A[i, (i+1)%m] = -ring.one
    else:
        A = rand(m, n, p-1, p-1)
    if argv.transpose:
        A = A.transpose()

    print("A:")
    print(shortstr(A))


    n, m = A.shape

    In = identity(n)
    Im = identity(m)

    H1s = kron(Im, A), -kron(A, Im)
    H1 = numpy.concatenate(H1s, axis=0) # horizontal concatenate

    H0s = kron(A, In), kron(In, A)
    H0 = numpy.concatenate(H0s, axis=1) # horizontal concatenate

    assert is_zero(dot(H0, H1))

    assert H1.shape == (n*m+m*n, m*m)
    assert H0.shape == (n*n, n*m+m*n)

    f0 = -tensor_swap(n, n)
    a = direct_sum( -tensor_swap(m, n), -tensor_swap(n, m))
    b = sum_swap(n*m, m*n)
    assert is_identity(compose(b, b))
    f1 = compose(a, b)
    assert is_identity(compose(f1, f1))
    f2 = tensor_swap(m, m)

    assert eq(compose(f2, H1), compose(H1, f1))
    lhs, rhs = ( (compose(f1, H0), compose(H0, f0)) )
    #print("lhs:")
    #print(shortstr(lhs))
    #print("rhs:")
    #print(shortstr(rhs))
    assert eq(compose(f1, H0), compose(H0, f0))

    g0 = coequalizer(
        f0, identity(f0.shape[0]))

    assert eq(compose(H0, g0), compose(f1, H0, g0))

    e = compose(H0, g0)
    g1, J0 = coequalizer(f1, identity(f1.shape[0]), e)

    assert eq(compose(H0, g0), compose(g1, J0))

    e = compose(H1, g1)
    g2, J1 = coequalizer(f2, identity(f2.shape[0]), e)

    assert eq(compose(H1, g1), compose(g2, J1))

    assert is_zero(compose(J1, J0))

    n = J1.shape[0]
    J1t = J1.transpose()
    mz = rank(J1t)
    mx = rank(J0)
    print("J1t:", J1t.shape, rank(J1t))
    print(shortstr(J1t))
    print("J0:", J0.shape)
    print(shortstr(J0))

    print("n:", n)
    print("mz:", mz)
    print("mx:", mx)
    print("k =", n-mx-mz)



if __name__ == "__main__":


    fn = argv.next() or "test"
    fn = eval(fn)
    fn()

    print("OK")


