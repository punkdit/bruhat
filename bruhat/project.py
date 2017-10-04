#!/usr/bin/env python

from __future__ import print_function

import sys, os
from fractions import Fraction
from random import shuffle

import numpy
from numpy.linalg import norm

from bruhat.argv import argv 
from bruhat.util import cross, factorial, write
from bruhat import sparse

from bruhat.sparse import Sparse as Op
from bruhat.sparse import Subspace


EPSILON = 1e-8

def mulclose(gen, maxsize=None, verbose=False):

    found = set(gen)
    els = list(gen)
    changed = True 
    while changed:
        if verbose:
            print("mulclose:", len(els))
        changed = False
        _els = list(els)
        for A in _els:
            for B in _els:
                C = A*B
                if C not in found:
                    els.append(C)
                    found.add(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
                    changed = True 
    return els



class Tensor(object):

    def __init__(self, ops):
        assert ops
        self.ops = list(ops)
        self.n = len(ops)
        for op in ops:
            assert op.shape[0] == op.shape[1]
        self.shape = tuple(len(op) for op in ops)

    def __mul__(self, other):
        assert self.n == other.n
        assert self.shape == other.shape
        ops = [A*B for A, B in zip(self.ops, other.ops)]
        return Tensor(ops)

    def __pow__(A, n):
        assert n>=1
        B = A
        while n>1:
            B = B*A
            n -= 1
        return B

    def tensor(self, other):
        ops = self.ops + other.ops
        return Tensor(ops)

    def __rmul__(self, phase):
        assert phase==1 or r==-1
        ops = [phase*self.ops[0]] + self.ops[1:]
        return Tensor(ops)

    def __neg__(self):
        return (-1)*self

    def __eq__(self, other):
        phase = +1
        for A, B in zip(self.ops, other.ops):
            r = A.eq_phase(B)
            if r is None:
                return False
            phase *= r
        return (phase==1)

    def __ne__(self, other):
        return not (self==other)

    def get_canonical(self):
        "put the phase in the first operator"
        ops = list(self.ops)
        phase = 1
        for i in range(1, self.n):
            op = self.ops[i]
            #r = op.signs[0]
            idxs = op.get_cols(0)
            r = op[idxs[0], 0]
            if r == -1:
                phase *= -1
                ops[i] = -op
        ops[0] = phase*self.ops[0]
        return Tensor(ops)

    _hash = None
    def __hash__(self):
        if self._hash is not None:
            return self._hash
        op = self.get_canonical()
        self._hash = hash(tuple(op.ops))
        return self._hash

    def todense(self):
        A = self.ops[0].tensor(*self.ops[1:])
        return A.todense()



def build_An(n):
    "Return list of generators for orthogonal reflection group A_{n-1}"
    assert n>=2

    gen = []

    # basis-swaps, "_controlled bitflips"
    for i in range(n-1):
        perm = list(range(n))
        perm[i], perm[i+1] = perm[i+1], perm[i]
        X = Op(n, n)
        for j in range(n):
            if j == i:
                X[j, j+1] = 1
            elif j == i+1:
                X[j, j-1] = 1
            else:
                X[j, j] = 1
        gen.append(X)

    return gen


def build_Bn(n):
    "Return list of generators for orthogonal reflection group B_n"
    assert n>=2

    gen = []

    # basis-swaps, "_controlled bitflips"
    for i in range(n-1):
        perm = list(range(n))
        perm[i], perm[i+1] = perm[i+1], perm[i]
        X = Op(n, n)
        for j in range(n):
            if j == i:
                X[j, j+1] = 1
            elif j == i+1:
                X[j, j-1] = 1
            else:
                X[j, j] = 1
        gen.append(X)

    # sign-swap, "_controlled phase-flip"
    Z = Op.identity(n)
    Z[n-1, n-1] = -1
    gen.append(Z)

    return gen


def build_Dn(n):
    "Return list of generators for orthogonal reflection group D_n"
    assert n>=2

    gen = []

    for i in range(n-1):
        perm = list(range(n))
        perm[i], perm[i+1] = perm[i+1], perm[i]
        X = Op(n, n)
        for j in range(n):
            if j == i:
                X[j, j+1] = 1
            elif j == i+1:
                X[j, j-1] = 1
            else:
                X[j, j] = 1
        gen.append(X)

    Z = Op.identity(n)
    Z[n-1, n-1] = -1
    Z[n-2, n-2] = -1
    gen.append(Z)

    return gen


def test():

    for n in [2, 3, 4]:

        ops = build_An(n)
        assert len(mulclose(ops))==factorial(n)
    
        ops = build_Bn(n)
        assert len(mulclose(ops))==2**n*factorial(n)
    
        ops = build_Dn(n)
        assert len(mulclose(ops))==2**(n-1)*factorial(n)
    

    from numpy import kron as tensor

    def allclose(A, B):
        return numpy.abs(A - B).sum()==0

    for A in ops:
      for B in ops:

        assert (A==B) == (hash(A)==hash(B))
        assert (A!=B) != (A==B)

        lhs = (A*B).todense()
        rhs = ((numpy.dot(A.todense(), B.todense())))
        assert allclose(lhs, rhs)

        lhs = (A.tensor(B)).todense()
        rhs = ((tensor(A.todense(), B.todense())))
        assert allclose(lhs, rhs)

    for A in ops:
        assert allclose(A.transpose().todense(), A.todense().transpose()), str(A)

    ops = mulclose(build_Bn(2))
    tops = []
    for A in ops:
      for B in ops:
        C = Tensor([A, B])
        tops.append(C)

    #print(len(tops))
    for A in tops:
        assert A.get_canonical() == A

    for A in tops:
      for B in tops:
        assert (A==B) == allclose(A.todense(), B.todense())
        assert (A==B) == (A.get_canonical()==B.get_canonical())
        assert (A==B) == (hash(A)==hash(B))

    print("OK")


def tensor(*ops):
    return ops[0].tensor(*ops[1:])


def order(g):
    n = 1
    a = g*g
    while a!=g:
        a = a*g
        n += 1
    #print(n, end=' ')
    return n


def inverse(g):
    n = order(g)
    if n==1:
        return g
    return g**(n-1)


def posproj(g, I):
    assert (g*g) == I
    p = Fraction(1, 2) * (I + g)
    assert (p*p) == p
    return p


def negproj(g, I):
    assert (g*g) == I
    p = Fraction(1, 2) * (I - g)
    assert (p*p) == p
    return p
    

def search(projs, n):

    i = len(projs)
    for p in projs:
      for q in projs:
        pq = p*q
        if pq == (q*p):
            yield [p, q]


def main():


    desc = argv.next() or "B2"
    n = int(desc[-1:])

    if desc[0] == "A":
        gen = build_An(n)
    elif desc[0] == "B":
        gen = build_Bn(n)
    elif desc[0] == "D":
        gen = build_Dn(n)
    else:
        print("no group desc found matching %s"%desc)
        return

    I = Op.identity(n)
    zero = Op(n)

    G = mulclose(gen)
    print(len(G))

    projs = []
    for g in G:
        if g==I:
            continue
        if g*g==I:
            projs.append(posproj(g, I))

    assert(len(set(projs))==len(projs))

    #for H in search(projs):
    #    print("found")

    pairs = []
    for p in projs:
      for q in projs:
        if p*q != q*p:
            pairs.append((p, q))
    print("pairs:", len(pairs))
    shuffle(pairs)

    for (p1, q1) in pairs:
      for (p2, q2) in pairs:
        A = tensor(p1, p2)
        B = tensor(q1, q2)
        if A*B == B*A:
            write(".")
    print()
      



if __name__ == "__main__":

    main()

