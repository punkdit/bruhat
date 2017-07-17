#!/usr/bin/env python

from __future__ import print_function

import sys, os

import numpy
from numpy import kron as tensor

from bruhat.argv import argv 
from bruhat.util import cross, factorial
#from bruhat import action
#from bruhat.action import Perm, Group, mulclose
from bruhat.gelim import solve, array, identity, zeros, dot, shortstr, eq, dotx, kernel
from bruhat.gelim import Subspace


def mulclose(gen, maxsize=None, verbose=False):

    found = set(str(g) for g in gen)
    els = list(gen)
    changed = True 
    while changed:
        if verbose:
            print("mulclose:", len(els))
        changed = False
        _els = list(els)
        for A in _els:
            for B in _els:
                C = dot(A, B)
                if str(C) not in found:
                    els.append(C)
                    found.add(str(C))
                    if maxsize and len(els)>=maxsize:
                        return list(els)
                    changed = True 
    return els


class Op(object):
    """
        B_n operator: a signed perumutation matrix.
    """
    def __init__(self, perm, signs):
        self.perm = list(perm)
        self.signs = list(signs)
        self.n = len(self.perm)
        assert len(self.signs)==self.n

    def __str__(self):
        return str(shortstr(self.todense()))

    def __call__(self, v):
        assert len(v) == self.n
        u = numpy.zeros(self.n, dtype=v.dtype)
        for i, idx in enumerate(self.perm):
            sign = self.signs[i]
            u[i] = v[idx]
        return u

    def __mul__(self, other):
        assert isinstance(other, Op)
        perm = []
        signs = []
        for i, idx in enumerate(self.perm):
            perm.append(other.perm[idx])
            signs.append(self.signs[i] * other.signs[idx])
        return Op(perm, signs)

    def todense(self):
        assert self.n < 2**10
        A = numpy.zeros((self.n, self.n))
        for i, idx in enumerate(self.perm):
            A[i, idx] = self.signs[i]
        return A


def test():

    perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 1, 0], [2, 0, 1]]
    signss = [(a, b, c) for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
    ops = [Op(perm, signs) for perm in perms for signs in signss]

    for A in ops:
      for B in ops:
        lhs = (A*B).todense()
        rhs = ((dot(A.todense(), B.todense())))
        assert eq(lhs, rhs)
    print(len(ops))

#test()
#sys.exit(0)



def build(n):
    "Return list of generators for orthogonal reflection group B_n"
    assert n>=2

    gen = []

    # basis-swaps, "_controlled bitflips"
    for i in range(n-1):
        X = zeros(n, n)
        for j in range(n):
            if j==i:
                X[j, j+1] = 1
            elif j==i+1:
                X[j, j-1] = 1
            else:
                X[j, j] = 1
        gen.append(X)

    # sign-swap, "_controlled phase-flip"
    Z = identity(n)
    Z[n-1, n-1] = -1
    gen.append(Z)

    return gen


def order(g):
    n = 1
    a = dot(g, g)
    while not eq(a, g):
        a = dot(a, g)
        n += 1
    #print(n, end=' ')
    return n


def main():

    X, Z = build(2)
    I = dot(X, X)

    #print(I)
    #print(X)
    #print(Z)

    II = tensor(I, I)

    P2 = [tensor(X, I), tensor(Z, I), tensor(I, X), tensor(I, Z)]
    found = set()
    for g in mulclose(P2):
        s = shortstr(g)
        found.add(s)
    assert len(found)==32
    assert shortstr(II) in found

    if 0:
        G = mulclose(P2)
        for g in G:
            A = numpy.array(g, dtype=numpy.float)
            vals = numpy.linalg.eigvals(A)
            print(vals) # all +1, -1

    d = argv.get("d", 4)
    gen = build(d)
    n = (2**d) * factorial(d)
    assert len(mulclose(gen)) == n # slow...
    G = mulclose(gen, n)

    print("orders:")
    for g in G:
        print(order(g), end=' ')
    print()

    if 0:
        for g in G:
            A = numpy.array(g, dtype=numpy.float)
            vals = numpy.linalg.eigvals(A)
            print(vals) # all kinds of stuff..

    if d==4:
        assert len([g for g in G if shortstr(g) in found]) == 32

    pairs = []

    c_comm = 0
    c_anti = 0
    total = 0
    for i in range(len(G)):
      for j in range(len(G)):
        g = G[i]
        h = G[j]

        gh = dot(g, h)
        hg = dot(h, g)

        total += 1
        if eq(gh, -hg):
            c_anti += 1
        elif eq(gh, hg):
            c_comm += 1

        if shortstr(g) in found or shortstr(h) in found:
            continue

        if eq(gh, -hg) and i<j:

            #assert order(g)==2
            print(order(g), end=' ')
            #print(g)
            #print(h)
            #print(eq(g, g.transpose()), end=' ')
            #print(eq(h, h.transpose()))
            if eq(g, g.transpose()) and eq(h, h.transpose()) and i < j:
                pairs.append((g, h))

        #print(".", end='')
    #print()
    #print(c_comm, c_anti, total)

    print("pairs:", len(pairs))

    ops = {}
    for (g, h) in pairs:
        ops[shortstr(g)] = g
        ops[shortstr(h)] = h
    ops = list(ops.values())
    print(len(ops))

    for g in ops:
        for h in ops:
            a = eq(dot(g, h), dot(h, g))
            b = eq(dot(g, h), -dot(h, g))
            if a==b==False:
                s = '  '
            elif a:
                s = '+ '
            elif b:
                s = '- '
            print(s, end=' ')
        print()


    return


    A, B, C, D = gen
    pairs = [s.split() for s in ("ABDA CBDC", "BACB DCDA", "BACB DCDC",
        "BACDBCBADCBC DCABDCDABACD" )]
    for l, r in pairs:
        left = [eval(s) for s in l]
        right = [eval(s) for s in r]
        left = dotx(*left)
        right = dotx(*right)

        if 1:
            print(shortstr(left))
            print(shortstr(left) in found)
            print()
            print(shortstr(right))
            print(shortstr(right) in found)
            print()
            print(shortstr(dotx(left, right)))
            print()
            print()
        assert eq((dotx(right, left)), -(dotx(left, right)))




if __name__ == "__main__":

    main()


