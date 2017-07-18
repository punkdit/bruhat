#!/usr/bin/env python

from __future__ import print_function

import sys, os

import numpy
from numpy import kron as tensor
from numpy.linalg import norm

from bruhat.argv import argv 
from bruhat.util import cross, factorial

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


class Op(object):
    """
        B_n operator: a signed perumutation matrix.
    """
    def __init__(self, perm, signs):
        self.perm = list(perm)
        self.signs = list(signs)
        self.n = len(self.perm)
        assert len(self.signs)==self.n
        self._hash = hash(str((self.perm, self.signs)))

    @classmethod
    def identity(cls, n):
        perm = range(n)
        signs = [1]*n
        return cls(perm, signs)

    def __str__(self):
        return str((self.todense()))

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

    def transpose(self):
        n = self.n
        perm = [None]*n
        signs = [None]*n
        for row, col in enumerate(self.perm):
            signs[col] = self.signs[row]
            perm[col] = row
        return Op(perm, signs)

    def __neg__(self):
        signs = [-sign for sign in self.signs]
        return Op(self.perm, signs)

    def __eq__(self, other):
        return self.perm==other.perm and self.signs==other.signs

    def __ne__(self, other):
        return self.perm!=other.perm or self.signs!=other.signs

    def __hash__(self):
        return self._hash

    def tensor(self, other, *rest):
        if rest:
            head, tail = rest[0], rest[1:]
            return self.tensor(other.tensor(head), *tail)

        m, n = self.n, other.n
        perm = []
        signs = []
        for i in range(m):
          for j in range(n):
            #k = i*n + j
            k = self.perm[i]*n + other.perm[j]
            perm.append(k)
            signs.append(self.signs[i] * other.signs[j])
        return Op(perm, signs)

    def todense(self):
        assert self.n < 2**10
        A = numpy.zeros((self.n, self.n))
        for i, idx in enumerate(self.perm):
            A[i, idx] = self.signs[i]
        return A

    def stab(self):
        n = self.n
        found = [0]*n
        orbits = []
        perm = self.perm
        signs = self.signs
        for i in range(n):
            if found[i]:
                continue
            parity = 1
            j = i
            orbit = [i]
            while 1:
                found[j] = 1
                parity *= signs[j]
                j = perm[j]
                if j==i:
                    break
                orbit.append(j)
            if parity == 1:
                orbit.sort()
                orbit = tuple(orbit)
                orbits.append(orbit)
        return orbits


def test():

    perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 1, 0], [2, 0, 1]]
    signss = [(a, b, c) for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
    ops = [Op(perm, signs) for perm in perms for signs in signss]

    for A in ops:
      for B in ops:

        assert (A==B) == (hash(A)==hash(B))
        assert (A!=B) != (A==B)

        lhs = (A*B).todense()
        rhs = ((numpy.dot(A.todense(), B.todense())))
        assert numpy.allclose(lhs, rhs)

        lhs = (A.tensor(B)).todense()
        rhs = ((tensor(A.todense(), B.todense())))
        assert numpy.allclose(lhs, rhs)

    for A in ops:
        assert numpy.allclose(A.transpose().todense(), A.todense().transpose()), str(A)


test()



def build_Bn(n):
    "Return list of generators for orthogonal reflection group B_n"
    assert n>=2

    gen = []

    # basis-swaps, "_controlled bitflips"
    for i in range(n-1):
        perm = list(range(n))
        perm[i], perm[i+1] = perm[i+1], perm[i]
        X = Op(perm, [1]*n)
        gen.append(X)

    # sign-swap, "_controlled phase-flip"
    signs = [1]*n
    signs[n-1] = -1
    Z = Op(range(n), signs)
    gen.append(Z)

    return gen


def toric(X, Z):

    n = 8
    N = X.n
    I = Op.identity(N)

    print("toric: dim=%d"%(N**n))

    plaqs = [(0, 1, 2, 5), (0, 2, 3, 7), (4, 5, 6, 1)] #, (6, 7, 3, 4)]
    stars = [(0, 1, 3, 4), (1, 2, 3, 6), (0, 4, 5, 7)] #, (2, 5, 6, 7)]

    for a in plaqs:
      for b in stars:
        assert len(set(a).intersection(set(b))) % 2 == 0

    ops = []
    for idxs in plaqs:
        op = [I] * n
        for idx in idxs:
            op[idx] = Z
        op = tensor(*op)
        ops.append(op)

    for idxs in stars:
        op = [I] * n
        for idx in idxs:
            op[idx] = X
        op = tensor(*op)
        ops.append(op)

    for A in ops:
      for B in ops:
        assert A*B == B*A

    v = numpy.zeros(N**n)
    v[0] = 1.
    for A in ops:
        v = v + A(v)
    assert(abs(norm(v)) > 0.1)
    v /= norm(v)
    #print(v)

    for A in ops:
        orbits = A.stab()
        print(A.n)
        for orbit in orbits:
            print("\t", orbit)
        u = A(v)
        r = numpy.dot(u, v)
        err = norm(u - r*v)
        print(r, err, end=' ')
    print()
    


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


def main():

    X, Z = build_Bn(2)
    I = Op.identity(2)

    #toric(X, Z)
    #return

    #print(I)
    #print(X)
    #print(Z)

    II = Op.identity(4)

    P2 = [X.tensor(I), Z.tensor(I), I.tensor(X), I.tensor(Z)]
    P2 = mulclose(P2)
    
    assert len(P2)==32
    assert II in P2

    errors = []
    for A in P2:
        if A == II or A == -II or -A in errors:
            continue
        errors.append(A)
    assert len(errors)==15

    if 0:
        G = mulclose(P2)
        for g in G:
            A = g.todense()
            vals = numpy.linalg.eigvals(A)
            print(vals) # all +1, -1

    d = argv.get("d", 4)
    gen = build_Bn(d)
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
        assert len([g for g in G if g in P2]) == 32

    pairs = []

    c_comm = 0
    c_anti = 0
    total = 0
    for i in range(len(G)):
      for j in range(len(G)):
        g = G[i]
        h = G[j]

        gh = g*h
        hg = h*g

        total += 1
        if (gh == -hg):
            c_anti += 1
        elif (gh == hg):
            c_comm += 1

        if (g) in P2 or (h) in P2:
            continue

        if (gh == -hg) and i<j:

            #assert order(g)==2
            print(order(g), end=' ')
            #print(g)
            #print(h)
            #print((g == g.transpose()), end=' ')
            #print((h == h.transpose()))
            if (g == g.transpose()) and (h == h.transpose()) and i < j:
                pairs.append((g, h))

        #print(".", end='')
    #print()
    #print(c_comm, c_anti, total)

    print("pairs:", len(pairs))
    for (g, h) in pairs:
        ok = False
        for E in errors:
            if g*E != E*g or h*E != E*h:
                break
        else:
            print("fail")
            continue
        toric(g, h)
        break

    return

    ops = set()
    for (g, h) in pairs:
        ops.add(g)
        ops.add(h)
    print(len(ops))

    for g in ops:
        for h in ops:
            a = (g*h == h*g)
            b = (g*h == -h*g)
            if a==b==False:
                s = '  '
            elif a:
                s = '+ '
            elif b:
                s = '- '
            print(s, end=' ')
        print()

    
    h = tensor(*([g]*8))
    print(h.n)


if __name__ == "__main__":

    main()


