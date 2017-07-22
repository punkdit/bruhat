#!/usr/bin/env python

from __future__ import print_function

import sys, os

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
    "Return list of generators for orthogonal reflection group A_n"
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



class Toric(object):
    def __init__(self, l, X, Z):
        assert l >= 2

        N = X.n
        I = Op.identity(N)
        assert X*X == I
        assert Z*Z == I
        assert Z*X == -X*Z

        keys = []
        lookup = {}
        for i in range(l):
          for j in range(l):
            for k in range(2):
                key = (i, j, k)
                idx = len(keys)
                keys.append(key)
                for di in (-l, 0, l):
                  for dj in (-l, 0, l):
                    lookup[i+di, j+dj, k] = idx

        n = 2*(l**2) # qubits

        II = Tensor([I]*n)

        ops = []
        for i in range(l):
          for j in range(l):

            # plaqs
            op = [I]*n
            for idx in [lookup[i, j, 0], lookup[i, j+1, 0], lookup[i, j, 1], lookup[i+1, j, 1]]:
                op[idx] = Z
            op = Tensor(op)
            assert op*op==II
            ops.append(op)
        
            # stars
            op = [I]*n
            for idx in [lookup[i, j, 0], lookup[i-1, j, 0], lookup[i, j, 1], lookup[i, j-1, 1]]:
                op[idx] = X
            op = Tensor(op)
            assert op*op==II
            ops.append(op)

        for A in ops:
          for B in ops:
            assert A*B==B*A

        if l==2:
            assert len(mulclose(ops))==64
        

def toric(Xs, Zs):

    n = 8
    N = Xs[0].n
    I = Op.identity(N)

#    assert X*X == I
#    assert Z*Z == I
#    assert Z*X == -X*Z
#    ZZ = tensor(Z, Z)
#    XX = tensor(X, X)
#    assert ZZ*XX == XX*ZZ

    X1, X2 = Xs
    Z1, Z2 = Zs
    assert X1*X2 == X2*X1
    assert Z1*Z2 == Z2*Z1
    assert X1*Z2 == Z2*X1
    assert Z1*X2 == X2*Z1
    assert X1*Z1 == -Z1*X1
    assert X2*Z2 == -Z2*X2

    print("toric: dim=%d"%(N**n))
    II = Op.identity(N**n)

    plaqs = [(0, 1, 2, 5), (0, 2, 3, 7), (4, 5, 6, 1), (6, 7, 3, 4)]
    stars = [(0, 1, 3, 4), (1, 2, 3, 6), (0, 4, 5, 7), (2, 5, 6, 7)]

    for a in plaqs:
      for b in stars:
        assert len(set(a).intersection(set(b))) % 2 == 0

    ops = []
    for idxs in plaqs:
        for Z in Zs:
            op = [I] * n
            for idx in idxs:
                op[idx] = Z
            op = tensor(*op)
            #assert op*op == II
            ops.append(op)

    for idxs in stars:
        for X in Xs:
            op = [I] * n
            for idx in idxs:
                op[idx] = X
            op = tensor(*op)
            #assert op*op == II
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

    spaces = []
    for A in ops:
        S = A.get_stab()
        print("get_stab:", len(S))
        spaces.append(S)

        u = A(v)
        r = numpy.dot(u, v)
        err = norm(u - r*v)
        assert abs(r-1.) < EPSILON
        assert err < EPSILON
        #print(r, err, end=' ')
    #print()
    
    S = spaces[0]
    for S1 in spaces[1:]:
        print("codespace:", len(S))
        S = S.intersect(S1)
    print("codespace:", len(S))


def multi_toric(Xs, Zs):

    n = 8
    N = Xs[0].n
    I = Op.identity(N)
    for X in Xs:
      for Z in Zs:
        assert X*X == I
        assert Z*Z == I
        assert Z*X == -X*Z

    print("toric: dim=%d"%(N**n))
    II = Op.identity(N**n)

    plaqs = [(0, 1, 2, 5), (0, 2, 3, 7), (4, 5, 6, 1), (6, 7, 3, 4)]
    stars = [(0, 1, 3, 4), (1, 2, 3, 6), (0, 4, 5, 7), (2, 5, 6, 7)]

    for a in plaqs:
      for b in stars:
        assert len(set(a).intersection(set(b))) % 2 == 0

    ops = []
    for idxs in stars:
        for X in Xs:
            op = [I] * n
            for idx in idxs:
                op[idx] = X
            op = tensor(*op)
            assert op*op == II
            ops.append(op)

    for idxs in plaqs:
        for Z in Zs:
            op = [I] * n
            for idx in idxs:
                op[idx] = Z
            op = tensor(*op)
            assert op*op == II
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

    spaces = []
    for A in ops:
        S = A.get_stab()
        #print("get_stab:", len(S))
        spaces.append(S)

        u = A(v)
        r = numpy.dot(u, v)
        err = norm(u - r*v)
        assert abs(r-1.) < EPSILON
        assert err < EPSILON
        #print(r, err, end=' ')
    #print()
    
    S = spaces[0]
    for S1 in spaces[1:]:
        print("codespace:", len(S))
        S = S.intersect(S1)
    print("codespace:", len(S))

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


def main():

    X, Z = build_Bn(2)
    I = Op.identity(2)

    l = argv.get("l", 2)
    if argv.B2:
        #toric(X, Z)
        code = Toric(l, X, Z)
        print("OK")
        return

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


    H = [g for g in G if g not in P2]

    for ops in all_quads(H):
        X1, X2, Z1, Z2 = ops
        ss = set()
        print("quad") # no quads found ...
        for E in errors:
            #s = tuple(order(E*A) for A in ops)
            #ss.append(s)
            ss.add(tuple(E*A*inverse(E)*inverse(A) for A in ops))
            #ss.add(tuple(order(E*A) for A in ops))
        #print(len(set(ss)))
        if len(set(ss))==len(errors):
            write("/")
            toric([X1, X2], [Z1, Z2])
        else:
            write("\\")


def all_quads(H):
    n = len(H)
    quads = []
    for i1 in range(n):
      for j1 in range(i1+1, n):
        X1 = H[i1]
        Z1 = H[j1]
#        if order(X1*Z1)!=4:
#            continue
        if X1*Z1 != -Z1*X1:
            continue

        for i2 in range(j1+1, n):
            X2 = H[i2]
#            if order(X1*X2)!=2:
#                continue
#            if order(Z1*X2)!=2:
#                continue
            if X1*X2 != X2*X1:
                continue
            if Z1*X2 != X2*Z1:
                continue

            for j2 in range(i2+1, n):
                Z2 = H[j2]
#                if order(X1*Z2)!=2:
#                    continue
#                if order(Z1*Z2)!=2:
#                    continue
#                if order(X2*Z2)!=4:
#                    continue
                if Z1*Z2 != Z2*Z1:
                    continue
                if X1*Z2 != Z2*X1:
                    continue
                if X2*Z2 != -Z2*X2:
                    continue
                yield (X1, X2, Z1, Z2)


def show_table(H):
    for g in H:
        assert order(g)==2
        for h in H:
            i = order(g*h)
            if i==2:
                s = '+ '
            elif i==4:
                s = '- '
            elif i==1:
                s = '  '
            else:
                s = '%d '%i
            print(s, end=' ')
        print()

    return

def scrap():

    for (X, Z) in pairs:
        ok = False
        for E in errors:
#            if X*E != E*X or Z*E != E*Z:
            if X*E == -E*X or Z*E == -E*Z:
                break
        else:
            print("fail")
            continue
        #toric(X, Z)

        code = Toric(l, X, Z)

        print("syndromes:")
        for E in errors:
            #print(code.get_syndrome(E))
            print("  ", end='')
            for op in [X, Z]:
                if op*E==E*op:
                    s = '.'
                if op*E==-E*op:
                    s = '-'
                else:
                    s = '?'
                print(s, end='')
        print()

        #return

    

if __name__ == "__main__":

    if argv.test:
        test()

    elif argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:
        main()


