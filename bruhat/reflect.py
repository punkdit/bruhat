#!/usr/bin/env python

from __future__ import print_function

import sys, os

import numpy
from numpy.linalg import norm

from bruhat.argv import argv 
from bruhat.util import cross, factorial
from bruhat import sparse
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

    def __rmul__(self, r):
        assert r==1 or r==-1
        if r==1:
            return self
        else:
            return -self

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

    def eq_phase(self, other):
        assert self.n == other.n
        if self.perm != other.perm:
            return None
        phase = self.signs[0] * other.signs[0]
        for a, b in zip(self.signs, other.signs):
            if a*b != phase:
                return None
        return phase

    def get_stab(self):
        "Get the stabilized subspace (+1 eigenspace)"
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
            if parity != 1:
                continue
            orbit.sort()
            orbit = tuple(orbit)
            orbits.append(orbit)
        A = sparse.zeros(len(orbits), n)
        for i, orbit in enumerate(orbits):
            for j in orbit:
                A[i, j] = 1
        S = Subspace(A)
        return S


class Tensor(object):

    def __init__(self, ops):
        assert ops
        self.ops = list(ops)
        self.n = len(ops)
        self.shape = tuple(op.n for op in ops)

    def __mul__(self, other):
        assert self.n == other.n
        assert self.shape == other.shape
        ops = [A*B for A, B in zip(self.ops, other.ops)]
        return Tensor(ops)

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
            r = op.signs[0]
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


def test():

    from numpy import kron as tensor
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
        assert (A==B) == numpy.allclose(A.todense(), B.todense())
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
        

def toric(X, Z):

    n = 8
    N = X.n
    I = Op.identity(N)
    assert X*X == I
    assert Z*Z == I
    assert Z*X == -X*Z
    ZZ = tensor(Z, Z)
    XX = tensor(X, X)
    assert ZZ*XX == XX*ZZ

    print("toric: dim=%d"%(N**n))
    II = Op.identity(N**n)

    plaqs = [(0, 1, 2, 5), (0, 2, 3, 7), (4, 5, 6, 1), (6, 7, 3, 4)]
    stars = [(0, 1, 3, 4), (1, 2, 3, 6), (0, 4, 5, 7), (2, 5, 6, 7)]

    for a in plaqs:
      for b in stars:
        assert len(set(a).intersection(set(b))) % 2 == 0

    ops = []
    for idxs in plaqs:
        op = [I] * n
        for idx in idxs:
            op[idx] = Z
        op = tensor(*op)
        assert op*op == II
        ops.append(op)

    for idxs in stars:
        op = [I] * n
        for idx in idxs:
            op[idx] = X
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

#    print("orders:")
#    for g in G:
#        print(order(g), end=' ')
#    print()

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
            if (g == g.transpose()) and (h == h.transpose()) and i < j:
                assert order(g)==2
                pairs.append((g, h))

    #print(c_comm, c_anti, total)
    print("pairs:", len(pairs))
    
    ops = set()
    for (g, h) in pairs:
        ops.add(g)
        ops.add(h)
    print("ops:", len(ops))

    assert len([g for g in G if g==g.transpose() and order(g)==2])==75

    ops = list(ops)
    X = ops[0]
    H = [g for g in ops if order(g*X) in [1, 2, 4]]
    assert len(H)==8

    Xs = [g for g in H if order(g*X) in [1, 2]]
    Zs = [g for g in H if order(g*X)==4]
    assert len(Xs)==4
    assert len(Zs)==4

    #multi_toric(Xs, Zs)

    print("syndromes:")
    bag = set()
    for E in errors:
    #for E in G:
        #print(code.get_syndrome(E))
        syndrome = ''
        for op in Xs+Zs:
            #s = str(order(op*E))
            bag.add(op*E*op*E)
            #syndrome += s
        #print(syndrome)
    print(len(bag))


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


