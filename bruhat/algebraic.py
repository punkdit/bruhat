#!/usr/bin/env python3

"""
Algebraic groups: matrix groups over Z/pZ.

"""


import sys, os
import random

import numpy

scalar = numpy.int64

from bruhat.gset import Group, Perm, GSet
from bruhat.action import mulclose
from bruhat.spec import isprime
from bruhat.argv import argv
from bruhat.solve import parse, enum2, row_reduce, span
from bruhat.dev import geometry
from bruhat.util import cross


DEFAULT_P = 2

class Matrix(object):
    def __init__(self, A, p=DEFAULT_P, shape=None):
        if type(A) == list or type(A) == tuple:
            A = numpy.array(A, dtype=scalar)
        else:
            A = A.astype(scalar) # makes a copy
        if shape is not None:
            A.shape = shape
        self.A = A
        #n = A.shape[0]
        #assert A.shape == (n, n)
        assert int(p) == p
        assert p>=0
        self.p = p
        #self.n = n
        if p>0:
            self.A %= p
        self.key = (self.p, self.A.tostring())
        self._hash = hash(self.key)
        self.shape = A.shape

    def __str__(self):
        return str(self.A)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        assert self.p == other.p
        return self.key == other.key

    def __ne__(self, other):
        assert self.p == other.p
        return self.key != other.key

    def __add__(self, other):
        A = self.A + other.A
        return Matrix(A, self.p)

    def __sub__(self, other):
        A = self.A - other.A
        return Matrix(A, self.p)

    def __neg__(self):
        A = -self.A
        return Matrix(A, self.p)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            A = numpy.dot(self.A, other.A)
            return Matrix(A, self.p)
        else:
            return NotImplemented

    def mask(self, A):
        return Matrix(self.A * A, self.p)


# https://math.stackexchange.com/questions/34271/
# order-of-general-and-special-linear-groups-over-finite-fields

def order_gl(n, q):
    order = 1
    for i in range(n):
        order *= (q**n - q**i)
    return order

def order_sl(n, q):
    order = order_gl(n, q)
    assert order % (q-1) == 0
    return order//(q-1)


def order_sp(n, q):
    # n = 2*m
    assert n%2==0
    m = n//2
    N = q**(m**2)
    for i in range(m):
        N *= (q**(2*(i+1)) - 1)
    return N

assert order_sp(2, 2)==6     # 3!
assert order_sp(4, 2)==720   # 6!


def SL(n, p=DEFAULT_P):
    "special linear group"
    assert int(n)==n
    assert int(p)==p
    assert n>0
    assert isprime(p)

    I = numpy.identity(n, scalar)
    gen = []
    for i in range(n):
        for j in range(n):
            if i==j:
                continue
            A = I.copy()
            A[i, j] = 1
            gen.append(Matrix(A, p))
    order = order_sl(n, p)
    G = mulclose(gen, maxsize=order)
    #G = mulclose(gen)
    #assert len(G)==order
    return list(G)



def GL(n, p=DEFAULT_P):
    "general linear group"
    assert int(n)==n
    assert int(p)==p
    assert n>0
    assert isprime(p)

    H = SL(n, p)
    nI = Matrix(-numpy.identity(n, scalar), p)
    if p>2:
        G = H + [nI*g for g in H]
    else:
        G = H
    order = order_gl(n, p)
    assert len(G)==order
    return list(G)


        
def get_subgroup(G, desc, check=False):
    A = parse(desc)
    H = []
    for op in G:
        if op.mask(A) == op:
            H.append(op)

    if check:
        for a in H:
          for b in H:
            assert (a*b).mask(A) == (a*b)
    return H


def get_permrep(G):
    assert len(G)
    op = G[0]
    n, _ = op.shape
    p = op.p
    space = list(numpy.array(v) for v in enum2(n))
    lookup = dict((v.tostring(), idx) for (idx, v) in enumerate(space))
    perms = []
    for op in G:
        idxs = [lookup[(numpy.dot(op.A, v)%p).tostring()] for v in space]
        perms.append(Perm(idxs))
    G = Group(perms)
    return G


def test_dynkin():

    G = GL(4, 2)
    subgroups = []

    n = len(G)
    print(".--.--. =", len(G), n//len(G))

    POINT = "1111 .111 .111 .111"
    H = get_subgroup(G, POINT)
    print("*--.--. =", len(H), n//len(H))
    subgroups.append(H)

    LINE = "1111 1111 ..11 ..11"
    H = get_subgroup(G, LINE)
    print(".--*--. =", len(H), n//len(H))
    subgroups.append(H)

    PLANE = "1111 1111 1111 ...1"
    H = get_subgroup(G, PLANE)
    print(".--.--* =", len(H), n//len(H))
    subgroups.append(H)

    PonL = "1111 .111 ..11 ..11"
    H = get_subgroup(G, PonL)
    print("*--*--. =", len(H), n//len(H))
    subgroups.append(H)

    PonA = "1111 .111 ..11 ..11"
    H = get_subgroup(G, PonA)
    print("*--.--* =", len(H), n//len(H))
    subgroups.append(H)

    LonA = "1111 1111 ..11 ...1"
    H = get_subgroup(G, LonA)
    print(".--*--* =", len(H), n//len(H))
    subgroups.append(H)

    FLAG = "1111 .111 ..11 ...1"
    H = get_subgroup(G, FLAG)
    print("*--*--* =", len(H), n//len(H))
    subgroups.append(H)

#    G = get_permrep(G)
#    for H in subgroups:
#        H = get_permrep(H)
#        X = G.action_subgroup(H) # XX far too slow...
#        print(X)
    

def normal_form(A, p=DEFAULT_P):
    "reduced row-echelon form"
    assert p==2
    #print("normal_form")
    #print(A)
    A = row_reduce(A)
    #print(A)
    m, n = A.shape
    j = 0
    for i in range(m):
        while A[i, j] == 0:
            j += 1
        i0 = i-1
        while i0>=0:
            if A[i0, j]:
                A[i0, :] += A[i, :]
                A %= p
            i0 -= 1
        j += 1
    #print(A)
    return A

CHECK = argv.get("check", False)

def _qchoose_2(n, m, p=DEFAULT_P):
    assert m<=n
    col = m
    row = n-m
    for A in geometry.get_cell(row, col, p):
        yield A


class Figure(object):
    "A partial flag"
    def __init__(self, items, p=DEFAULT_P, normal=False, check=CHECK):
        self.p = p
        if check:
            self.items = items
            self.check()
        if not normal:
            items = [normal_form(A, p) for A in items]
        items = [A.copy() for A in items]
        self.items = items
        #self._str = self.__str__()
        self._str = b' '.join(A.tostring() for A in items)
        if check:
            self.check()

    def check(self):
        items = self.items
        if len(items)<2:
            return
        for idx in range(len(items)-1):
            u = items[idx]
            v = items[idx+1]
            v = set(str(x) for x in span(v))
            #print(v)
            for x in span(u):
                #print("\t", x)
                assert str(x) in v, self

    def __str__(self):
        s = ",".join(str(item) for item in self.items)
        s = s.replace("\n", "")
        return "Figure(%s, %s)"%(s, self.p)
    __repr__ = __str__

    def __hash__(self):
        return hash(self._str)

    def __eq__(self, other):
        return self._str == other._str

    def __ne__(self, other):
        return self._str != other._str

    def __lt__(self, other):
#        print(self._str == other._str)
#        print(self)
#        print(other)
#        assert (self._str == other._str) == (str(self) == str(other))
        return self._str < other._str

    def __le__(self, other):
        return self._str <= other._str

    def __add__(self, other):
        assert isinstance(other, Figure)
        items = self.items + other.items
        return Figure(items, self.p)

    def __rmul__(self, g):
        assert isinstance(g, Matrix)
        A = g.A
        #print("__rmul__")
        #print("\t", self)
        #print(g)
        p = self.p
        items = [(numpy.dot(B, A))%p for B in self.items]
        #print("items:", items)
        return Figure(items, p)
        
    @classmethod
    def get_atomic(cls, m, n, p=DEFAULT_P):
        v = numpy.zeros((m, n), dtype=scalar)
        for i in range(m):
            v[i, i] = 1
        fig = Figure([v], p)
        return fig

    @classmethod
    def qchoose(cls, dims, p=DEFAULT_P):
        assert len(dims)>1
        d0 = dims[0]
        items = []
        for d in dims[1:]:
            assert d<=d0
            items.append(list(_qchoose_2(d0, d)))
            d0 = d
        n = len(items)
        for select in cross(items):
            A = select[0]
            flag = [A]
            for i in range(n-1):
                A = numpy.dot(select[i+1], A) % p
                flag.append(A)
            flag = list(reversed(flag))
            flag = cls(flag)
            yield flag


class Orbit(object):
    def __init__(self, figures):
        figures = list(figures)
        figures.sort()
        self.figures = tuple(figures)

    def __eq__(self, other):
        return self.figures == other.figures

    def __ne__(self, other):
        return self.figures != other.figures

    def __hash__(self):
        return hash(self.figures)



def test():

    figures = set()
    for fig in Figure.qchoose([4,3,1]):
        assert fig not in figures
        figures.add(fig)
    assert len(figures) == 105

    n = argv.get("n", 3)
    G = GL(n)
    print("|G| =", len(G))

    if n==3:
        items = [3, 2, 1]
    elif n==4:
        items = [4, 3, 2, 1]

    orbits = set()
    figures = list(Figure.qchoose(items))
    for p in figures:
      for q in figures:
        fig = p+q
        orbit = set()
        for g in G:
            orbit.add(g*fig)
        orbit = Orbit(orbit)
        orbits.add(orbit)
    print("orbits:", len(orbits))


def test_bruhat():
    n = argv.get("n", 4)

    if n==4:
        point = Figure.get_atomic(1, n)
        line = Figure.get_atomic(2, n)
        plane = Figure.get_atomic(3, n)
    
        figures = [
            point,
            line,
            plane,
            point + line,
            line + plane,
            point + plane,
            point + line + plane,
        ]

    elif n==3:
        point = Figure.get_atomic(1, n)
        line = Figure.get_atomic(2, n)
        figures = [point, line, point + line]

    else:
        return

    G = GL(n)

    Hs = []
    Xs = []
    repG = get_permrep(G)
    #for H in rep.subgroups():
    #    print(len(H), len(G)//len(H))
    for fig in figures:
        #fig = Figure([line, plane])
        orbit = set()
        for g in G:
            u = g*fig
            orbit.add(u)
        print(len(orbit))

        X = list(orbit)
        lookup = dict((v, idx) for (idx, v) in enumerate(X))
        perms = []
        for g in G:
            perm = [lookup[g*v] for v in X]
            #print(perm)
            perms.append(Perm(perm))
        tgt = Group(perms)
        #send_perms = [tgt.lookup[perm] for perm in perms]
        #assert send_perms == list(range(len(send_perms)))
        send_perms = list(range(len(perms)))
        X = GSet(repG, tgt, send_perms)
        Xs.append(X)

        #print(X.send_perms)
        H = X.get_stabilizer()
        Hs.append(H)
        assert len(G)//len(H) == len(orbit)

    if 0:
        found = set(Hs)
        n = len(G)
        print("Hs:", [len(G)//len(H) for H in found])
        for H in list(Hs):
            for g in repG:
                K = Group([g*h*~g for h in H])
                found.add(K)
        print("Hs:", [len(G)//len(H) for H in found])
        while 1:
            pairs = [(H, K) for H in found for K in found]
            for (H, K) in pairs:
                HK = H.intersect(K)
                if len(HK)>1 and HK not in found:
                    found.add(HK)
                    break
            else:
                break
            continue
        Hs = list(found)
        Hs.sort(key = len)
        print("Hs:", [len(G)//len(H) for H in Hs])
        return

    for X in Xs:
        print(X.signature(Hs))
        for Y in Xs:
            XY = X*Y
            Bs = [hom.src for hom in XY.get_atoms()]
            print('\t', [B.signature(Hs) for B in Bs])
    



def main():

    n = argv.get("n", 2)
    p = argv.get("p", 2)

    print("|GL(%d, %d)| = %d"%(p, n, order_gl(n, p)))
    print("|SL(%d, %d)| = %d"%(p, n, order_sl(n, p)))
    if n%2 == 0:
        print("|Sp(%d, %d)| = %d"%(p, n, order_sp(n, p)))

    if argv.SL:
        G = SL(n, p)
        print("|G| =", len(G))

    if argv.GL:
        G = GL(n, p)
        print("|G| =", len(G))



if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()


    

    




