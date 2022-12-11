#!/usr/bin/env python
"""
q polynomials for grassmanians, etc.

See also: coxeter.py

"""

import string, os
from time import sleep, time
from functools import reduce
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import matmul

import numpy
from numpy import alltrue, zeros, dot

from bruhat.argv import argv
from bruhat.smap import SMap
from bruhat.element import Z
from bruhat.poly import Poly
from bruhat.todd_coxeter import Schreier

# See:
# https://en.wikipedia.org/wiki/Q-analog
# https://math.ucr.edu/home/baez/week187.html


ring = Z
zero = Poly({}, ring)
one = Poly({():1}, ring)
q = Poly("q", ring)



class Coxeter(object):
    """
        Coxeter reflection group.
    """
    def __init__(self, ngen, rel):
        """
        rel : map pairs (i,j) of generators to m_{ij}
            (this is the Coxeter matrix)
        """
        for (i, j) in rel.keys():
            assert 0<=i<ngen
            assert 0<=j<ngen
        gen = list(range(ngen))
        for i in gen:
          for j in gen:
            if rel.get((i, j)) and not rel.get((j, i)):
                rel[j, i] = rel[i, j]
        for i in gen:
          for j in gen:
            m = rel.get((i, j))
            if m is None:
                rel[i, j] = 2 # default
        for i in gen:
          for j in gen:
            assert rel[i, j] == rel[j, i]
            #assert rel[i, j] in (2, 3, 4, 6)
        self.ngen = ngen
        self.rel = rel
        key = list(rel.items())
        key.sort()
        self.key = tuple(key)

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    @cache
    def get_graph(self):
        return Schreier.make_reflection(self.ngen, self.rel)

    @cache
    def get_group(self):
        graph = self.get_graph()
        G = graph.get_group()
        return G

    def get_residual_graph(self, i):
        assert 0 <= i < self.ngen
        graph = Schreier.make_reflection(self.ngen, self.rel, False)
        hgens = [(j,) for j in range(self.ngen) if i!=j]
        graph.build(hgens)
        return graph

    def get_order(self):
        return len(self.get_group())

    @cache
    def get_poincare(self, ring, q):
        "the poincare polynomial"
        graph = self.get_graph()
        return graph.get_poincare(ring, q)

    def slow_get_poincare(self, ring, q):
        G = self.get_group()
        gen = G.gen
        I = G.identity
        found = {I:0}
        bdy = set([I])
        dist = 0
        while bdy:
            # bdy is the last thing we found
            dist += 1
            _bdy = set()
            for h in bdy:
              for g in gen:
                gh = g*h
                if gh in found:
                    continue
                found[gh] = dist
                _bdy.add(gh)
            bdy = _bdy
        found = list(found.values())
        p = ring.zero
        for i in range(dist):
            p = p + found.count(i) * q**i
        return p

    def subgroup(self, gen):
        lookup = dict((g, i) for (i,g) in enumerate(gen))
        rel = {}
        for i in gen:
          for j in gen:
            rel[lookup[i], lookup[j]] = self.rel[i, j]
        return Coxeter(len(gen), rel)

    def residual(self, i):
        assert 0<=i<self.ngen
        gen = [j for j in range(self.ngen) if j!=i]
        return self.subgroup(gen)

    def __mul__(self, other):
        ngen = self.ngen + other.ngen
        #lookup = dict((i, ngen+i))
        rel = dict(self.rel)
        for key,val in other.rel.items():
            i, j = key
            rel[i+self.ngen, j+self.ngen] = val
        return Coxeter(ngen, rel)

    @classmethod
    @cache
    def A(cls, n):
        rel = {(i, i+1):3 for i in range(n-1)}
        return cls(n, rel)

    @classmethod
    @cache
    def B(cls, n):
        rel = {}
        for i in range(n-1):
            rel[i, i+1] = 4 if i==n-2 else 3
        return cls(n, rel)
    C = B

    @classmethod
    @cache
    def D(cls, n):
        rel = {}
        for i in range(n-2):
            rel[i, i+1] = 3
        if n>2:
            rel[n-3, n-1] = 3
        return cls(n, rel)

    @classmethod
    @cache
    def E(cls, n):
        rel = {}
        for i in range(1, n-1):
            rel[i, i+1] = 3
        rel[0, 3] = 3
        return cls(n, rel)

    @classmethod
    @cache
    def F(cls):
        n = 4
        rel = {(0, 1):3, (1, 2):4, (2, 3): 3}
        return cls(n, rel)

    @classmethod
    @cache
    def G(cls):
        return cls(2, {(0,1):6})


A0 = Coxeter.A(0)
A1 = Coxeter.A(1)
A2 = Coxeter.A(2)
A3 = Coxeter.A(3)
A4 = Coxeter.A(4)

assert A0 is Coxeter.A(0)

G = A4.get_group()
assert len(G) == 120
assert G is A4.get_group()

G = (A2*A3).get_group()
assert len(G) == 6*24

assert A2 == A2
assert A2 != A3
assert A4.subgroup([1, 2]) == A2
assert A4.subgroup([0, 2, 3]) == A1 * A2
assert A4.residual(2) == A2*A1


B0 = Coxeter.B(0)
B1 = Coxeter.B(1)
B2 = Coxeter.B(2)
B3 = Coxeter.B(3)
B4 = Coxeter.B(4)

assert len(B2.get_group()) == 8
assert len(B3.get_group()) == 48
assert len(B4.get_group()) == 384

assert B3.residual(0) == B2
assert B3.residual(1) == A1*A1
assert B3.residual(2) == A2

D0 = Coxeter.D(0)
D1 = Coxeter.D(1)
D2 = Coxeter.D(2)
D3 = Coxeter.D(3)
D4 = Coxeter.D(4)
D5 = Coxeter.D(5)
D6 = Coxeter.D(6)
assert len(D4.get_group()) == 192
#assert D3 == A3 # needs re-indexing
assert D2 == A1*A1
assert D4.residual(1) == A1*A1*A1

E6 = Coxeter.E(6)
E7 = Coxeter.E(7)
E8 = Coxeter.E(8)
F4 = Coxeter.F()
G2 = Coxeter.G()


def shortstr(p):
    assert isinstance(p, Poly)
    keys = p.keys
    items = {}
    m = 0
    for key in keys:
        if len(key)==0:
            items[0] = p.cs[key]
        else:
            c, n = key[0]
            items[n] = p.cs[key]
            m = max(n, m)
    items = [items.get(i, 0) for i in range(m+1)]
    if max(items)>9:
        return "(%s)"%','.join(str(c) for c in items)
    else:
        return ''.join(str(c) for c in items)
    


class Poincare(object):
    def __init__(self, ring, q):
        self.ring = ring
        self.zero = ring.zero
        self.one = ring.one
        self.q = q

    def qbracket(self, n):
        q = self.q
        p = self.zero
        for i in range(n):
            p = p + q**i
        return p
    
    def qfactorial(self, n):
        q = self.q
        p = one
        for i in range(n):
            p = p * self.qbracket(i+1)
        return p
    
    def A(self, n):
        return self.qfactorial(n+1)
    
    def B(self, n):
        p = one
        for i in range(n):
            p = p * self.qbracket(2*(i+1))
        return p
    C = B
    
    def D(self, n):
        p = one
        if n==1:
            return 1+self.q # arghh...
        for i in range(n-1):
            p = p * self.qbracket(2*(i+1))
        if n>0: 
            p = p * self.qbracket(n)
        return p


def test():

    print("test()")

    poincare = Poincare(ring, q)
    A, B, D = poincare.A, poincare.B, poincare.D

    assert q**0 == one

    p = (1+q)**5

    assert poincare.qbracket(5) == 1 + q + q**2 + q**3 + q**4

    assert poincare.qfactorial(3) == 1 + 2*q + 2*q**2 + q**3

    # These count points of the maximal flag variety over the
    # field with q elements:
    assert shortstr(A(0)) == "1"
    assert shortstr(A(1)) == "11"
    assert shortstr(A(2)) == "1221"
    assert shortstr(A(3)) == "1356531"

    assert A(3) == A3.get_poincare(ring, q)

    assert shortstr(B(1)) == "11"
    assert shortstr(B(2)) == "12221"
    assert shortstr(B(3)) == "1357887531"

    assert B(3) == B3.get_poincare(ring, q)

    assert shortstr(D(1)) == "11"
    assert shortstr(D(2)) == "121"
    assert shortstr(D(3)) == "1356531"

    assert D(1) == D1.get_poincare(ring, q)
    assert D(2) == D2.get_poincare(ring, q)
    assert D(3) == D3.get_poincare(ring, q)
    assert D(4) == D4.get_poincare(ring, q)

    p = D(4)
    p = p / (A(1)**3)
    assert shortstr(p) == "1133443311"

    #print(p)

    get = lambda p : p.substitute((('q',2),))
    assert get(p) == 1575

    assert get(B(1) / A(0)) == 3

    assert get(B(2) / A(1)) == 15

    assert get(B(3) / B(2)) == 63
    assert get(B(3) / (A(1)**2)) == 315
    assert get(B(3) / A(2)) == 135

    assert get(B(4) / B(3)) == 255
    assert get(B(4) / (A(1) * B(2))) == 5355
    assert get(B(4) / (A(2) * A(1))) == 11475
    assert get(B(4) / A(3)) == 2295


if __name__ == "__main__":

    start_time = time()


    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")


