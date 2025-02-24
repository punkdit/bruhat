#!/usr/bin/env python

"""

"""

from random import choice
from operator import mul, matmul, add
from functools import reduce
#from functools import cache
from functools import lru_cache
cache = lru_cache(maxsize=None)

import numpy

from sage.all_cmdline import (FiniteField, CyclotomicField, latex, block_diagonal_matrix,
    PolynomialRing, ZZ)
from sage import all_cmdline 

from bruhat.argv import argv
from bruhat.matrix_sage import Matrix
from bruhat.gset import Group, mulclose


def GL(n,p):
    from bruhat.algebraic import Algebraic, get_permrep
    G = Algebraic.GL(n,p)
    G = get_permrep(G)
    return G


class Rep:
    def __init__(self, ring, G, rep, dim):
        assert isinstance(G, Group)
        self.ring = ring
        self.G = G
        self.rep = rep
        self.dim = dim

    def __str__(self):
        return "Rep(%s, group of order %s, dim=%s)"%(self.ring, len(self.G), self.dim)

    def check(self):
        dim = self.dim
        G = self.G
        for g in G:
          rep_g = self(g)
          assert rep_g.shape == (dim,dim)
          for h in G:
            rep_h = self(h)
            rep_gh = self(g*h)
            assert rep_g*rep_h == rep_gh

    @classmethod
    def regular(cls, ring, G):
        n = len(G)
        lookup = {g:i for (i,g) in enumerate(G)}
        rep = {}
        for i,g in enumerate(G):
            rows = []
            for j,h in enumerate(G):
                k = lookup[g*h]
                row = [0]*n
                row[k] = 1
                rows.append(row)
            M = Matrix(ring, rows)
            rep[g] = M.t
        return Rep(ring, G, rep, n)

    @classmethod
    def trivial(cls, ring, G):
        M = Matrix(ring, [[1]])
        rep = {g:M for g in G}
        return Rep(ring, G, rep, 1)

    def __call__(self, g):
        return self.rep[g]

    def __add__(self, other):
        assert self.G is other.G
        assert self.ring is other.ring
        rep = {g : self(g).direct_sum(other(g)) for g in self.G}
        return Rep(self.ring, self.G, rep, self.dim+other.dim)

    def __matmul__(self, other):
        assert self.G is other.G
        assert self.ring is other.ring
        rep = {g : self(g)@other(g) for g in self.G}
        return Rep(self.ring, self.G, rep, self.dim*other.dim)

    def restrict(self, H):
        rep = {h:self.rep[h] for h in H}
        return Rep(self.ring, H, rep, self.dim)
    

def test_repr():
    print("test()")

    ring = ZZ

    G = GL(3,2)
    n = len(G)
    rep = Rep.regular(ZZ, G)
    print(rep)

    rep.check()




def test():
    print("test()")

    ring = ZZ

    G = GL(3,2)
    n = len(G)
    rep = Rep.regular(ZZ, G)
    print(rep)

    I = G.identity
    for g in G:
        if g != I and g*g*g==I:
            break
    H = Group.generate([g])
    print(len(H))

    res = rep.restrict(H)
    print(res)

    res.check()

    v = Rep.trivial(ring, G)
    v.check()

    vv = v@v
    vv.check()
    
    vv = v+v
    vv.check()

    



if __name__ == "__main__":

    from time import time
    start_time = time()

    profile = argv.profile
    name = argv.next() or "test"
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
    print("OK! finished in %.3f seconds\n"%t)




