#!/usr/bin/env python

"""
practice representation_theory, induce / restrict, etc.

good notes here:
https://dec41.user.srcf.net/h/II_L/representation_theory/10

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
from bruhat.action import mulclose, mulclose_hom
from bruhat.gset import Perm, Group, Coset

ring = CyclotomicField()

def GL(n,p):
    from bruhat.algebraic import Algebraic, get_permrep
    G = Algebraic.GL(n,p)
    G = get_permrep(G)
    return G


def colcat(col):
    M = None
    for block in col:
        M = block if M is None else M.stack(block)
    return M

def rowcat(row):
    M = None
    for block in row:
        M = block if M is None else M.augment(block)
    return M


class Char:
    def __init__(self, G, chi):
        assert len(chi) == len(G)
        self.n = len(chi)
        self.chi = chi
        self.G = G

    def __str__(self):
        return "Char(%s)"%(self.chi,)

    def __eq__(self, other):
        assert self.G is other.G
        return self.chi == other.chi

    def __getitem__(self, g):
        i = self.G.lookup[g]
        return self.chi[i]

    def __add__(self, other):
        assert self.G is other.G
        chi = [self.chi[i]+other.chi[i] for i in range(self.n)]
        return Char(self.G, chi)

    def __sub__(self, other):
        assert self.G is other.G
        chi = [self.chi[i]-other.chi[i] for i in range(self.n)]
        return Char(self.G, chi)

    def __mul__(self, other):
        assert self.G is other.G
        chi = [self.chi[i]*other.chi[i] for i in range(self.n)]
        return Char(self.G, chi)

    def hom(self, other):
        assert self.G is other.G
        u = ring.zero()
        for i in range(self.n):
            u += self.chi[i].conjugate() * other.chi[i]
        u /= self.n
        return u


class Rep:
    def __init__(self, G, rep, dim):
        assert isinstance(G, Group)
        self.G = G
        self.rep = rep
        self.dim = dim
        self.chi = Char(self.G, [rep[g].trace() for g in G])

    def __str__(self):
        return "Rep(group of order %s, dim=%s)"%(len(self.G), self.dim)

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

    def dump(self):
        for g in self.G:
            print(g)
            print(self(g))
    
    def __eq__(self, other):
        "equality on-the-nose"
        assert self.G is other.G
        if self.dim != other.dim:
            return False
        return self.rep == other.rep

    @classmethod
    def regular(cls, G):
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
        return Rep(G, rep, n)

    @classmethod
    def permutation(cls, G, H):
        X = G.action_subgroup(H)
        dim = X.rank
        rep = {}
        for i,g in enumerate(G):
            j = X.send_perms[i]
            h = X.tgt[j]
            #print(h)
            M = [[0]*dim for _ in range(dim)]
            for (src,tgt) in enumerate(h):
                M[src][tgt] = 1
            M = Matrix(ring, M)
            M = M.t
            rep[g] = M
        return Rep(G, rep, dim)

    @classmethod
    def trivial(cls, G):
        M = Matrix(ring, [[1]])
        rep = {g:M for g in G}
        return Rep(G, rep, 1)

    @classmethod
    def generate(cls, G, gs, Ms):
        assert len(gs) == len(Ms)
        assert len(gs)
        rep = mulclose_hom(gs, Ms)
        assert len(rep) == len(G)
        dim = len(Ms[0])
        return Rep(G, rep, dim)

    @classmethod
    def fourier(cls, A, i=0):
        assert A.is_abelian()
        for a in A:
            # find a generator
            if len(mulclose([a])) == len(A):
                break
        else:
            assert 0, "not implemented" # XXX
        n = len(A)
        u = CyclotomicField(n).gen()
        u = ring(u)
        rep = {}
        for j in range(n):
            M = Matrix(ring, [[u**(j*i)]])
            rep[a**j] = M
        return Rep(A, rep, 1)

    def __call__(self, g):
        return self.rep[g]

    def __add__(self, other):
        assert self.G is other.G
        rep = {g : self(g).direct_sum(other(g)) for g in self.G}
        return Rep(self.G, rep, self.dim+other.dim)

    def __matmul__(self, other):
        assert self.G is other.G
        rep = {g : self(g)@other(g) for g in self.G}
        return Rep(self.G, rep, self.dim*other.dim)

    def dual(self):
        rep = {g : self.rep[~g].t for g in self.G}
        return Rep(self.G, rep, self.dim)

    def restrict(self, H):
        rep = {h:self.rep[h] for h in H}
        return Rep(H, rep, self.dim)

    def induce(self, G):
        H = self.G
        dim = self.dim
        cosets = G.left_cosets(H)
        n = len(cosets)
        #print("left_cosets:", n)
        ts = [gH[0] for gH in cosets]
        #print(ts)
        for t,gH in zip(ts, cosets):
            assert Coset([t*h for h in H]) == gH
    
        H = set(H)
    
        rep = {}
        for g in G:
            #print(g)
            lookup = {}
            for i in range(n):
              for j in range(n):
                if (~ts[j])*g*ts[i] in H:
                    assert lookup.get(i) is None
                    lookup[i] = j
            #print(lookup)
            cols = []
            for i in range(n):
                #blocks = [Matrix.zero(ring,dim,dim) for _ in range(n)]
                blocks = [Matrix.zero(ring,dim,dim)]*n
                j = lookup[i]
                U = self((~ts[j])*g*ts[i])
                #print(U, U.shape)
                blocks[j] = U
                cols.append(colcat(blocks))
            M = rowcat(cols)
            #print(M)
            rep[g] = M
        return Rep(G, rep, dim*n)


def hom(v, w):
    assert isinstance(v, Rep)
    assert isinstance(w, Rep)
    assert v.G is w.G




def test_repr():
    print("test()")

    G = GL(3,2)
    n = len(G)
    rep = Rep.regular(G)
    print(rep)

    rep.check()




def test():
    print("test()")

    A = Matrix(ring, [[1,1,1],[1,0,1],[0,0,1]])
    B = A.t
    C = A.solve(B)
    assert C is not None
    assert A*C == B

    #G = GL(3,2)
    G = Group.symmetric(3)

    n = len(G)
    reg = Rep.regular(G)
    reg.check()

    I = G.identity
    for g in G:
        if g != I and g*g*g==I:
            break
    H = Group.generate([g])
    assert len(H) == 3

    res = reg.restrict(H)
    res.check()

    v = Rep.trivial(G)
    v.check()

    vv = v@v
    vv.check()
    
    vv = v+v
    vv.check()

    v = Rep.trivial(H)
    r = Rep.fourier(H, 1)
    r.check()
    assert r==r
    assert r@r != r
    assert r@r@r == v

    assert Rep.fourier(H, 0) == v
    assert Rep.fourier(H, 2) == r@r

    s = r.induce(G)
    s.check()

    #reg.dump()
    #print("permutation:")
    for H in G.subgroups():
        r = Rep.permutation(G, H)
        r.check()
        if len(H) == 1:
            assert r.chi == reg.chi
        #    r.dump()
            print(r.chi)

    G = Group.symmetric(4)
    H = Group([g for g in G if g[3] == 3]) # S3 in G

    a = Perm([0,2,1,3])
    b = Perm([1,2,0,3])
    assert a in H
    assert b in H

    # 2d rep of S3
    A = Matrix(ring, [[0,1],[1,0]])
    B = Matrix(ring, [[0,-1],[1,-1]])
    rep = Rep.generate(H, [a,b], [A,B])
    rep.check()
    #rep.dump()

    r = rep.induce(G)
    assert r.dim == 8
    r.check()

#    for g in G:
#        print(g.order())

    # -------------------------------

    # https://ncatlab.org/nlab/show/Gram-Schmidt+process#CategorifiedGramSchmidtProcess

    reps = []
    #for H in G.subgroups():
    U4 = G
    U31 = Group([g for g in G if g[3]==3])
    U22 = Group([g for g in G if g[0] in [0,1] and g[1] in [0,1]])
    U211 = Group([g for g in G if g[3]==3 and g[2]==2])
    U1111 = Group([g for g in G if g[3]==3 and g[2]==2 and g[1]==1])

    Hs = [U4, U31, U22, U211, U1111]
    for H in Hs:
        r = Rep.permutation(G, H)
        reps.append(r)
        v = Rep.trivial(H)
        v = v.induce(G)
        assert v==r # too strong ?
        assert v.chi == r.chi
        #v.dump()
    N = len(reps)

    for i in range(N):
      for j in range(N):
        v = reps[i]
        w = reps[j]
        hom(v,w)

    return

    chis = [r.chi for r in reps]

    for i in range(N):
      for j in range(N):
        u = chis[i].hom(chis[j])
        #dim = (reps[i].dual() @ reps[j]).dim
        #print("%s:%s"%(u,dim), end=" ", flush=True)
        print(u, end=" ", flush=True)
      print()

#    print()
#    for i in range(N):
#      for j in range(N):
#        dim = (reps[i].dual() @ reps[j]).dim
#        print(dim, end=" ", flush=True)
#      print()



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




