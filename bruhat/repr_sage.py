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

from sage.all_cmdline import (
    FiniteField, CyclotomicField, latex, block_diagonal_matrix,
    PolynomialRing, ZZ, QQ)
from sage import all_cmdline 

from bruhat.argv import argv
from bruhat.matrix_sage import Matrix
from bruhat.action import mulclose, mulclose_hom
from bruhat.gset import Perm, Group, Coset
from bruhat.smap import SMap


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

    def dot(other, self): # other <---- self
        assert self.G is other.G
        #u = ring.zero()
        u = 0
        for i in range(self.n):
            u += self.chi[i].conjugate() * other.chi[i]
        u /= self.n
        return u


class Rep:
    ring = CyclotomicField()

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
            M = Matrix(cls.ring, rows)
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
            M = Matrix(cls.ring, M)
            M = M.t
            rep[g] = M
        return Rep(G, rep, dim)

    @classmethod
    def trivial(cls, G):
        M = Matrix(cls.ring, [[1]])
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
        u = cls.ring(u)
        rep = {}
        for j in range(n):
            M = Matrix(cls.ring, [[u**(j*i)]])
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

    def is_irrep(self):
        chi = self.chi
        r = chi.dot(chi)
        return r==1

    def restrict(self, H):
        rep = {h:self.rep[h] for h in H}
        return Rep(H, rep, self.dim)

    def induce(self, G):
        ring = self.ring
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
    
        Ms = []
        for g in G.get_gens():
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
            #rep[g] = M
            Ms.append(M)
        return Rep.generate(G, G.get_gens(), Ms)

    def hom(w, v): 
        """
            find basis for space of intertwiners
            w <--- v
        """
        assert isinstance(v, Rep)
        assert isinstance(w, Rep)
        assert v.G is w.G
    
        G = v.G
        Iv = Matrix.identity(v.ring, v.dim)
        Iw = Matrix.identity(w.ring, w.dim)
        blocks = []
        #for g in G:
        G.get_gens()
        for g in G.gens:
            lhs = v(~g) @ Iw
            rhs = Iv @ w(g)
            f = lhs-rhs
            blocks.append(f)
        M = rowcat(blocks)
        K = M.cokernel()
        homs = []
        for f in K:
            f = f.reshape(v.dim, w.dim).t
            hom = Hom(w, v, f)
            homs.append(hom)
        return homs
    

class Hom:
    "intertwiner of Rep's"
    def __init__(self, tgt, src, M):
        assert isinstance(tgt, Rep)
        assert isinstance(src, Rep)
        assert isinstance(M, Matrix)
        assert tgt.G is src.G
        self.tgt = tgt
        self.src = src
        self.G = tgt.G
        assert M.shape == (tgt.dim, src.dim)
        self.M = M

    def __str__(self):
        return "(%s<---%s)"%(self.tgt, self.src)

    def check(self):
        tgt = self.tgt
        src = self.src
        M = self.M
        for g in self.G:
            assert M*src(g) == tgt(g)*M

    def __mul__(lhs, rhs):
        assert isinstance(rhs, Hom)
        assert lhs.src == rhs.tgt
        M = lhs.M * rhs.M
        return Hom(lhs.tgt, rhs.src, M)

    def cokernel(self):
        tgt = self.tgt
        M = self.M
        G = self.G
        K = M.cokernel()
        Ki = K.pseudoinverse()
        rep = {}
        for g in G:
            rep[g] = K * tgt(g) * Ki
        dim = K.shape[0]
        r = Rep(G, rep, dim)
        return Hom(r, tgt, K)



# https://ncatlab.org/nlab/show/Gram-Schmidt+process#CategorifiedGramSchmidtProcess
class GramSchmidt:
    def __init__(self, reps):
        self.reps = reps

    def showtable(self):
        reps = self.reps
        N = len(reps)
        print("   ", "===="*N)
        chis = [r.chi for r in reps]
        for i in range(N):
          print("%3d:"%i, end="")
          for j in range(N):
            u = chis[j].dot(chis[i])
            print("%3s "%u, end="", flush=True)
          print()
        print("   ", "===="*N)
        print()

    def subtract(self, i, j):
        reps = self.reps
        print("subtract", i, j)
        fs = reps[i].hom(reps[j])
        assert len(fs)
        print("\thoms:", len(fs))
        f = fs[0]
        reps[i] = f.cokernel().tgt

    def __getitem__(self, idx):
        return self.reps[idx]



def test_rep():
    print("test_rep()")

    A = Matrix(Rep.ring, [[1,1,1],[1,0,1],[0,0,1]])
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

    G = Group.symmetric(4)
    H = Group([g for g in G if g[3] == 3]) # S3 in G

    a = Perm([0,2,1,3])
    b = Perm([1,2,0,3])
    assert a in H
    assert b in H

    # 2d rep of S3
    A = Matrix(Rep.ring, [[0,1],[1,0]])
    B = Matrix(Rep.ring, [[0,-1],[1,-1]])
    rep = Rep.generate(H, [a,b], [A,B])
    rep.check()
    #rep.dump()

    r = rep.induce(G)
    assert r.dim == 8
    r.check()


def test_gram_schmidt():
    #Rep.ring = QQ # slower than CyclotomicField() !

    G = Group.symmetric(4)
    reps = []
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
    #for i in range(N):
    #  for j in range(N):
    #    v = reps[i]
    #    w = reps[j]
    #    homs = w.hom(v)
    #    if v.dim*w.dim > 100:
    #        continue
    #    print(len(homs), end=" ", flush=True)
    #    for f in homs:
    #        f.check()
    #  print()

    fs = reps[2].hom(reps[1])
    for f in fs:
        f.check()
    f = fs[0]
    k = f.cokernel()
    k.tgt.check()
    k.check()

    gs = GramSchmidt(reps)
    showtable = gs.showtable
    subtract = gs.subtract
    showtable()

    subtract(1, 0)
    subtract(2, 0)
    subtract(3, 0)
    subtract(4, 0)
    showtable()

    subtract(2, 1)
    subtract(3, 1)
    subtract(3, 1)
    subtract(4, 1)
    subtract(4, 1)
    subtract(4, 1)
    showtable()

    subtract(3, 2)
    subtract(4, 2)
    subtract(4, 2)
    showtable()

    subtract(4, 3)
    subtract(4, 3)
    subtract(4, 3)
    showtable()

    for r in gs.reps:
        print(r)



def GL(n,p):
    from bruhat import algebraic
    from bruhat.algebraic import Algebraic, get_subgroup, parse

    G = Algebraic.GL(n,p)
    #for g in G.gen:
    #    print(g)

    torus = []
    for M in [
        algebraic.Matrix([[0,0,1],[1,0,0],[0,1,1]]),
        algebraic.Matrix([[0,0,1],[1,0,1],[0,1,0]]),
    ]:
        #H = Algebraic([M])
        #assert len(H) == 7
        H = [g for g in G if g*M==M*g]
        H = get_permrep(H)
        torus.append(H)

    parabolics = []
    assert n==3

    FLAG = "111 .11 ..1"
    FLAG = parse(FLAG).reshape(n,n)
    flag = get_subgroup(G, FLAG)
    print(len(flag), len(G)//len(flag))
    parabolics.append(flag)

    point = "111 .11 .11"
    point = parse(point).reshape(n,n)
    H = get_subgroup(G, point)
    #print("*--. =", len(H), len(G)//len(H))
    parabolics.append(H)

    LINE = "111 111 ..1"
    LINE = parse(LINE).reshape(n,n)
    H = get_subgroup(G, LINE)
    #print(".--* =", len(H), len(G)//len(H))
    parabolics.append(H)

    child = "1.. .11 .11"
    child = parse(child).reshape(n,n)
    child = get_subgroup(G, child)
    child = get_permrep(child)

    parabolics = [get_permrep(H) for H in parabolics]
    G = get_permrep(G)
    parabolics.append(G)

    for H in parabolics:
        assert G.is_subgroup(H)
        #print(H)
    G.parabolics = parabolics
    G.torus = torus
    G.child = child

    return G


def get_permrep(G):
    """
    permutation action of G on the non-zero vectors (in lexicographic order)
    """
    from bruhat import gset
    from bruhat.algebraic import Matrix, enum2
    assert len(G)
    op = G[0]
    n, _ = op.shape
    p = op.p
    space = [Matrix(numpy.array(v).reshape(n,1)) for v in enum2(n) if sum(v)!=0]
    lookup = dict((v, idx) for (idx, v) in enumerate(space))
    perms = []
    rep = {}
    for g in G:
        idxs = [lookup[g*v] for v in space]
        perm = gset.Perm(idxs)
        rep[g] = perm
        perms.append(perm)
    X = gset.Group(perms)
    X.get_gens()
    X.rep = rep
    X.space = space
    X.G = G
    return X


def test_parabolic_induction():
    from bruhat import algebraic
    from bruhat.algebraic import Algebraic, get_subgroup, parse

    n, p = 3, 2
    GLn = Algebraic.GL(n,p)
    X = get_permrep(GLn)
    r = get_parabolic_induction(X)

    print(r)
    r.check()
    print(r.is_irrep())



def sort_sign(items):
    sign = +1
    N = len(items)
    for n in range(N-1, 0, -1): # bubble sort
        swapped = False
        for i in range(n):
            if items[i] > items[i + 1]:
                items[i], items[i + 1] = items[i + 1], items[i]
                swapped = True
                sign *= -1
        if not swapped:
            break
    jtems = list(items)
    jtems.sort()
    assert jtems==items
    return sign
    

def get_parabolic_induction(X):

    points = X.space
    lines = set()
    for a in points:
      for b in points:
        if a==b:
            continue
        c = a+b
        l = [a,b,c]
        l.sort()
        l = tuple(l)
        lines.add(l)
    lines = list(lines)
    lines.sort()

    rep = {}
    for g in X.G:
        #print(g)
        rows = []
        for i,l in enumerate(lines):
            a, b, c = l
            l1 = [g*a, g*b, g*c]
            l2 = list(l1)
            #l2.sort()
            sign = sort_sign(l2)
            l2 = tuple(l2)
            j = lines.index(l2)
            #print(i, "-->", j, sign)
            row = [0]*len(lines)
            row[j] = sign
            rows.append(row)
        M = Matrix(Rep.ring, rows)
        #print(M.t)
        perm = X.rep[g]
        rep[perm] = M.t

    r = Rep(X, rep, len(lines))
    return r




def test_gl():
    Rep.ring = QQ
    #Rep.ring = CyclotomicField(7)

    G = GL(3,2)
    n = len(G)

    print("test_gl()")

    print(G)

    H = G.child
    print(H)
    rep = {}
    I = G.identity
    for g in H:
        if g != I and g*g == I:
            rep[g] = Matrix(Rep.ring, [[-1]])
        else:
            rep[g] = Matrix(Rep.ring, [[1]])
    r = Rep(H, rep, 1)
    r.check()
    r = r.induce(G)
    print(r)
    #return

    Hs = G.parabolics
    reps = [Rep.permutation(G, H) for H in [Hs[0], Hs[1], Hs[3]]]

    reps.append(get_parabolic_induction(G))
    #reps.append(r)

    for r in reps:
        print(r)
    gs = GramSchmidt(reps)
    gs.showtable()

    gs.subtract(0, 2)
    gs.subtract(1, 2)
    gs.showtable()
    gs.subtract(0, 1)
    gs.showtable()
    gs.subtract(0, 1)
    gs.showtable()
#    gs.subtract(4, 0)
#    gs.showtable()
#    gs.subtract(4, 3)
#    gs.showtable()
#    gs.subtract(4, 3)
#    gs.showtable()

    for r in gs.reps:
        print(r)
        print("is_irrep:", r.is_irrep())

    #return

    #reps.append(Rep.permutation(G, G.torus[0]))
    Rep.ring = CyclotomicField(7)

    for H in G.torus:
        for j in [1,3]: #range(1,7):
            r = Rep.fourier(H, j)
            print(r, "induce:")
            r_torus = r.induce(G)
            print("\t", r_torus)
            gs.reps.append(r_torus)
            #break
        break
    #return

    gs.showtable()

    gs.subtract(4, 0)
    gs.showtable()
    gs.subtract(4, 3)
    gs.showtable()
    #gs.subtract(4, 1)
    #gs.showtable()

    for r in reps:
        print(r)
    return

    reps = [Rep.permutation(G, H) for H in G.parabolics[1:]]
    for r in reps:
        print(r)
    gs = GramSchmidt(reps)

    subtract = gs.subtract
    showtable = gs.showtable

    showtable()
    subtract(0,2)
    subtract(1,2)
    #subtract(0, 3)
    #subtract(1, 3)
    #subtract(2, 3)

    showtable()

    r0 = (gs.reps[0])
    r1 = (gs.reps[1])

    assert r0.is_irrep()

    fs = r0.hom(r1)
    assert len(fs)==1
    f = fs[0]
    print(f)
    print(f.M)


def test():
    test_rep()
    test_gram_schmidt()
    test_gl()


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




