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

from sage.all_cmdline import FiniteField, CyclotomicField, latex, ZZ, QQ
from sage import all_cmdline 

from bruhat.argv import argv
from bruhat.matrix_sage import Matrix
from bruhat.action import mulclose, mulclose_hom
from bruhat.gset import Perm, Group, Coset
from bruhat.smap import SMap
from bruhat import gset
from bruhat.algebraic import Matrix as FMatrix


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
        G = self.G
        items = []
        for cls in G.conjugacy_classes():
            i = G.lookup[cls[0]]
            items.append(self.chi[i])
        return "Char(%s)"%(items,)

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
        return "Rep(group of order %s, dim=%s, irrep=%s)"%(len(self.G), self.dim, self.is_irrep())

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
    def sign(cls, G):
        M = Matrix(cls.ring, [[1]])
        rep = {g:g.sign()*M for g in G}
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

    def hom_fail(w, v): 
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
            #lhs = v(~g) @ Iw
            #rhs = Iv @ w(g)

            #lhs = v(g) @ Iw
            #rhs = Iv @ w(~g)

            #lhs = Iw @ v(~g).t
            #rhs = w(g) @ Iv

            lhs = Iw @ v(g)
            rhs = w(~g).t @ Iv

            f = lhs-rhs
            blocks.append(f.t)
        M = rowcat(blocks)
        K = M.cokernel()
        assert (K*M).is_zero()
        print("Rep.hom()")
        print("\tK:", K)
        homs = []
        for f in K:
            #f = f.reshape(w.dim, v.dim) #??
            f = f.reshape(v.dim, w.dim).t
            hom = Hom(w, v, f)
            #hom.other = Hom(w,v,f0)
            homs.append(hom)
        return homs
    
    def hom(w, v): 
        """
            find basis for space of intertwiners
            w <--- v
        """
        assert isinstance(v, Rep)
        assert isinstance(w, Rep)
        assert v.G is w.G

        #print("hom", w, "<--", v)
        m, n = w.dim, v.dim
    
        G = v.G
        Iv = Matrix.identity(v.ring, v.dim)
        Iw = Matrix.identity(w.ring, w.dim)
        blocks = []
        #for g in G:
        G.get_gens()

        from bruhat.system import Unknown, System, array, dot
        f = Unknown(m, n)
        system = System(v.ring, f)
        for g in G.gens:
            wg = (w(g)).to_numpy()
            vg = (v(g)).to_numpy()
            lhs = dot(wg, f)
            rhs = dot(f, vg)
            system.append(lhs-rhs, 0)

        #f = system.solve()
        
        homs = []
        #for f in system.all_solutions():
        for f in system.solve_homogeneous():
            #print("f =")
            #print(f)
            f = f.reshape(w.dim, v.dim)
            hom = Hom(w, v, f)
            hom.check()
            #hom.other = Hom(w,v,f0)
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
        assert self.is_valid()

    def is_valid(self):
        tgt = self.tgt
        src = self.src
        M = self.M
        for g in self.G.gens:
            if M*src(g) != tgt(g)*M:
                return False
        return True

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
class Basis:
    def __init__(self, reps):
        self.reps = reps

    def __getitem__(self, idx):
        return self.reps[idx]

    def __len__(self):
        return len(self.reps)

    def append(self, rep):
        assert isinstance(rep, Rep)
        self.reps.append(rep)

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
          print(" ", reps[i], reps[i].chi)
        print("   ", "===="*N)
        print()

    def subtract(self, i, j, idx=0):
        reps = self.reps
        print("subtract: %d-%d" %( i, j ), end=", ")
        fs = reps[i].hom(reps[j])
        assert len(fs)
        print("homs:", len(fs))
        f = fs[idx]
        reps[i] = f.cokernel().tgt



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
    Rep.ring = QQ # ARGHHH much slower than CyclotomicField() !

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

    basis = Basis(reps)
    showtable = basis.showtable
    subtract = basis.subtract
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

    for r in basis.reps:
        print(r)


def test_induce():

    print("\n\ntest_induce()")

    n = 3
    G = Group.symmetric(n)
    print(G)

    for g in G:
        if g.order() == n:
            break
    H = Group.generate([g])

#    rep = Rep.fourier(H, 1)
#    rep = rep.induce(G)
#    print(rep)
#    print(rep.chi)
#    rep.dump()
#
#    fs = rep.hom(rep)
#    print(len(fs))
#    print(rep.is_irrep())

    #return
    Rep.ring = QQ

    r0 = Rep.trivial(G)
    r1 = Rep.sign(G)

    reg = Rep.regular(G)

    for f in reg.hom(r0):
        f.check()
    for f in reg.hom(r1):
        f.check()

    a = Perm([0,2,1])
    b = Perm([1,2,0])
    assert a in G
    assert b in G

    # 2d rep of S3
    A = Matrix(Rep.ring, [[0,1],[1,0]])
    B = Matrix(Rep.ring, [[0,-1],[1,-1]])
    r2 = Rep.generate(G, [a,b], [A,B])
    r2.check()
    #rep.dump()

    homs = reg.hom(r2)
    assert len(homs) == 2
    for f in homs:
        print("is_valid:", f.is_valid()) #, f.other.is_valid())

    k = homs[0].cokernel()
    k.check()
    r = k.tgt
    r.check()

    print(r)
    print("\n\nDEBUG:")

    homs = r.hom(r2)
    assert len(homs) == 1
    for f in homs:
        f.src.check()
        f.tgt.check()
        print("src:", f.src)
        #f.src.dump()
        print("tgt:", f.tgt)
        #f.tgt.dump()
        print("f:")
        print(f.M)
        print("is_valid:", f.is_valid()) #, f.other.is_valid())
        f.check()


class Space:
    "vector space over finite field F_p"
    def __init__(self, n, p):
        self.n = n
        self.p = p
        shape = (p,)*n
        space = [v for v in numpy.ndindex(shape) if sum(v)] # skip zero vector ?
        assert len(space) == p**n - 1
        space = [FMatrix(numpy.array(v).reshape(n,1), p) for v in space]
        lookup = dict((v, idx) for (idx, v) in enumerate(space))
        self.space = space
        self.lookup = lookup

    def __str__(self):
        s = str(self.space)
        s = s.replace(" ", "").replace("\n", "")
        return s

    def __getitem__(self, idx):
        return self.space[idx]

    def __len__(self):
        return len(self.space)

    def get_permrep(self, G):
        """
        permutation action of G on the non-zero vectors (in lexicographic order)
        """
        assert len(G)
        space, lookup = self.space, self.lookup
        op = G[0]
        n, _ = op.shape
        p = op.p
        perms = []
        rep = {}
        #print(lookup)
        #print("lookup:")
        #for k,v in lookup.items():
        #    print(k)
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


def GL32():
    from bruhat import algebraic
    from bruhat.algebraic import Algebraic, get_subgroup, parse

    n, p = 3, 2
    G = Algebraic.GL(n,p)
    #for g in G.gen:
    #    print(g)

    space = Space(n, p)
    get_permrep = space.get_permrep

    torus = []
    for M in [
        FMatrix([[0,0,1],[1,0,0],[0,1,1]]),
        FMatrix([[0,0,1],[1,0,1],[0,1,0]]),
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


def test_gl2():
    from bruhat import algebraic
    from bruhat.algebraic import Algebraic, get_subgroup, parse

    p = 3
    GL1 = Algebraic.GL(1,p)
    GL2 = Algebraic.GL(2,p)

    s1 = Space(1, p)
    s2 = Space(2, p)

    X1 = s1.get_permrep(GL1)
    X2 = s2.get_permrep(GL2)
    print(X1, s1)
    print(X2, s2)
    #return

    point = parse("11 .1").reshape(2,2)
    point = get_subgroup(GL2, point)
    point = s2.get_permrep(point)

    print(point)

    r0 = Rep.trivial(X2)
    r1 = Rep.permutation(X2, point)

    basis = Basis([r0, r1])
    basis.showtable()
    basis.subtract(1,0)
    basis.showtable()


    dim = len(s2)
    rep = {}
    for g in GL2:
        #print(g)
        cols = []
        for v in s2:
            col = [0]*dim
            gv = g*v
            idx = s2.lookup[gv]
            col[idx] = 1
            cols.append(col)
        M = Matrix(Rep.ring, cols).t
        rep[X2.rep[g]] = M
    rep = Rep(X2, rep, dim)
    print(rep)
    #rep.check()

    basis.append(rep)
    basis.subtract(2,0)
    basis.subtract(2,1)
    basis.showtable()

    #basis[2].dump()

    remain = set(s2)
    orbits = []
    lookup = {}
    sign = {}
    while remain:
        v = remain.pop()
        u = 2*v
        assert u in remain
        remain.remove(u)
        o = (u,v)
        orbits.append(o)
        lookup[u] = o
        lookup[v] = o
        sign[u] = +1
        sign[v] = -1

    rep = {}
    for g in GL2:
        s = 1
        for (u,v) in orbits:
            assert lookup[g*u] == lookup[g*v]
            u1 = g*u
            s *= sign[u]*sign[u1]
        x = s*Matrix(Rep.ring, [[1]])
        rep[X2.rep[g]] = x
        #print(s, end=" ")
    #print()
    rep = Rep(X2, rep, 1)
    print(rep)
    rep.check()

    basis.append(rep)
    basis.showtable()


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


def test_parabolic_induction():
    from bruhat import algebraic
    from bruhat.algebraic import Algebraic, get_subgroup, parse

    G = GL32()
    r = get_parabolic_induction(G)

    print(r)
    r.check()
    assert r.is_irrep()





def test_gl():
    Rep.ring = QQ
    #Rep.ring = CyclotomicField(7)

    G = GL32()
    n = len(G)

    clss = G.conjugacy_classes()
    assert len(clss) == 6

    print("test_gl()")

    print(G)

    #r = Rep.permutation(G, Group([G.identity]))
    rep = {}
    for g in G:
        M = Matrix.get_perm(Rep.ring, g)
        rep[g] = M
    cnot = Rep(G, rep, G.rank)

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
    basis = Basis(reps)
    basis.showtable()

    basis.subtract(0, 2)
    basis.subtract(1, 2)
    basis.showtable()
    basis.subtract(0, 1)
    basis.showtable()
    basis.subtract(0, 1)
    basis.showtable()
#    basis.subtract(4, 0)
#    basis.showtable()
#    basis.subtract(4, 3)
#    basis.showtable()
#    basis.subtract(4, 3)
#    basis.showtable()

    for r in basis.reps:
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
            basis.reps.append(r_torus)
            #break
        break
    for i,r in enumerate(basis.reps):
        print(i, r)
        #print(r.chi)

    chis = [rep.chi for rep in basis.reps]
    chis[4] = chis[4] - chis[0] - chis[1] - chis[3]
    chis[5] = chis[5] - chis[0] - chis[1] - chis[3]

    for c in chis:
      for d in chis:
        print("%5s"%(c.dot(d)), end=" ")
      print()
    print()

    for c in chis:
        print(c)

    basis.showtable()
    basis.subtract(4,0)
    basis.showtable()
    basis.subtract(4,1)
    basis.showtable()
    basis.subtract(4,3)
    basis.showtable()

    basis.subtract(5,0)
    basis.showtable()
    basis.subtract(5,1)
    basis.showtable()
    basis.subtract(5,3)
    basis.showtable()

    for r in basis.reps:
        print(r)

    #basis.reps.append(cnot)
    #basis.showtable()


def test_structure():

    from bruhat.algebraic import qchoose_2

    n = 4
    p = 3

    for m in range(n+1):
        count = 0
        for item in qchoose_2(n, m, p):
            count += 1
        print(count)

    print( p**0 + p + p**2 + p**3 )
    print( p**0 + p + 2*p**2 + p**3 + p**4 )



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




