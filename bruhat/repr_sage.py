#!/usr/bin/env python

"""
practice representation_theory, induce / restrict, etc.

good notes here:
https://dec41.user.srcf.net/h/II_L/representation_theory/10

see also: bruhat.cuspforms 

https://ncatlab.org/nlab/show/Gram-Schmidt+process#CategorifiedGramSchmidtProcess
"""

from math import lcm
from random import choice, shuffle
from string import ascii_uppercase, ascii_lowercase
from operator import mul, matmul, add
from functools import reduce
#from functools import cache
from functools import lru_cache
cache = lru_cache(maxsize=None)

ascii_names = ascii_uppercase + ascii_lowercase

import numpy

from sage.all_cmdline import FiniteField, CyclotomicField, latex, ZZ, QQ
from sage import all_cmdline as sage
sage.Parallelism().set(nproc=1)

from bruhat.argv import argv
from bruhat.matrix_sage import Matrix
from bruhat.action import mulclose, mulclose_hom
from bruhat.gset import Perm, Group, Coset
from bruhat.smap import SMap
from bruhat import gset
from bruhat.algebraic import Algebraic
from bruhat.algebraic import Matrix as FMatrix
from bruhat.util import cross, all_primes


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
        self._chi = chi
        self.G = G
        self.ring = Rep.ring
        self.name = None

    def __str__(self):
        G = self.G
        items = []
        for cls in G.conjugacy_classes():
            i = G.lookup[cls[0]]
            items.append(self._chi[i])
        return "Char(%s)"%(items,)

    def __eq__(self, other):
        assert self.G is other.G
        return self._chi == other._chi

    def __getitem__(self, g):
        i = self.G.lookup[g]
        return self._chi[i]

    def __add__(self, other):
        assert self.G is other.G
        chi = [self._chi[i]+other._chi[i] for i in range(self.n)]
        return Char(self.G, chi)

    def __sub__(self, other):
        assert self.G is other.G
        chi = [self._chi[i]-other._chi[i] for i in range(self.n)]
        return Char(self.G, chi)

    def __mul__(self, other):
        assert self.G is other.G
        chi = [self._chi[i]*other._chi[i] for i in range(self.n)]
        return Char(self.G, chi)

    def dot(other, self): # other <---- self
        assert self.G is other.G
        #u = ring.zero()
        u = 0
        for i in range(self.n):
            u += self._chi[i].conjugate() * other._chi[i]
        u /= self.n
        return u


class Rep:
    #ring = CyclotomicField()
    ring = QQ
    name = None

    def __init__(self, G, rep, dim):
        assert isinstance(G, Group)
        self.G = G
        self.rep = rep
        self.dim = dim
        self.chi = Char(self.G, [rep[g].trace() for g in G])

    def __str__(self):
        extra = ", name=%r"%self.name if self.name else ""
        return "Rep(group of order %s, dim=%s, irrep=%s%s)"%(
            len(self.G), self.dim, self.is_irrep(), extra)
    __repr__ = __str__

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

    def some_check(self, count=100):
        dim = self.dim
        G = self.G
        gs = list(G)
        for i in range(count):
            g = choice(gs)
            rep_g = self(g)
            assert rep_g.shape == (dim,dim)
            h = choice(gs)
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

    def is_iso(self, other):
        if self.chi != other.chi:
            return False
        for f in self.hom(other):
            if f.is_iso():
                return True
        return False # ???

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
        ring = CyclotomicField(n)
        u = ring.gen()
        #u = cls.ring(u)
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

    def __mul__(self, other):
        G = self.G
        H = other.G
        GH = []
        rep = {}
        for g in G:
          for h in H:
            gh = g*h
            assert gh not in rep
            rep[gh] = self(g) @ other(h)
            GH.append(gh)
        GH = Group(GH, G.gens + H.gens)
        return Rep(GH, rep, self.dim*other.dim)

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

    @classmethod
    def irreps(cls, G):
        Hs = G.conjugacy_subgroups()
        reps = [Rep.permutation(G, H) for H in Hs]
        basis = Basis(reps)
        basis.reduce()
        return basis
    

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

    def __eq__(self, other):
        assert self.src is other.src
        assert self.tgt is other.tgt
        return self.M == other.M

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

    def __add__(lhs, rhs):
        assert isinstance(rhs, Hom)
        assert lhs.src == rhs.src
        assert lhs.tgt == rhs.tgt
        M = lhs.M + rhs.M
        return Hom(lhs.tgt, lhs.src, M)

    def __sub__(lhs, rhs):
        assert isinstance(rhs, Hom)
        assert lhs.src == rhs.src
        assert lhs.tgt == rhs.tgt
        M = lhs.M - rhs.M
        return Hom(lhs.tgt, lhs.src, M)

    def __rmul__(lhs, u):
        M = u*lhs.M
        return Hom(lhs.tgt, lhs.src, M)

    def __invert__(self):
        M = self.M.pseudoinverse()
        return Hom(self.src, self.tgt, M)

    def is_iso(self):
        tgt, src = self.tgt, self.src
        if tgt.dim != src.dim:
            return False
        iself = ~self
        a = self*iself
        b = iself*self
        I = Matrix.identity(Rep.ring, src.dim)
        return a.M==I and b.M==I

    def cokernel(self):
        tgt = self.tgt
        M = self.M
        G = self.G
        K = M.cokernel()
        if not len(K):
            return None # ???
        Ki = K.pseudoinverse()
        #print(K, K.shape)
        #print(Ki, Ki.shape)
        rep = {}
        for g in G:
            rep[g] = K * tgt(g) * Ki
        dim = K.shape[0]
        r = Rep(G, rep, dim)
        return Hom(r, tgt, K)

    def kernel(self):
        src = self.src
        M = self.M
        G = self.G
        K = M.kernel()
        if not len(K):
            return None # ???
        Ki = K.pseudoinverse()
        #print(K, K.shape)
        #print(Ki, Ki.shape)
        rep = {}
        for g in G:
            rep[g] = Ki * src(g) * K
        dim = K.shape[1]
        r = Rep(G, rep, dim)
        return Hom(src, r, K)



class Table:
    "Character Table"
    def __init__(self, reps, verbose=True):
        chis = []
        for rep in reps:
            chi = rep.chi
            chi.name = rep.name
            chis.append(chi)
        self.chis = chis
        self.verbose = verbose

    def __getitem__(self, idx):
        return self.chis[idx]

    def __setitem__(self, idx, value):
        if isinstance(value, Rep):
            chi = value.chi
            chi.name = value.name
            value = chi
        self.chis[idx] = value

    def __delitem__(self, idx):
        del self.chis[idx]

    def __len__(self):
        return len(self.chis)

    def append(self, rep):
        assert isinstance(rep, Rep)
        chi = rep.chi
        chi.name = rep.name
        self.chis.append(chi)

    def get(self, i, j):
        chis = self.chis
        c, d = chis[i], chis[j]
        return d.dot(c)

    def show(self, fmt="%s", showreps=True, showchars=False):
        chis = self.chis
        N = len(chis)
        print("   ", "="*N)
        for i in range(N):
          print("%3d:"%i, end="")
          for j in range(N):
            u = chis[j].dot(chis[i])
            print(fmt%("." if u==0 else u), end="", flush=True)
          if showchars:
            print(" ", chis[i])
          else:
            print()
        print("   ", "="*N)
        print("rank:", self.rank())
        print()

    def subtract(self, i, j, idx=0):
        chis = self.chis
        chis[i] = chis[i] - chis[j]

    def swap(self, i, j):
        chis = self.chis
        #print("Basis.swap", i, j)
        chis[i], chis[j] = chis[j], chis[i]

    def pop(self, i):
        #print("Basis.pop", i)
        return self.chis.pop(i)

    def reduce(self, verbose=False):
        row = 0
        while row < len(self):
            if verbose:
                print("reduce", row)
                self.show()
            if self.get(row, row) != 1:
                #break # ??
                #for sow in range(row+1, len(self)):
                for sow in range(len(self)):
                    if sow==row:
                        continue
                    u = self.get(row, sow)
                    if u==1:
                        if verbose:
                            print("\tsubtract", row, sow)
                        self.subtract(row, sow)
                        break
                row += 1
            else:
                for sow in range(row+1, len(self)):
                    u = self.get(row, sow)
                    if self.get(sow, sow) == 1 and u:
                        assert u == 1, u
                        if verbose:
                            print("\tpop", sow)
                        self.pop(sow)
                        break
                    if u:
                        if verbose:
                            print("\tsubtract", row, sow)
                        self.subtract(sow, row)
                        break
                else:
                    row += 1
        
    def rank(self):
        if not len(self):
            return 0
        
        r = self[0]
        G = r.G
        rows = [[chi[g] for g in G] for chi in self]
        M = Matrix(r.ring, rows)
        #print(M)
        return M.rank()



# https://ncatlab.org/nlab/show/Gram-Schmidt+process#CategorifiedGramSchmidtProcess
class Basis:
    "categorified Character Table"
    def __init__(self, reps, verbose=True):
        self.reps = reps
        self.verbose = verbose

    def __getitem__(self, idx):
        return self.reps[idx]

    def __setitem__(self, idx, value):
        self.reps[idx] = value

    def __delitem__(self, idx):
        del self.reps[idx]

    def __len__(self):
        return len(self.reps)

    def append(self, rep):
        assert isinstance(rep, Rep)
        self.reps.append(rep)

    def get(self, i, j):
        reps = self.reps
        c, d = reps[i].chi, reps[j].chi
        return d.dot(c)

    def show(self, fmt="%s", showreps=True, showchars=False):
        reps = self.reps
        N = len(reps)
        print("   ", "="*N)
        chis = [r.chi for r in reps]
        for i in range(N):
          print("%3d:"%i, end="")
          for j in range(N):
            u = chis[j].dot(chis[i])
            print(fmt%("." if u==0 else u), end="", flush=True)
          if showreps:
            print(" ", reps[i], end="")
          if showchars:
            print(" ", reps[i].chi)
          else:
            print()
        print("   ", "="*N)
        print("rank:", self.rank())
        print()

    def subtract(self, i, j, idx=0):
        reps = self.reps
        #print("Basis.subtract: %d-%d" %( i, j ), end=", ")
        fs = reps[i].hom(reps[j])
        assert len(fs)
        #print("homs:", len(fs))
        f = fs[idx]
        reps[i] = f.cokernel().tgt

    def swap(self, i, j):
        #print("Basis.swap", i, j)
        reps = self.reps
        reps[i], reps[j] = reps[j], reps[i]

    def pop(self, i):
        #print("Basis.pop", i)
        return self.reps.pop(i)

    def reduce(self, verbose=False):
        row = 0
        while row < len(self):
            #print("Basis.reduce", row)
            #self.show()
            if self.get(row, row) != 1:
                #break # ??
                #for sow in range(row+1, len(self)):
                for sow in range(len(self)):
                    if sow==row:
                        continue
                    u = self.get(row, sow)
                    if u==1:
                        self.subtract(row, sow)
                        break
                row += 1
            else:
                for sow in range(row+1, len(self)):
                    u = self.get(row, sow)
                    if self.get(sow, sow) == 1 and u:
                        assert u == 1, u
                        self.pop(sow)
                        break
                    if u:
                        self.subtract(sow, row)
                        break
                else:
                    row += 1
        
    def rank(self):
        if not len(self):
            return 0
        
        r = self[0]
        G = r.G
        rows = [[r.chi[g] for g in G] for r in self]
        M = Matrix(r.ring, rows)
        #print(M)
        return M.rank()


class Algebra:
    def __init__(self, fs):
        self.fs = fs
        rows = []
        for f in fs:
            M = f.M
            n = M.shape[0]*M.shape[1]
            M = M.reshape(1,n)
            rows.append(M.M[0])
        self.A = Matrix(Rep.ring, rows).t
    def getname(self, M):
        A = self.A
        M = M.reshape(M.shape[0]*M.shape[1], 1)
        a = A.solve(M).t
        return a
    def __str__(self):
        return "Algebra(dim=%d)"%(len(self.fs),)
    def display(self):
        print(self)
        fs = self.fs
        for f in fs:
          for g in fs:
            gf = g*f
            #if gf in fs:
            #    print(fs.index(gf), end=" ")
            #else:
            #    print("?", end=" ")
            a = self.getname(gf.M)
            print(a, end=" ")
          print()

    
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
        hom = {}
        for g in G:
            idxs = [lookup[g*v] for v in space]
            perm = gset.Perm(idxs)
            hom[g] = perm
            perms.append(perm)
        X = gset.Group(perms)
        X.get_gens()
        X.hom = hom
        X.space = space
        X.G = G
        return X

    def get_parabolic(self, G, figure):
        n = self.n
        from bruhat.algebraic import parse, get_subgroup
        figure = parse(figure).reshape(n,n)
        H = get_subgroup(G, figure)
        H = self.get_permrep(H)
        return H



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

    assert rep.is_iso(rep)

    U = Matrix(Rep.ring, [[1,1],[0,1]])
    Ui = ~U
    A = U*A*Ui
    B = U*B*Ui
    sep = Rep.generate(H, [a,b], [A,B])
    sep.check()

    assert rep != sep
    assert rep.is_iso(sep)
    



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

#    basis = Basis(reps)
    basis = Table(reps)
    show = basis.show
    subtract = basis.subtract
    show()

    subtract(1, 0)
    subtract(2, 0)
    subtract(3, 0)
    subtract(4, 0)
    show()

    subtract(2, 1)
    subtract(3, 1)
    subtract(3, 1)
    subtract(4, 1)
    subtract(4, 1)
    subtract(4, 1)
    show()

    subtract(3, 2)
    subtract(4, 2)
    subtract(4, 2)
    show()

    subtract(4, 3)
    subtract(4, 3)
    subtract(4, 3)
    show()

    for r in basis:
        print(r)

    assert basis.rank() == 5


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
    basis.show()
    basis.subtract(1,0)
    basis.show()


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
    basis.show()

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
    basis.show()


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
        perm = X.hom[g]
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
    print(cnot)
    return

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
    basis.show()

    basis.subtract(0, 2)
    basis.subtract(1, 2)
    basis.show()
    basis.subtract(0, 1)
    basis.show()
    basis.subtract(0, 1)
    basis.show()
#    basis.subtract(4, 0)
#    basis.show()
#    basis.subtract(4, 3)
#    basis.show()
#    basis.subtract(4, 3)
#    basis.show()

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

    basis.show()
    basis.subtract(4,0)
    basis.show()
    basis.subtract(4,1)
    basis.show()
    basis.subtract(4,3)
    basis.show()

    basis.subtract(5,0)
    basis.show()
    basis.subtract(5,1)
    basis.show()
    basis.subtract(5,3)
    basis.show()

    for r in basis.reps:
        print(r)

    if argv.cnot:
        print("cnot:")
        basis.reps.append(cnot)
        basis.show()


def test_structure():

    p = 2
    n = 4
    m = 2
    k = n-m
    # n == k+n
    items = list(FMatrix.qchoose(n, m, p))
    assert len(items) == p**0 + p + 2*p**2 + p**3 + p**4

    for item in items:
        assert item.normal_form() == item
    lookup = {item:i for (i,item) in enumerate(items)}

    GL = Algebraic.GL(n, p)
    print(len(GL), len(items))
    dim = len(items)

    Rep.ring = ZZ
    GLm = Algebraic.GL(m, p)
    print(GLm)
    Xm = Space(2, p).get_permrep(GLm)
    r = Rep.sign(Xm)
    print(r, r.chi)

    space = Space(n, p)
    Xn = space.get_permrep(GL)

    ring = ZZ

    rep = {}
    for g in GL.gens:
        cols = []
        for idx,item in enumerate(items):
            jtem = (g*item.t).t
            ktem = jtem.normal_form()
            h = ktem.t.solve(jtem.t)
            assert h in GLm
            u = r(Xm.hom[h])
            jdx = lookup[ktem]
            col = [0]*dim
            col[jdx] = u.M[0,0]
            cols.append(col)
        rep[Xn.hom[g]] = Matrix(ring, cols).t
    gs = list(rep.keys())
    rep = Rep.generate(Xn, gs, [rep[g] for g in gs])
    #rep = Rep(G, rep, dim)
    print(rep, rep.chi)
    rep.some_check()
            

def get_2d(G):
    assert isinstance(G, Group)
    assert len(G) == 6

    a = Perm([0,2,1])
    b = Perm([1,2,0])

    for a in G:
      for b in G:
        # 2d rep of S3
        A = Matrix(Rep.ring, [[0,1],[1,0]])
        B = Matrix(Rep.ring, [[0,-1],[1,-1]])
        hom = mulclose_hom([a,b], [A,B])
        if len(hom) != len(G):
            continue
        rep = Rep.generate(G, [a,b], [A,B])
        try:
            rep.check()
            return rep
        except AssertionError:
            pass
    assert 0


def test_levi():

    n, p = 3, 3
    space = Space(n, p)
    GL = Algebraic.GL(n,p)

    levis, uni, parabolic = GL.levi_decompose([1,2])
    assert len(levis) == 2
    assert len(levis[0]) == 2
    assert len(levis[1]) == 48
    assert len(uni) == 3**2, len(uni)
    assert len(parabolic) == 2*48 * 3**2, len(parabolic)

    n, p = 4, 2
    space = Space(n, p)
    GL = Algebraic.GL(n,p)

    levis, uni, parabolic = GL.levi_decompose([1,3])
    assert len(levis) == 2
    assert len(levis[0]) == 1
    assert len(levis[1]) == 168
    assert len(uni) == 8
    assert len(parabolic) == 1344, len(parabolic)

    levis, uni, parabolic = GL.levi_decompose([2,2])
    assert len(levis) == 2
    assert len(levis[0]) == 6
    assert len(levis[1]) == 6
    assert len(uni) == 16
    assert len(parabolic) == 576 

    #parabolic = GL.get_parabolic([0,1,0])

    G = space.get_permrep(GL)
    P = Group([G.hom[g] for g in parabolic]) # == get_permrep(parabolic)

    #l0 = GL.get_levi(0,2)
    #l1 = GL.get_levi(2,2)
    l0, l1 = levis
    GL2 = Algebraic.GL(2,p)
    for g in l0:
        h = l0.hom[g]
        assert g in GL, str(h)
        assert h in GL2
        assert g in parabolic
    for g in GL2:
        assert l0.ihom[g] in l0

    #uni = GL.get_unipotent()
    #assert len(uni) == p**6
    #uni = [g for g in uni if g[0,1]==g[2,3]==0]
    #assert len(uni) == p**4
    Uni = space.get_permrep(uni)
    Uni.do_check()
    assert len(Uni) == p**4

    L0 = space.get_permrep(l0)
    assert P.is_subgroup(L0)
    reps0 = Rep.irreps(L0)
    assert len(reps0) == 3

    L1 = space.get_permrep(l1)
    assert P.is_subgroup(L1)
    reps1 = Rep.irreps(L1)

    print( P )
    print( Uni )
    print( L0 )
    print( L1 )

    L = set(l0*l1 for l0 in L0 for l1 in L1)
    L = Group(L)
    assert P.is_subgroup(L)
    assert P.is_subgroup(Uni)
    assert P.is_normal(Uni)
    assert len(L) * len(Uni) == len(P)

    basis = []
    for r0 in reps0:
      for r1 in reps1:
        r01 = r0*r1
        r01.check()
        print(r01)
        L = r01.G
        #print( len(L) * len(Uni) , len(P) )
        rep = {}
        for l in L:
          for u in Uni:
            lu = u*l
            assert lu not in rep
            rep[lu] = r01(l)
        rep = Rep(P, rep, r01.dim)
        print(rep)
        rep.check()
        rep = rep.induce(G)
        print(rep)
        basis.append(rep)
      print()

    #basis = Basis(basis)
    #basis.show()


def test_permutation():

    G = Group.symmetric(3)
    r = Rep.regular(G)
    print(r)

    homs = r.hom(r)
    print(len(homs))

    #for hom in homs:
    #    print(hom.M)

    perms = []
    for f in homs:
      perm = []
      for g in homs:
        gf = g*f
        i = homs.index(gf)
        perm.append(i)
      perms.append(Perm(perm))
    E = Group.generate(perms)
    I = E.identity
    refls = [a for a in E if a*a==I and a!=I]
    a, b = refls[:2]
    assert a*b*a == b*a*b

    for H in G.subgroups():
        if len(H) == 2:
            break
    r = Rep.permutation(G, H)
    print(r)
    homs = r.hom(r)
    print(len(homs))
    homs = [f.M for f in homs]
    for f in homs:
        print(f)
    f = homs[1]
    print(f*f)


def get_torus(gl):
    from bruhat.cuspforms import slow_get_conjugacy_classes
    cgys = slow_get_conjugacy_classes(gl)

    print(len(cgys))

    F = sage.GF(p)
    for cgy in cgys:
        g = cgy[0]
        A = sage.matrix(F, n, n, g.A)
        char = A.characteristic_polynomial()
        #char.factor()
        if not char.is_irreducible():
            continue
        assert p**2-1 == 24
        print([g.order() for g in cgy])
        #for g in cgy:
        #    if g.order() == p**2-1:
        #        break
        #else:
        #    assert 0
        #print(g, char.factor(), g.order())


def test_irr():

    #n = argv.get("n", 5)
    #G = Group.symmetric(n)
    #G = Group.dihedral(n)
    #G = Group.alternating(n)

#    chis = burnside_irr(G)
#    print()
#    for chi in chis:
#        print(chi)
#
#    #return

    n = argv.get("n", 2)
    p = argv.get("p", 7)
    space = Space(n, p)
    gl = Algebraic.GL(n,p)
    G = space.get_permrep(gl)
    print(G)

    chis = dixon_irr(G)
    for chi in chis:
        print(chi)

    assert len(G) == sum(chi[0]**2 for chi in chis)
    cgys = G.conjugacy_classes()
    N = len(cgys)
    for i in range(N):
      for j in range(N):
        x = reduce(add, [len(cgys[k]) * chis[i][k] * chis[j][k].conjugate() for k in range(N)])
        if i==j:
            assert x == len(G)
        else:
            assert x == 0




def dixon_irr(G):
    #print(G)

    m = 1
    for g in G:
        k = g.order()
        m = lcm(m, k)
    #print("m =", m)

    for p in all_primes(10*len(G)):
        if p <= 2*(len(G)**0.5):
            continue
        if (p-1)%m == 0:
            break
    else:
        assert 0, "whoops, ran out of primes, call Euler"
    #print("p =", p)
    assert p % len(G) != 0

    ring = FiniteField(p)
    one = ring.one()
    zero = ring.zero()

    # choose z having multiplicative order m when viewed as an element of Z*_p.
    for z in range(1, p):
        for i in range(1, m):
            if (z**i)%p == 1:
                break
        else:
            if (z**m)%p == 1:
                break
        continue
    else:
        assert 0

    z = one*z
    #print("z =", z)
    assert z**m == 1


    K = G.conjugacy_classes()
    N = len(K)
    print("dixon_irr: cgys =", N)
    #for cgy in K:
    #    print(len(cgy), end=" ")
    #print()

    lookup = {}
    for i,cgy in enumerate(K):
        for g in cgy:
            lookup[g] = i

    rev = {}
    for i in range(N):
      for g in K[i]:
        rev[g] = i
    inv = []
    for cgy in K:
        j = None
        for g in cgy:
            assert j is None or j == rev[~g]
            j = rev[~g]
        inv.append(j)

    if 0:
        cmats = [numpy.zeros((N,N), dtype=int) for r in range(N)]
        idxs = {K[t][0]:t for t in range(N)}
        for r in range(N):
          for s in range(N):
            #print(".", end="", flush=True)
            for x in K[r]:
              for y in K[s]:
                xy = x*y
                t = idxs.get(xy)
                if t is None: # TODO too slow !!
                    continue
                key = (r,s,t)
                #matrix[key] = matrix.get(key, 0) + 1
                cmats[r][s,t] += 1
          #print()

    mats = [numpy.zeros((N,N), dtype=int) for r in range(N)]
    for r in range(N):
      for t in range(N):
        xy = K[t][0] # arbitrary choice
        for x in K[r]:
            y = (~x)*xy # remove inverse here ? HOTSPOT
            s = lookup[y]
            mats[r][s,t] += 1
    #for (mat,cmat) in zip(mats, cmats):
    #    assert numpy.all(mat==cmat)
    cmats = mats

    cmats = [Matrix(ring, cmat) for cmat in cmats]
    for a in cmats:
      for b in cmats:
        assert a*b==b*a

    best = None
    for A in cmats:
        evecs = A.M.eigenvectors_right()
        if best is None:
            best = evecs
        w = max([item[2] for item in best])
        if max([item[2] for item in evecs]) < w:
            best = evecs

    while len(best) < N:
        #print("best:", len(best))
        dims = []
        basis = []
        for val,vecs,dim in best:
            dims.append(dim)
            basis += vecs
        U = Matrix.promote(ring, basis).t
        #U = Matrix(ring, basis).t
        V = U.pseudoinverse()

        i = 0
        for idx,dim in enumerate(dims):
            if dim > 1:
                break
            i += 1
        else:
            assert 0
        
        for A in cmats:
            B = V*A*U
            b = B[i:i+dim, i:i+dim]
            assert b.shape == (dim,dim)
            if b.is_zero():
                continue
            evecs = b.M.eigenvectors_right()
            if len(evecs) == dim:
                break
        else:
            assert 0

        best.pop(idx)
        for (val, vecs, dim) in evecs:
            #print(val, vecs, dim)
            assert len(vecs)==dim==1
            v = vecs[0]
            u = [0]*N
            for row in range(len(v)):
                u[i+row] = v[row]
            u = Matrix(ring, u).t
            u = U*u
            assert A*u == val*u
            u = tuple(u.M[j,0] for j in range(N)) # UGRRRH!
            best.append( (val, [u,], dim) )

    omega = []
    for (val, vecs, dim) in best:
        assert dim==1
        vec = vecs[0]
        omega.append(vec)
        vec = Matrix(ring, vec).t
        assert vec.M[0,0] == 1

        for A in cmats:
            u = A*vec
            val = u.M[0,0]
            assert u == val*vec
        #assert A*vec == val*vec

    root = {}
    for u in range(p//2+1):
        u = u*one
        #print(u, u*u)
        root[u*u] = u
    #return

    tchis = []
    for i in range(N):
        chi1 = zero
        for r in range(N):
            chi1 += omega[i][r] * omega[i][inv[r]] / len(K[r])
        chi1 = 1/chi1
        chi1 = len(G)*chi1
        chi1 = root[chi1]
        chi = [chi1]
        for j in range(1, N):
            c = omega[i][j] * chi1 / len(K[j])
            chi.append(c)
        tchis.append(chi)

    tchis.sort( key = lambda chi : (int(chi[0]),str(chi)) )
    #for chi in tchis:
    #    print(chi)
    assert len(G) == sum(int(chi[0])**2 for chi in tchis)

    F = CyclotomicField()
    zero = F.zero()
    one = F.one()
    zeta = F.gen(m)
    #print(zeta)
    chis = []
    for tchi in tchis:
      d = int(tchi[0]) # _dimension
      chi = []
      for i in range(N):
        cgy = K[i]
        x = cgy[0]
        #k = x.order()
        xls = [G.identity]
        g = x
        while g != G.identity:
            xls.append(g)
            g = x*g
        k = len(xls)
        assert x.order() == k
        eta = zeta ** (m//k)
        #print("eta:", eta)
        value = zero
        for s in range(k):
            #u = reduce(add, [tchi[lookup[x**l]]/(z**(s*l*m//k)) for l in range(k)])
            u = reduce(add, [tchi[lookup[xls[l]]]/(z**(s*l*m//k)) for l in range(k)]) # HOTSPOT
            u = int(u/k)
            value += u * (eta ** s)
        chi.append(value)
      chis.append(chi)

    if len(chis)>1:
        chis.sort( key = lambda chi : (chi[0],str(chi).count("-"),str(chi)) )

    return chis
    
        


def burnside_irr(G):
    #print(G)

    K = G.conjugacy_classes()
    N = len(K)
    #print("cgys:", N)
    #for cgy in K:
    #    print(len(cgy), end=" ")
    #print()

    rev = {}
    for i in range(N):
      for g in K[i]:
        rev[g] = i
    inv = []
    for cgy in K:
        j = None
        for g in cgy:
            assert j is None or j == rev[~g]
            j = rev[~g]
        inv.append(j)

    cmats = [numpy.zeros((N,N), dtype=int) for r in range(N)]
    idxs = {K[t][0]:t for t in range(N)}
    for r in range(N):
     for s in range(N):
        for x in K[r]:
          for y in K[s]:
            z = x*y
            t = idxs.get(z)
            if t is None:
                continue
            key = (r,s,t)
            #matrix[key] = matrix.get(key, 0) + 1
            cmats[r][s,t] += 1

    #ring = QQ
    ring = CyclotomicField()

    cmats = [Matrix(ring, cmat) for cmat in cmats]
    for a in cmats:
      for b in cmats:
        assert a*b==b*a
    for A in cmats:
        #print(A)
        evs = A.M.eigenvalues()
        #print(evs)
        if len(set(evs)) != len(evs):
            continue
        evecs = A.M.eigenvectors_right()
        break
    else:
        print("need simultaneous eigenvectors_right")
        return

    omega = []
    for (val, vecs, dim) in evecs:
        assert dim==1
        vec = vecs[0]
        omega.append(vec)
        vec = Matrix(ring, vec).t
        assert vec.M[0,0] == 1
        #print(val, vec.t)
        assert A*vec == val*vec

    one = ring.one()
    zero = ring.zero()
    chis = []
    for i in range(N):
        chi1 = zero
        for r in range(N):
            chi1 += omega[i][r] * omega[i][inv[r]] / len(K[r])
        chi1 = 1/chi1
        chi1 = len(G)*chi1
        chi1 = chi1**(one/2)
        chi = [chi1]
        for j in range(1, N):
            c = omega[i][j] * chi1 / len(K[j])
            chi.append(c)
        chis.append(chi)

    if len(chis)>1:
        chis.sort( key = lambda chi : (chi[0],str(chi).count("-"),str(chi)) )
        
    return chis
        





    
    
        
    


def test_cuspidal(n=2, p=5):

    n = argv.get("n", n)
    p = argv.get("p", p)

    space = Space(n, p)
    gl = Algebraic.GL(n,p)
    gl1 = Algebraic.GL(1,p)
    print("|GL(%d,%d)| = %d"%(n, p, len(gl)))

    Rep.ring = CyclotomicField(p**n - 1)

    if (n,p) == (2,5):
        ts = [gl.get_torus(i) for i in range(10)]
    
        for a in ts:
          for b in ts:
            c = set(a).intersection(b)
            print(len(c), end=" ")
          print()

    torus = gl.get_torus()
    assert len(torus) == p**n-1

    stab = []
    inv = gl.get_inv()
    for g in gl:
        gi = inv[g]
        for h in torus:
            if g*h*gi not in torus:
                break
        else:
            stab.append(g)
    print(len(stab))

    G = space.get_permrep(gl)
    Torus = space.get_permrep(torus)
    Stab = space.get_permrep(stab)

#    for X in [
#        X = G.action_subgroup(Stab),
#        X = G.action_subgroup(Torus)
#    ]:
#    
#        XX = X*X
#        for o in XX.get_atoms():
#            print(o)
#    
#        r = Rep.permutation(G, Torus)
#        print(r)
#        fs = r.hom(r)
#        print(len(fs))
#        A = Algebra(fs)
#        A.display()

    
    diag = [1]*n
    levis, uni, parabolic = gl.levi_decompose(diag)
    P = space.get_permrep(parabolic)

#    for j in range(p**n-1):
#        r = Rep.fourier(Torus, j)
#        r = r.induce(G)
#        print(r)
#        fs = r.hom(r)
#        A = Algebra(fs)
#        A.display()
#    return


    l0, l1 = levis
    for g in l0:
        h = l0.hom[g]
        assert g in gl, str(g)
        assert h in gl1
        assert g in parabolic
    for g in gl1:
        assert l0.ihom[g] in l0

    Uni = space.get_permrep(uni)
    Uni.do_check()
    assert len(Uni) == p

    L0 = space.get_permrep(l0)
    assert P.is_subgroup(L0)
    assert L0.is_abelian()
    assert len(L0) == p-1
    reps0 = [Rep.fourier(L0, i) for i in range(p-1)]
    #for r in reps0:
    #    print(r, r.chi)

    L1 = space.get_permrep(l1)
    assert P.is_subgroup(L1)
    reps1 = [Rep.fourier(L1, i) for i in range(p-1)]

    print( P )
    print( Uni )
    print( L0 )
    print( L1 )

    L = set(l0*l1 for l0 in L0 for l1 in L1)
    L = Group(L)
    assert P.is_subgroup(L)
    assert P.is_subgroup(Uni)
    assert P.is_normal(Uni)
    assert len(L) * len(Uni) == len(P)

    basis = []
    N = p-1

    for i in range(N):
        reps0[i].name = ascii_names[i]
        reps1[i].name = ascii_names[i]

    for i in range(N):
        r01 = reps0[i]*reps1[i]
        #r01.check()
        #print(r01)
        L = r01.G
        #print( len(L) * len(Uni) , len(P) )
        rep = {}
        for l in L:
          for u in Uni:
            lu = u*l
            assert lu not in rep
            rep[lu] = r01(l)
        rep = Rep(P, rep, r01.dim)
        #print(rep)
        #rep.check()
        rep = rep.induce(G)
        #print(rep)
        a, b = rep.hom(rep)

        r = (a+b).cokernel()
        c = r.tgt
        assert c.is_irrep()
        c.name = r"\yng(1,1)(%s)"%reps0[i].name
        basis.append(c)
        f, = rep.hom(c)
        #print("f:")
        #print(f)
        d = f.cokernel().tgt
        #print(d)
        assert d.is_irrep()
        d.name = r"\yng(2)(%s)"%reps0[i].name
        basis.append(d)

    for i in range(N):
      for j in range(i+1, N):
        r01 = reps0[i]*reps1[j]
        #r01.check()
        print(r01)
        L = r01.G
        #print( len(L) * len(Uni) , len(P) )
        rep = {}
        for l in L:
          for u in Uni:
            lu = u*l
            assert lu not in rep
            rep[lu] = r01(l)
        rep = Rep(P, rep, r01.dim)
        print(rep)
        #rep.check()
        rep = rep.induce(G)
        print(rep)
        assert rep.is_irrep()
        rep.name = r"%s\otimes %s"%(reps0[i].name, reps1[j].name)
        basis.append(rep)
      print()

    for r in basis:
        print(r)
        #print(r.chi)
        r.some_check()
    #return

    #basis = Table(basis)
    basis = Basis(basis)
    basis.show()

    #return

    cuspidals = []
    jdxs = list(range(p**2-1))
    #for j in [0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19]:
    #for j in [1, 2, 3, 4, 7, 8, 9, 13, 14, 19]:
    if p==5:
        jdxs = [1, 2, 3, 4, 7, 8, 9, 13, 14, 19]
    if p==7:
        jdxs = [1,2,3,4,5,6,9,10,11,12,13,17,18,19,20,25,26,27,33,34,41]

    for jdx in jdxs:
        r = Rep.fourier(Torus, jdx)
        r = r.induce(G)
        print(r)
        #r.jdx = jdx
        #cuspidals.append(r)

        basis.append(r)
        basis.show()
        basis.reduce()
        basis.show()

    basis.show()
    basis.reduce(verbose=True)
    basis.show()

    #Basis(cuspidals).show()

#    #return
#    basis = basis.reps + cuspidals
#    basis = Basis(basis)
#    basis.show()
#
#    basis.reduce()
#    basis.show()

    i = p-1
    for chi in basis:
        if chi.name is None:
            chi.name = ascii_names[i]
            i += 1
        print(chi)

    #for cusp in cuspidals:
    for jdx in jdxs:
        r = Rep.fourier(Torus, jdx)
        r = r.induce(G)
        print(r)
        name = []
        for s in basis:
            if r.chi.dot(s.chi):
                name.append(s.name)
        print(jdx, name)

    return basis



def test_monoidal():

    n, p = 3, 2
    space = Space(n, p)
    gl = Algebraic.GL(n,p)

    levis, uni, parabolic = gl.levi_decompose([1,1,1])
    P = space.get_permrep(parabolic)
    G = space.get_permrep(gl)
    rep = Rep.permutation(G, P)
    print(rep)

    fs = rep.hom(rep)
    assert len(fs) == 6 # S_3 acts on A@A@A

    N = rep.dim**2
    rows = []
    for f in fs:
        M = f.M
        M = M.reshape(1,N)
        rows.append(M.M[0])
    A = Matrix(QQ, rows).t
    #print(A.shape)

    for f in fs:
      for g in fs:
        gf = g*f
        if gf in fs:
            print(fs.index(gf), end=" ")
        else:
            print("?", end=" ")
      print()


    def getname(M):
        M = M.reshape(M.shape[0]*M.shape[1], 1)
        a = A.solve(M).t
        return a

    for f in fs:
      for g in fs:
        gf = g*f
        a = getname(gf.M)
        print(a, end=" ")
      print()
    #levis, uni, parabolic = GL.levi_decompose([1,2])
    return

    print()

    I = Matrix.identity(Rep.ring, rep.dim)
    found = []
    #for cs in cross([tuple(range(-4,5))]*len(fs)):
    for cs in cross([tuple(range(-2,2))]*len(fs)):
        if cs[1:] == (0,)*(len(fs)-1):
            continue
        f = reduce(add, [cs[i]*fs[i].M for i in range(len(fs))])
        ff = f*f
        if ff == I:
            print("ff=I:", cs)
        elif ff == f:
            print("ff=f:", cs)
        elif ff == f+2*I:
            print("ff=f+2*I:", cs)
            found.append(f)

    N = len(found)
    print("found:", N)
    pairs = []
    for i in range(N):
      for j in range(i+1, N):
        a = found[i]
        b = found[j]
        if a*b*a == b*a*b:
            print("aba == bab", i, j)
            pairs.append((a,b))

    for (a,b) in pairs:
        alg = [I, a, b, a*b, b*a, a*b*a]
        for g in alg:
            print(getname(g))
        print()
        



class Levi:
    def __init__(self, space, GL, ms):
        levis, uni, parabolic = GL.levi_decompose(ms)
        self.space = space
        self.levis = levis
        self.uni = uni
        self.parabolic = parabolic

    #def induce(self, reps):
            




def test():
    test_rep()
    test_gram_schmidt()
    test_gl()
    test_cuspidal(2,3)


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
        print("%s()"%name)
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        print("%s()"%name)
        fn = eval(name)
        fn()

    else:
        test()


    t = time() - start_time
    print("OK! finished in %.3f seconds\n"%t)




