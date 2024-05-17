#!/usr/bin/env python
"""
Use group theory to build surface codes, color codes, etc.
"""

import string, os
from random import randint, choice, random
from time import sleep, time
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import matmul
import json

import numpy
from numpy import alltrue, zeros, dot

from bruhat.util import cross
from bruhat import solve 
from bruhat.solve import (
    array2, zeros2, shortstr, dot2, solve2, linear_independent, row_reduce, find_kernel,
    span, intersect, rank, enum2, shortstrx, identity2, eq2, pseudo_inverse)
from bruhat.action import Perm, Group, Coset, mulclose, close_hom, is_hom
from bruhat.todd_coxeter import Schreier
from bruhat.argv import argv
from bruhat.smap import SMap
from bruhat import lins_db

infty = "\u221E"

def parse(s):
    for c in "XZY":
        s = s.replace(c, '1')
    s = s.replace("I", "0")
    for c in " [],":
        s = s.replace(c, '')
    return solve.parse(s)

def css_to_symplectic(Lx, Lz):
    assert Lz.shape == Lx.shape
    m, n = Lx.shape
    L = zeros2(2*m, n, 2)
    L[0::2, :, 0] = Lx
    L[1::2, :, 1] = Lz
    return L


def css_to_isotropic(Hx, Hz):
    mx, n = Hx.shape
    mz, n1 = Hz.shape
    assert n==n1
    H = zeros2(mx+mz, n, 2)
    H[:mx, :, 0] = Hx
    H[mx:, :, 1] = Hz
    return H

def flatten(H):
    if H is not None and len(H.shape)==3:
        H = H.view()
        m, n, _ = H.shape
        H.shape = m, 2*n
    return H

def complement(H):
    H = flatten(H)
    H = row_reduce(H)
    m, nn = H.shape
    #print(shortstr(H))
    pivots = []
    row = col = 0
    while row < m:
        while col < nn and H[row, col] == 0:
            #print(row, col, H[row, col])
            pivots.append(col)
            col += 1
        row += 1
        col += 1
    while col < nn:
        pivots.append(col)
        col += 1
    W = zeros2(len(pivots), nn)
    for i, ii in enumerate(pivots):
        W[i, ii] = 1
    #print()
    return W

def get_weight_slow(v):
    count = 0
    for i in range(len(v)//2):
      if v[2*i] or v[2*i+1]:
        count += 1
    c1 = get_weight_fast(v)
    if count != c1:
        print("v =")
        print(v)
        print(count, c1)
        assert 0
    return count

def get_weight_fast(v): # not much faster for n=18
    n = len(v)//2
    v.shape = n,2
    w = v[:,0] + v[:,1]
    v.shape = 2*n,
    return numpy.count_nonzero(w)

get_weight = get_weight_fast


def monte_carlo(H, v, p=0.5, trials=10000):
    H = H.view()
    if len(H.shape) == 3:
        m, n, o = H.shape
        assert o==2
        nn = 2*n
        H.shape = m, nn
    else:
        m, nn = H.shape
    assert v.shape == (nn,), v.shape
    d0 = get_weight_fast(v)
    #print("[",d0, end=",", flush=True)
    p0 = p**d0
    #randint = numpy.random.randint
    for trial in range(trials):
        #u = randint(2, size=m)
        #h = dot2(u, H)
        i = randint(0, m-1)
        h = H[i]
        w = (v+h)%2
        d1 = get_weight_fast(w)
        p1 = p**d1
        a = random()
        if (p0/p1) < a:
            v = w
            d0 = d1
            p0 = p**d0
            #print(d0, end=",", flush=True)
    #print("]")
    return d0

def strop(H):
    assert len(H.shape) == 2, H.shape
    smap = SMap()
    m, nn = H.shape
    for i in range(m):
      for j in range(nn//2):
        x, z = H[i, 2*j:2*j+2]
        c = '.'
        if x and z:
            c = 'Y'
        elif x:
            c = 'X'
        elif z:
            c = 'Z'
        smap[i,j] = c
    return str(smap)

@cache
def symplectic_form(n):
    F = zeros2(2*n, 2*n)
    for i in range(n):
        F[2*i:2*i+2, 2*i:2*i+2] = [[0,1],[1,0]]
    return F


class QCode(object):
    def __init__(self, H, T=None, L=None, J=None, d=None, check=True):
        assert H.max() <= 1
        H = flatten(H) # stabilizers
        m, nn = H.shape
        assert nn%2 == 0
        n = nn//2
        T = flatten(T)
        L = flatten(L)
        J = flatten(J)
        self.H = H # stabilizers
        self.T = T # destabilizers
        self.L = L # logicals
        self.J = J # gauge generators
        self.m = m
        self.n = n
        self.k = n - m # logicals
        self.d = d
        self.shape = m, n
        if L is not None:
            assert L.shape == (2*self.k, nn)
        if T is not None:
            assert T.shape == (self.m, nn)
        if check:
            self.check()

    @classmethod
    def build_css(cls, Hx, Hz, Lx=None, Lz=None, Jx=None, Jz=None, check=True):
        H = css_to_isotropic(Hx, Hz)
        L = css_to_symplectic(Lx, Lz) if Lx is not None else None
        J = css_to_symplectic(Jx, Jz) if Jx is not None else None
        T = None
        return QCode(H, T, L, J, check=check)

    @classmethod
    def build_gauge(cls, J):
        m, n, _ = J.shape
        J1 = J.copy()
        J1[:, :, :] = J1[:, :, [1,0]]
        J1.shape = (m, 2*n)
        K = find_kernel(J1)
        #print("K:", K.shape)
        #print(shortstr(K))
        J = J.view()
        J.shape = m, 2*n
        H = intersect(J, K)
        #print("H:", H.shape)
        #H.shape = len(H), n, 2
        #J.shape = len(J), n, 2
        code = QCode(H, None, None, J)
        return code

    @classmethod
    def build_gauge_css(cls, Jx, Jz):
        J = css_to_isotropic(Jx, Jz)
        code = cls.build_gauge(J)
        return code

    @classmethod
    def fromstr(cls, Hs, Ls=None, check=True):
        stabs = Hs.split()
        H = []
        I, X, Z, Y = [0,0], [1,0], [0,1], [1,1]
        lookup = {'X':X, 'Z':Z, 'Y':Y, 'I':I, '.':I}
        for stab in stabs:
            row = [lookup[c] for c in stab]
            H.append(row)
        H = array2(H)
        return QCode(H, check=check)

    def __str__(self):
        s = "H =\n%s"%strop(self.H)
        T = self.T
        if T is not None and len(T):
            s += "\nT =\n%s"%strop(T)
        L = self.L
        if L is not None and len(L):
            s += "\nL =\n%s"%strop(L)
        return s
    __repr__ = __str__

    def shortstr(self):
        H = self.H.view()
        H.shape = (self.m, 2*self.n)
        return shortstr(H)

    @property
    def deepH(self):
        H = self.H.view()
        m, n = self.shape
        H.shape = m, n, 2
        return H

    @property
    def deepL(self):
        L = self.L.view()
        L.shape = 2*self.k, self.n, 2
        return L

    @property
    def deepT(self):
        T = self.T.view()
        T.shape = self.m, self.n, 2
        return T

    def check(self):
        H = self.H
        m, n = self.shape
        H1 = dot2(H, symplectic_form(n))
        R = dot2(H, H1.transpose())
        if R.sum() != 0:
            assert 0, "not isotropic"
        L = self.get_logops()
        R = dot2(L, H1.transpose())
        if R.sum() != 0:
            print("R:")
            print(shortstr(R))
            assert 0
        R = dot2(L, symplectic_form(n), L.transpose())
        if not eq2(R, symplectic_form(self.k)):
            assert 0
        T = self.get_destabilizers()
        HT = array2(list(zip(H, T)))
        HT.shape = 2*m, 2*n
        A = numpy.concatenate((HT, L))
        F = dot2(A, symplectic_form(n), A.transpose())
        assert eq2(F, symplectic_form(n))

    def dual(self):
        D = array2([[0,1],[1,0]])
        H = dot(self.H, D) % 2
        return QCode(H)

    def overlap(self, other):
        F = symplectic_form(self.n)
        A = dot2(dot2(self.L, F), other.L.transpose())
        #A = dot2(self.L, other.L.transpose())
        A = dot2(symplectic_form(self.k), A)
        return A

    def permute(self, f):
        n = self.n
        H = self.deepH[:, [f[i] for i in range(n)], :].copy()
        L = self.deepL
        if L is not None:
            L = L[:, [f[i] for i in range(n)], :].copy()
        T = self.deepT
        if T is not None:
            T = T[:, [f[i] for i in range(n)], :].copy()
        return QCode(H, T, L)

    def apply(self, idx, gate):
        if idx is None:
            code = self
            for idx in range(self.n):
                code = code.apply(idx, gate) # <--- recurse
            return code
        H = self.deepH.copy()
        H[:, idx] = dot(H[:, idx], gate) % 2
        T = self.T
        if T is not None:
            T = self.deepT.copy()
            T[:, idx] = dot(T[:, idx], gate) % 2
        L = self.L
        if L is not None:
            L = self.deepL.copy()
            L[:, idx] = dot(L[:, idx], gate) % 2
        return QCode(H, T, L)

    def apply_H(self, idx=None):
        # swap X<-->Z on bit idx
        H = array2([[0,1],[1,0]])
        return self.apply(idx, H)

    def apply_S(self, idx=None):
        # swap X<-->Y
        S = array2([[1,1],[0,1]])
        return self.apply(idx, S)

    def apply_SH(self, idx=None):
        # X-->Z-->Y-->X 
        SH = array2([[0,1],[1,1]])
        return self.apply(idx, SH)

    def apply_CZ(self, idx, jdx):
        assert idx != jdx
        n = self.n
        A = identity2(2*n)
        A[2*idx, 2*jdx+1] = 1
        A[2*jdx, 2*idx+1] = 1
        H = dot2(self.H, A)
        T = dot2(self.T, A)
        L = dot2(self.L, A)
        return QCode(H, T, L)

    def row_reduce(self):
        H = self.H.copy()
        m, n = self.shape
        H.shape = m, 2*n
        H = row_reduce(H)
        m, nn = H.shape
        return QCode(H)

#    def get_logops(self):
#        if self.L is not None:
#            return self.L
#        H = self.H
#        H1 = dot2(H, symplectic_form(self.n))
#        H = row_reduce(H)
#        m = len(H)
#        #print("H:", H.shape)
#        #print(shortstr(H))
#        W = complement(H)
#        #print("W:", W.shape)
#        #print(shortstr(W))
#        K = find_kernel(H1)
#        #print("find_kernel(H):", K.shape)
#        L = intersect(W, K)
#        #print("L:", L.shape)
#        #print("L")
#        #print(shortstr(L))
#        kk = len(L)
#        assert kk%2 == 0
#        k = kk//2
#        assert dot2(H1, L.transpose()).sum() == 0 # commutes w stabilizers
#        HL = numpy.concatenate((H, L))
#        assert rank(HL) == m + 2*k # linearly independant
#        J = self.J
#        if J is not None:
#            J1 = J.view()
#            J1.shape = len(J1), J1.shape[1]*2
#            Jc = complement(J)
#            #print("Jc:", Jc.shape)
#            #print("rank(J):", rank(J1))
#            L = intersect(L, Jc)
#            #print("L:", L.shape)
#        self.L = L.copy()
#        #self.L.shape = (2*k, self.n, 2)
#        return L

    def get_logops(self):
        if self.L is not None:
            return self.L
        m, n = self.shape
        k = self.k
        H = self.H
        L = []
        for i in range(k):
            M = dot2(H, symplectic_form(n))
            K = find_kernel(M) # XXX do this incrementally for more speed XXX
            J = dot2(K, symplectic_form(n), K.transpose())
            for (row, col) in zip(*numpy.where(J)):
                break
            lx = K[row:row+1]
            lz = K[col:col+1]
            L.append(lx)
            L.append(lz)
            H = numpy.concatenate((H, lx, lz))
    
        L = array2(L)
        L.shape = (2*k, 2*n)
    
        M = dot2(L, symplectic_form(n), L.transpose())
        assert eq2(M, symplectic_form(k))
        self.L = L
        return L

    def get_destabilizers(self):
        if self.T is not None:
            return self.T
        m, n = self.shape
        H = self.H
        L = self.get_logops()
        HL = numpy.concatenate((H, L))
        T = []
        for i in range(m):
            A = dot2(HL, symplectic_form(n))
            #B = numpy.concatenate((identity2(m), zeros2(2*self.k+i, m)))
            B = zeros2(len(A), 1)
            B[i] = 1
            #print("A =")
            #print(A)
            tt = solve2(A, B) # XXX do this incrementally for more speed XXX
            assert tt is not None
            t = tt.transpose()
            t.shape = (1, 2*n)
            T.append(t)
            HL = numpy.concatenate((HL, t))
        T = array2(T)
        T.shape = (m, 2*n)
        self.T = T

        HT = array2(list(zip(H, T)))
        HT.shape = 2*m, 2*n
        A = numpy.concatenate((HT, L))
        M = dot2(A, symplectic_form(n), A.transpose())
        assert eq2(M, symplectic_form(n))

        return T

    def build(self):
        self.get_logops()
        self.get_destabilizers()

    def get_all_logops(self):
        L = self.get_logops()
        H = self.H
        HL = numpy.concatenate((H, L))
        return HL

    def get_params(self, max_mk=22):
        L = self.get_logops()
        kk = len(L)
        H = self.H
        m = len(H)
        HL = numpy.concatenate((H, L))
        mk = len(HL)
        if mk > max_mk:
            return self.n, kk//2, None
        d = self.n
        for w in numpy.ndindex((2,)*mk):
            u, v = w[:m], w[m:]
            if sum(v) == 0:
                continue
            v = dot2(w, HL)
            count = get_weight(v) # destructive
            if count:
                d = min(count, d)
        return self.n, kk//2, d

    def bound_distance(self):
        L = self.get_logops()
        L = L.copy()
        H = self.H
        kk = len(L)
        #H = self.H
        d = self.n
        for u in L:
            d = min(d, u.sum())
            w = monte_carlo(H, u)
            d = min(w, d)
        return d

    def unwrap(self, check=True):
        H0 = self.deepH
        m, n, _ = H0.shape
        Sx = H0[:, :, 0]
        Sz = H0[:, :, 1]
        Sxz = numpy.concatenate((Sx, Sz), axis=1)
        Szx = numpy.concatenate((Sz, Sx), axis=1)
        H = zeros2(2*m, 2*n, 2)
        H[:m, :, 0] = Sxz
        H[m:, :, 1] = Szx
        code = QCode(H, check=check)
        return code

    
def find_code(H, css=False, yop=False, degen=False):
    "build a QCode from face--vert adjacency matrix"
    import z3
    from z3 import Bool, And, Or, Not, Implies, If, Solver

    m, n = H.shape
    vs = {}
    xvars = {}
    zvars = {}
    yvars = {}

    clauses = []
    solver = Solver()

    yop = False if css else yop

    add = solver.add

    for i in range(m):
      for j in range(n):
        if H[i, j] == 0:
            continue
        X = Bool("X_%d_%d"%(i,j))
        Z = Bool("Z_%d_%d"%(i,j))
        vs[i, j] = (X, Z)
        xvars[i,j] = X
        zvars[i,j] = Z
        yvars[i,j] = And(X, Z)
        # X or Z
        add(Or(X, Z))
        if not yop:
            add(Not(And(X, Z))) # not both

    # each qubit has X & Z check
    for j in range(n): # cols
        idxs = [i for i in range(m) if H[i, j]] # rows
        assert len(idxs) > 1, idxs
        add(Or(*[xvars[i,j] for i in idxs]))
        add(Or(*[zvars[i,j] for i in idxs]))
        add(Not(And(*[yvars[i,j] for i in idxs])))
        # yes?

    if degen:
      # product of all stabilizers is identity
      for j in range(n): # cols
        checks = [i for i in range(m) if H[i, j]] # rows
        assert len(checks)>1
        for k in [0, 1]: # X, Z
          clauses = []
          for bits in cross([(0,1)]*len(checks)):
            clause = []
            for idx, bit in enumerate(bits):
                term = vs[checks[idx], j][k] 
                term = term if bit else Not(term)
                clause.append(term)
            clause = And(*clause)
            parity = sum(bits)%2
            if not parity:
                clauses.append(clause)
          clause = Or(*clauses)
          solver.add(clause)


    if css:
      for i in range(m):
        xbits, zbits = [], []
        for j in range(n):
            if H[i,j] == 0:
                continue
            xbits.append(xvars[i,j])
            zbits.append(zvars[i,j])
        clause = Or(And(*xbits), And(*zbits))
        solver.add(clause)

    # commuting stabilizers
    for i0 in range(m):
      for i1 in range(i0+1, m):
        bits = [j for j in range(n) if H[i0,j]>0 and H[i1,j]>0]
        if not bits:
            continue
        #print(bits)
        for bits0 in cross([((1,0),(0,1),(1,1))]*len(bits)):
         for bits1 in cross([((1,0),(0,1),(1,1))]*len(bits)):
            syndrome = sum(l[0]*r[1]+l[1]*r[0] for (l,r) in zip(bits0, bits1))
            #print(bits0, bits1, syndrome)
            if syndrome%2 == 0:
                continue # allowed
            clause = []
            for j, jj in enumerate(bits):
                X0, Z0 = vs[i0, jj]
                X1, Z1 = vs[i1, jj]
                clause.append(X0 if bits0[j][0] else Not(X0))
                clause.append(Z0 if bits0[j][1] else Not(Z0))
                clause.append(X1 if bits1[j][0] else Not(X1))
                clause.append(Z1 if bits1[j][1] else Not(Z1))
            clause = Not(And(*clause)) # not allowed
            #print(str(clause).replace(" ", "").replace("\n", ""))
            solver.add(clause)

    result = solver.check()
    if result != z3.sat:
        #print(result)
        return None

    H = zeros2(m, n, 2)
    model = solver.model()
    for (i, j), (X, Z) in vs.items():
        X = model.evaluate(X)
        Z = model.evaluate(Z)
        assert X or Z
        if X:
            H[i, j, 0] = 1
        if Z:
            H[i, j, 1] = 1
    code = QCode(H)
    return code


def make_dot(graph, cmd="neato"):
    n = len(graph)
    gens = graph.get_gens()

    f = open("schreier.dot", "w")
    labels = string.ascii_lowercase
    if len(labels) < n:
        labels = []
        i = 0
        while len(labels) < n:
            for c in string.ascii_lowercase:
                labels.append("%s%d"%(c,i))
            i += 1
    print("graph\n{\n", file=f)
    #cls = "green blue red".split()
    cls = "red blue green".split()
    for k in range(3):
        gen = gens[k]
        for i in gen.items:
            if i <= gen(i):
                print("    %s -- %s [color=%s];" % (
                    labels[i], labels[gen(i)], cls[k]), file=f)
    print("}\n", file=f)
    f.close()

    os.system("%s -Tpdf schreier.dot > schreier.pdf" % cmd)
    #data = os.popen("dot -Tpdf schreier.dot").read()
    #print("data:", data)
    #open("schreier.pdf", "w").write(data)



class Geometry(object):
    "A geometry specified by a Coxeter reflection group"
    def __init__(self, orders, lins_idx=0, build=True):
        ngens = len(orders)+1
        a, b, c, d, e = [(i,) for i in range(5)]
        orders = tuple(orders)
        assert orders in lins_db.db, (list(lins_db.db.keys()))
        #print("oeqc.Geometry.__init__:", len(lins_db.db[orders]))
        rels = lins_db.db[orders][lins_idx]
        rels = lins_db.parse(rels, **locals())
        self.orders = orders
        self.ngens = ngens
        self.rels = rels
        self.dim = len(orders)
        if build:
            self.G = self.build_group()

    def build_graph(self, figure=None, hgens=None):
        gens = [(i,) for i in range(5)]
        ngens = self.ngens
        rels = [gens[i]*2 for i in range(ngens)]
        orders = self.orders
        for i in range(ngens-1):
            order = orders[i]
            if order is not None:
                rels.append( (gens[i]+gens[i+1])*order )
            for j in range(i+2, ngens):
                rels.append( (gens[i]+gens[j])*2 )
        rels = rels + self.rels
        #print(rels)
        graph = Schreier(ngens, rels)
        if figure is not None:
            assert len(figure) == len(gens)
            gens = [gens[i] for i, fig in enumerate(figure) if fig] or [G.identity]
            graph.build(gens)
        elif hgens is not None:
            graph.build(hgens)
        else:
            graph.build()
        return graph

    def build_group(self):
        graph = self.build_graph()
        G = graph.get_group()
        return G

    def get_cosets(self, figure):
        G = self.G
        gens = G.gens
        assert len(figure) == len(gens)
        gens = [gens[i] for i, fig in enumerate(figure) if fig] or [G.identity]
        #print("gens:", gens)
        H = Group.generate(gens)
        #pairs = G.left_cosets(H)
        cosets = G.left_cosets(H)
        return cosets


def get_adj(left, right):
    A = zeros2((len(left), len(right)))
    for i, l in enumerate(left):
      for j, r in enumerate(right):
        lr = l.intersection(r)
        A[i, j] = len(lr)>0
    return A


def build_code(geometry):
    dim = geometry.dim
    print("dim =", dim)

    if dim == 2:
        faces = geometry.get_cosets([0,1,1])
        edges = geometry.get_cosets([1,0,1])
        verts = geometry.get_cosets([1,1,0])
        print("faces=%d, edges=%d, verts=%d"%(len(faces), len(edges), len(verts)))
    else:
        bodis = geometry.get_cosets([0,1,1,1])
        faces = geometry.get_cosets([1,0,1,1])
        edges = geometry.get_cosets([1,1,0,1])
        verts = geometry.get_cosets([1,1,1,0])
        partial_flags = [
            bodis, faces, edges, verts,
            geometry.get_cosets([0,0,1,1]),
            geometry.get_cosets([0,1,0,1]),
            geometry.get_cosets([0,1,1,0]),
            geometry.get_cosets([1,0,0,1]),
            geometry.get_cosets([1,0,1,0]),
            geometry.get_cosets([1,1,0,0]),
            geometry.get_cosets([1,0,0,0]),
            geometry.get_cosets([0,1,0,0]),
            geometry.get_cosets([0,0,1,0]),
            geometry.get_cosets([0,0,0,1]),
        ]

    Hx = Hz = None
    Jx = Jz = None
    code = None

    if argv.homology == 1:
        Hz = get_adj(faces, edges)
        Hx = get_adj(verts, edges)

    elif argv.homology == 2:
        Hz = get_adj(bodis, faces)
        Hx = get_adj(edges, faces)

    elif argv.flag and dim==2:
        flags = geometry.get_cosets([0]*(dim+1))
        H0 = get_adj(faces, flags)
        H1 = get_adj(edges, flags)
        H2 = get_adj(verts, flags)
        #print(shortstr(H0), H0.shape, H0.sum(1)) # on 6,4 lattice weight 8
        #print(shortstr(H1), H1.shape, H1.sum(1)) # on 6,4 lattice weight 4
        #print(shortstr(H2), H2.shape, H2.sum(1)) # on 6,4 lattice weight 12
        Hx = numpy.concatenate((H0, H1, H2))
        Hz = Hx.copy()

    elif argv.selfdual and dim==2:
        A = get_adj(faces, verts)
        #print(A.shape)
        H = linear_independent(A)
        if argv.dump:
            print(shortstr(H))
        if argv.find_lower_distance:
            from bruhat.hecke import find_lower_distance, get_logops
            dd = argv.get("d", 4)
            L = get_logops(H)
            result = find_lower_distance(H, L, dd)
            if result:
                print("found weight %d logicals"%dd)
            else:
                print("no weight %d logicals"%dd)
            #return
        if argv.distance:
            from bruhat.hecke import get_lower_distance, find_lower_distance, get_logops
            from bruhat.hecke import find_distance_dynamic_big
            L = get_logops(H)
            d = find_distance_dynamic_big(A, L)
            print("d <=", d)
            for dd in [2, 3, 4, 5, 6, 7]:
                #if dd<7: continue
                result = find_lower_distance(H, L, dd)
                print("d = %d? %s" % (dd, result))
                if result:
                    break
        Hx = H.copy()
        Hz = H.copy()

    elif argv.flag and dim==3:
        flags = geometry.get_cosets([0]*(dim+1))
        H0 = get_adj(bodis, flags)
        H1 = get_adj(faces, flags)
        H2 = get_adj(edges, flags)
        H3 = get_adj(verts, flags)
        #print(shortstr(H0))
        #print(shortstr(H1))
        #print(shortstr(H2))
        Hx = numpy.concatenate((H0, H1, H2, H3))
        Hz = Hx.copy()

    elif argv.cells:
        A = get_adj(faces, edges)
        B = get_adj(edges, verts)
        print("A =")
        print(shortstr(A), A.shape)
        print("B =")
        print(shortstr(B), B.shape)

    elif argv.partial_flags and dim==3:
        flags = geometry.get_cosets([0]*(dim+1))
        Hs = [get_adj(p, flags) for p in partial_flags]
        N = len(partial_flags)
        #for i in range(N):
        # for j in range(N):
        #    print("%4s"%dot2(Hs[i], Hs[j].transpose()).sum(), end=" ")
        # print()
        Hx = numpy.concatenate(tuple(Hs[:4]))
        Hz = numpy.concatenate(tuple(Hs[:10]))

    elif argv.gauge and dim==3:
        flags = geometry.get_cosets([0]*(dim+1))
        Hs = [get_adj(p, flags) for p in partial_flags]
        N = len(partial_flags)
        for i in range(N):
         for j in range(N):
            print("%4s"%dot2(Hs[i], Hs[j].transpose()).sum(), end=" ")
         print()
        print()
        J = []
        for i in range(N):
            weights = list(Hs[i].sum(1))
            if weights.count(4) == len(weights) or weights.count(6) == len(weights):
                J.append(Hs[i])
        J = numpy.concatenate(J)
        #print(shortstr(J))
        Jx = J
        Jz = Jx.copy()
    else:
        assert 0

    if Hx is not None:
        Ax = Hx.copy()
        Az = Hz.copy()

        Hx = linear_independent(Hx)
        Hz = linear_independent(Hz)

        A = dot2(Hx, Hz.transpose())
        #print("chain condition:", A.sum() == 0)
    
        if A.sum() != 0:
            print("not commutative\n")
            return
    
        Hxs = Hx.sum(1)
        Hzs = Hz.sum(1)
        print("Hx weights:", Hxs.min(), "to", Hxs.max())
        print("Hz weights:", Hzs.min(), "to", Hzs.max())
    
        code = QCode.build_css(Hx, Hz, check=argv.check)
        print("rank(Hx) =", rank(Hx))
        print("rank(Hz) =", rank(Hz))

    elif Jx is not None:
        Jxs = Jx.sum(1)
        Jzs = Jz.sum(1)
        print("Jx weights:", Jxs.min(), "to", Jxs.max())
        print("Jz weights:", Jzs.min(), "to", Jzs.max())
    
        code = QCode.build_gauge_css(Jx, Jz)
    else:
        return

    n, k, d = code.get_params(max_mk=15)
    if d is None:
        #L = code.get_logops()
        #print(list(L.sum(1)))
        #print(shortstr(L))
        #print("L", L.shape)
        d = code.bound_distance()
        print("[[%d, %d, d<=%d]]" % (n, k, d))
    else:
        print("[[%d, %d, %s]]" % (n, k, d))

    if argv.show:
        print(code)
    print()

    if not argv.autos:
        return code

    #code.build()

    test_zx(Ax, Az)

    return code


def build_geometry():

    key = argv.get("key", (4,3,4))
    print("key:", key)

    idx = argv.get("idx")
    idxs = [idx] if idx else list(range(1000))

    for idx in idxs:
        try:
            geometry = Geometry(key, idx)
        except IndexError:
            break

        G = geometry.G
        print("|G| = %d, idx = %d" % (len(G), idx))
        if len(G)<10:
            continue

        code = build_code(geometry)

    #print("build_geometry: idx =", idx)


def test_zx(Ax, Az):
    ax, n = Ax.shape
    az, _ = Az.shape
    print("Ax =")
    print(shortstr(Ax))
    print("Az =")
    print(shortstr(Az))

    from bruhat.isomorph import Tanner, search
    print(Ax.shape, Az.shape)
    U = Tanner.build2(Ax, Az)
    V = Tanner.build2(Ax, Az)
    perms = []
    for f in search(U, V):
        perms.append(f)
    print("autos:", len(perms))
    assert len(f) == n+ax+az
    items = list(range(len(f)))
    # each perm acts on ax+az+n in that order
    perms = [Perm(perm, items) for perm in perms]
    G = Group(perms, items, check=True)
    print(G)
    for g in G:
        for i in range(ax):
            assert 0<=g[i]<ax
        for i in range(ax, ax+az):
            assert ax<=g[i]<ax+az

    if 1:
        # find a zx duality
        U = Tanner.build2(Ax, Az)
        V = Tanner.build2(Az, Ax)
        for duality in search(U, V):
            break
        else:
            assert 0
        print([duality[i] for i in range(ax+az+n)])
    
        dual_x = [duality[i] for i in range(ax)]
        dual_z = [duality[i]-ax for i in range(ax,ax+az)]
        dual_n = [duality[i]-ax-az for i in range(ax+az,ax+az+n)]
    
        print(dual_x, dual_z, dual_n)
        Ax = Ax[dual_z, :][:, dual_n]
        Az = Az[dual_x, :][:, dual_n]
        print("Ax =")
        print(shortstr(Ax))
        print("Az =")
        print(shortstr(Az))

    import bruhat.solve
    from qupy.ldpc import solve
    solve.int_scalar = bruhat.solve.int_scalar
    from bruhat.algebraic import Matrix, Algebraic

    from qupy.ldpc.css import CSSCode
    Hx = linear_independent(Ax)
    Hz = linear_independent(Az)
    code = CSSCode(Hx=Hx, Hz=Hz, check=True)
    print(code)
    k = code.k
    Lx, Lz = code.Lx, code.Lz
    #print("Lx =")
    #print(shortstr(Lx))
    #print("Lz =")
    #print(shortstr(Lz))
    assert eq2(dot2(Lx, Lz.transpose()), identity2(k))

    qcode = QCode.build_css(Hx, Hz, Lx, Lz)
    L = qcode.L
    Fn = symplectic_form(n)
    Fk = symplectic_form(k)
    assert eq2(dot2(dot2(L, Fn), L.transpose()), Fk)
    print(qcode)

    #tgt = qcode.apply_H()
    #print(tgt)
    #print(qcode.overlap(tgt))

    gen = []
    perms = [{i-ax-az:f[i]-ax-az for i in range(ax+az,ax+az+n)} for f in perms]
    items = list(range(n))
    perms = [Perm(f, items) for f in perms]
    zx0 = Perm(dual_n, items)
    phase_gates = []
    for f in perms:
        zx = zx0*f
        if (zx*zx).is_identity() and len(zx.fixed())%2==0:
            phase_gates.append(zx)
        tgt = qcode.permute(f)
        #print(tgt)
        S = qcode.overlap(tgt)
        #print(shortstr(S))
        #print()
        lhs = dot2(S, Fk, S.transpose())
        assert eq2(lhs, Fk)
        S = Matrix(S)
        gen.append(S)
    print("phase_gates:", len(phase_gates))

    gen = set(gen)
    print("|gen| =", len(gen))
    for g in gen:
      for h in gen:
        assert g*h in gen

    gen = list(gen)
    while 1:
        gens = [choice(gen), choice(gen), choice(gen)]
        G = mulclose(gens)
        if len(G) == len(gen):
            break

    fixed = [i for i in range(n) if dual_n[i]==i]
    #print("dual_n:", dual_n)
    #print("fixed:", fixed)

    # hadamard gate
    tgt = qcode.permute(dual_n)
    tgt = tgt.apply_H()
    #print(tgt)
    S = qcode.overlap(tgt)
    #print(shortstr(S))
    lhs = dot2(S, Fk, S.transpose())
    assert eq2(lhs, Fk)
    S = Matrix(S)
    gens.append(S)

    if phase_gates:
        # phase-type gate
        tgt = qcode
        zx = phase_gates[0]
        remain = set(range(n))
        for i in zx.fixed():
            tgt = tgt.apply_S(i)
            remain.remove(i)
        for i in range(n):
            if i not in remain:
                continue
            j = zx[i]
            assert zx[j] == i
            remain.remove(i)
            remain.remove(j)
            tgt = tgt.apply_CZ(i, j)
        print(tgt)
        S = qcode.overlap(tgt)
        #print(shortstr(S))
        lhs = dot2(S, Fk, S.transpose())
        assert eq2(lhs, Fk)
        S = Matrix(S)
        gens.append(S)

    else:
        G = Algebraic(list(gens))
        print("|G| =", len(G))

    from bruhat.orthogonal import gap_code
    gap_code(gens)


def main_zx():
    Ax = parse("""
    .1.....111..1..
    ......1...111.1
    1..1.....1.1.1.
    .11.1........11
    1.1..1..1.1....
    ...11111.......
    """)
    Az = parse("""
    .111.1...1.....
    .1......1.11.1.
    1...11.....1..1
    1..1...1..1.1..
    ......11.1...11
    ..1.1.1.1...1..
    """)
    test_zx(Ax, Az)


def main_unwrap():

    stabs = "XZZXI IXZZX XIXZZ ZXIXZ" # five qubit code
    code = QCode.fromstr(stabs)
    print(code)
    print(code.get_params())

    code = code.unwrap()
    print(code)
    print(code.get_params())

    # 13,1,5
    op = "XZIIIZXIIIIII"
    n = len(op)
    assert n==13
    stabs = ' '.join(''.join(op[(j+i)%n] for j in range(n)) for i in range(n-1))
    code = QCode.fromstr(stabs, check=True)
    print(code)
    print(code.get_params())

    code = code.unwrap()
    print(code)
    print(code.get_params())
    return

    f = open("codetables.txt")
    data = {}
    for line in f:
        record = json.loads(line)
        n = record['n']
        k = record['k']
        data[n,k] = record
        if n > 8:
            break
        stabs = record["stabs"]
        code = QCode.fromstr(stabs)
        print(code, code.get_params())
        code = code.unwrap()
        print(code, code.get_params())
        print()


def test_build():

    f = open("codetables.txt")
    data = {}
    for line in f:
        record = json.loads(line)
        n = record['n']
        k = record['k']
        data[n,k] = record
        if n > 10:
            break
        stabs = record["stabs"]
        print("stabs:")
        print(stabs)
        code = QCode.fromstr(stabs)
        code.build()
        print(code, code.get_params())
        print()


def main_symplectic_unwrap_5():
    from qupy.ldpc.css import CSSCode
    Hx = parse("""
    X..X..XX..
    .X..X..XX.
    X.X.....XX
    .X.X.X...X
    """)
    Tz = parse("""
    ...ZZ...ZZ
    Z..Z.Z..Z.
    .Z.Z..Z.Z.
    ..ZZ...ZZ.
    """)
    Hz = parse("""
    .ZZ..Z..Z.
    ..ZZ..Z..Z
    ...ZZZ.Z..
    Z...Z.Z.Z.
    """)
    Tx = parse("""
    ...XX...XX
    X..X.X..X.
    .X.X..X.X.
    ..XX...XX.
    """)
    Lx = parse("""
    XXXXX.....
    .....XXXXX
    """)
    Lz = parse("""
    ZZZZZ.....
    .....ZZZZZ
    """)

    code = CSSCode(Hx=Hx, Hz=Hz, Tx=Tx, Tz=Tz, Lx=Lx, Lz=Lz, check=True)
    #print(code.longstr())

    n = code.n

    H = ("""
    X..X..XX..
    .X..X..XX.
    X.X.....XX
    .X.X.X...X
    XXXXX.....
    .ZZ..Z..Z.
    ..ZZ..Z..Z
    ...ZZZ.Z..
    Z...Z.Z.Z.
    .....ZZZZZ
    """)

    T = ("""
    ...ZZ...ZZ
    Z..Z.Z..Z.
    .Z.Z..Z.Z.
    ..ZZ...ZZ.
    ZZZZZ.....
    ...XX...XX
    X..X.X..X.
    .X.X..X.X.
    ..XX...XX.
    .....XXXXX
    """)

    HC = QCode.fromstr(H, check=True)
    TC = QCode.fromstr(T, check=True)

    F = symplectic_form(n)

    H = HC.H
    T = TC.H
    for i in range(5):
        H[i] += H[i+5]
        T[i+5] += T[i]

    M = list(zip(H, T))
    M = numpy.array(M)
    M.shape = (2*n, 2*n)

    print(shortstr(M))

    #print()
    #for i in range(10):
    #    M[i] += M[i+10]
    #print(shortstr(M))

    # M is symplectic:
    A = (dot2(dot2(M, F), M.transpose()))
    assert eq2(A, F)

    Mi = dot2(dot2(F, M.transpose(), F)) # inverse
    assert eq2(dot2(Mi, M), identity2(2*n))

    print()

    # ---------------------------------------

    n = 5
    
    H = ("""
    XZZX.
    .XZZX
    X.XZZ
    ZX.XZ
    XXXXX
    """)
    T = ("""
    ...YY
    Y..Y.
    .Y.Y.
    ..YY.
    ZZZZZ
    """)

    Hx = QCode.fromstr(H, check=True)
    Tz = QCode.fromstr(T, check=True)
    Hz = Hx.apply_H()
    Tx = Tz.apply_H()

    F = symplectic_form(n)

    Hx = Hx.H
    Tz = Tz.H
    Hz = Hz.H
    Tx = Tx.H

    M = list(zip(Hx, Tz))
    M = numpy.array(M)
    M.shape = (2*n, 2*n)

    #print(shortstr(M))

    # M is symplectic:
    A = (dot2(dot2(M, F), M.transpose()))
    assert eq2(A, F)
    
    # ---------------------------------------

    n = 5
    #H = numpy.block([[H, zeros2(n,2*n)], [zeros2(n,2*n), H]])
    #T = numpy.block([[T, zeros2(n,2*n)], [zeros2(n,2*n), T]])
    #H = numpy.block([[Hx, Hx], [zeros2(n,2*n), Hx]])
    #T = numpy.block([[Tz, zeros2(n,2*n)], [Tz, Tz]])
    H = numpy.block([[Hx, Hz], [zeros2(n,2*n), Hz]])
    T = numpy.block([[Tz, zeros2(n,2*n)], [Tz, Tx]])

    n = 10
    F = symplectic_form(n)

    M1 = list(zip(H, T))
    M1 = numpy.array(M1)
    M1.shape = (2*n, 2*n)

    print(shortstr(M1))

    # M1 is symplectic:
    A = (dot2(dot2(M1, F), M1.transpose()))
    assert eq2(A, F)

    A = dot2(M1, Mi)
    print()
    print(shortstr(A), A.sum())
    print(A.shape)

    pairs = [(i,j) for (i,j) in zip(*numpy.where(A)) if i!=j]
    print(pairs)
    print(len(pairs))
    remain = list(pairs)
    for (i,j) in pairs:
        if i>j:
            continue
        if i%2 != j%2:
            continue
        assert i%2
        a,b = i//2, j//2
        print("CNOT(%s, %s)"%(a,b))
        assert (2*b, 2*a) in pairs
        remain.remove((i,j))
        remain.remove((2*b, 2*a))
    print(remain, len(remain))

    for (i,j) in list(remain):
        if i>j or i%2==0:
            continue
        assert j%2==0
        a, b = i//2, j//2
        print("CZ(%s, %s)"%(a,b))
        #print((a,b), (i,j), (2*a+1, 2*b) in remain, (2*b+1, 2*a) in remain)
        remain.remove((i, j))
        remain.remove((2*b+1, 2*a))
    print(remain, len(remain))

    for (i,j) in list(remain):
        assert i%2 == 0
        assert j%2 == 1
        if i>j:
            continue
        a, b = i//2, j//2
        print("CZH(%s, %s)"%(a,b))
        remain.remove((i, j))
        remain.remove((2*b, 2*a+1))
    assert not remain

    #B = A - identity2(2*n)
    #print(shortstr(B))
    #print()

#    W = parse("""
#    1...
#    11..
#    ..11
#    ...1
#    """)
#    print(W)
#    W = W[[0,2,1,3],:][:, [0,2,1,3]]
#    print(W)

    def cnot(i, j):
        assert i!=j
        W = identity2(2*n)
        W[2*j, 2*i] = 1
        W[2*i+1, 2*j+1] = 1
        return W

    def cz(i, j):
        assert i!=j
        W = identity2(2*n)
        W[2*j+1, 2*i] = 1
        W[2*i+1, 2*j] = 1
        return W
    
    def czh(i, j):
        assert i!=j
        W = identity2(2*n)
        W[2*j, 2*i+1] = 1
        W[2*i, 2*j+1] = 1
        return W
    
    n = 2
    #print(czh(0,1))

    #print(shortstr(cnot(0, 6)))
    
    # ---------------------------------------


def main_symplectic_unwrap_13():

    n = 13
    H = """
    XZ...ZX......
    Z...ZX......X
    ...ZX......XZ
    ..ZX......XZ.
    .ZX......XZ..
    ZX......XZ...
    X......XZ...Z
    ......XZ...ZX
    .....XZ...ZX.
    ....XZ...ZX..
    ...XZ...ZX...
    ..XZ...ZX....
    XXXXXXXXXXXXX
    """

    Hx = QCode.fromstr(H, check=True)
    HHx = Hx.H

    rows = []
    for idx in range(n-1):
      for bits in numpy.ndindex((2,)*n):
        if sum(bits)%2:
            continue
        v = zeros2((1,2*n))
        v[0, ::2] = bits
        v[0, 1::2] = bits
        #print(shortstr(v))
        A = dot2(HHx, v.transpose())
        if A[idx] == 1 and A.sum() == 1:
            print(idx, bits, sum(bits))
            rows.append(''.join('.Y'[ii] for ii in bits))
    print('\n'.join(rows))

    T = """
    ..Y..Y..YY.YY
    .YY.YY.Y.YY.Y
    YYYYYYY.....Y
    YY.YY...YY...
    Y..Y.Y.Y.Y.Y.
    ....YYY..YYY.
    ..YYY.....YYY
    .Y.Y.Y..Y.Y.Y
    Y...YY.YY...Y
    ..YYYYYYYY...
    .Y.YY.YY.Y.YY
    Y..Y..Y..YY.Y
    ZZZZZZZZZZZZZ
    """

    Tz = QCode.fromstr(T, check=True)
    TTz = Tz.unwrap()
    HHx = Hx.unwrap()

    Hz = Hx.apply_H()
    Tx = Tz.apply_H()

    F = symplectic_form(n)

    Hx = Hx.H
    Tz = Tz.H
    Hz = Hz.H
    Tx = Tx.H

    M = list(zip(Hx, Tz))
    M = numpy.array(M)
    M.shape = (2*n, 2*n)

    #print(shortstr(M))

    # M is symplectic:
    A = (dot2(dot2(M, F), M.transpose()))
    assert eq2(A, F)
    
    # ---------------------------------------
    # two copies of the 13 qubit code

    #H = numpy.block([[H, zeros2(n,2*n)], [zeros2(n,2*n), H]])
    #T = numpy.block([[T, zeros2(n,2*n)], [zeros2(n,2*n), T]])
    #H = numpy.block([[Hx, Hx], [zeros2(n,2*n), Hx]])
    #T = numpy.block([[Tz, zeros2(n,2*n)], [Tz, Tz]])
    H = numpy.block([[Hx, Hz], [zeros2(n,2*n), Hz]])
    T = numpy.block([[Tz, zeros2(n,2*n)], [Tz, Tx]])

    n = 26
    F = symplectic_form(n)

    M1 = list(zip(H, T))
    M1 = numpy.array(M1)
    M1.shape = (2*n, 2*n)

    #print(shortstr(M1))

    # M1 is symplectic:
    A = (dot2(dot2(M1, F), M1.transpose()))
    assert eq2(A, F)

    # ---------------------------------------

    n = 26
    F = symplectic_form(n)

    print(HHx)
    print(TTz)

    H = HHx.H
    T = TTz.H

    # swap X/Z ops
    T = numpy.concatenate((T[n//2:], T[:n//2]))

    # change basis
    for i in range(13):
        H[i] += H[i+13]
        T[i+13] += T[i]

    M = list(zip(H, T))
    M = numpy.array(M)
    M.shape = (2*n, 2*n)

    # M is symplectic:
    A = (dot2(dot2(M, F), M.transpose()))
    assert eq2(A, F)

    Mi = dot2(dot2(F, M.transpose(), F)) # inverse
    assert eq2(dot2(Mi, M), identity2(2*n))

    # ---------------------------------------

    A = dot2(M1, Mi)
    print()
    print(shortstr(A), A.sum())
    print(A.shape)

    pairs = [(i,j) for (i,j) in zip(*numpy.where(A)) if i!=j]
    print(pairs)
    print("2-qubit gates:", len(pairs)//2)
    
    # 162 2-qubit gates... too many 




def monte_carlo_css(H, v, p=0.5, trials=10000):
    m, n = H.shape
    assert v.shape == (n,)
    d0 = v.sum()
    #print("[",d0, end=",", flush=True)
    p0 = p**d0
    best = v
    #randint = numpy.random.randint
    for trial in range(trials):
        #u = randint(2, size=m)
        #h = dot2(u, H)
        i = randint(0, m-1)
        h = H[i]
        w = (v+h)%2
        d1 = w.sum()
        p1 = p**d1
        a = random()
        if (p0/p1) < a:
            v = w
            d0 = d1
            p0 = p**d0
            best = w
            #print(d0, end=",", flush=True)
    #print("]")
    return best


def make_ramified():
    import bruhat.solve
    from qupy.ldpc import solve
    solve.int_scalar = bruhat.solve.int_scalar
    from qupy.ldpc.css import CSSCode

    print("make_ramified")
    key = (3, 8)
    idx = argv.get("idx", 10)
    geometry = Geometry(key, idx, False)

    hgens = [(1,2)*4]
    #hgens = []
    graph = geometry.build_graph(hgens=hgens)

    graph = graph.compress()
    print(len(graph))

    faces = graph.components([(1,), (2,)])
    print(len(faces))
    print([len(c) for c in faces])

    verts = graph.components([(0,), (1,)])
    print(len(verts))
    print([len(c) for c in verts])

    H = get_adj(faces, verts)
    print(shortstr(H))
    print(dot2(H, H.transpose()))

    code = CSSCode(Hx=H, Hz=H)
    print(code)



def make_colour():
    print("make_colour")
    key = argv.get("key", (3, 8))
    idx = argv.get("idx")
    idxs = argv.get("idxs")
    max_idx = argv.get("max_idx", 35)

    if idxs is not None:
        pass
    elif idx is None:
        idxs = list(range(7, max_idx+1))
        #idxs = [22, 25, 26]
    else:
        idxs = [idx]
    print("idxs:", idxs)

    lookup = {}
    for idx in idxs:
        print()
        geometry = Geometry(key, idx, True)
        #graph = geometry.build_graph(desc)
        G = geometry.G
        print("|G| = %d, idx = %d" % (len(G), idx))
        lookup[idx] = G
        #continue
    
        faces = geometry.get_cosets([0,1,1])
        edges = geometry.get_cosets([1,0,1])
        verts = geometry.get_cosets([1,1,0])
        print("faces=%d, edges=%d, verts=%d"%(len(faces), len(edges), len(verts)))
    
        A = get_adj(faces, verts)
        #from bruhat.hecke import colour
        #colour(A)
        print(shortstr(A))
        print()
    
        G = geometry.G
        gens = G.gens
        f, e, v = gens # face, edge, vert
        H = Group.generate([f*e, f*v, e*v])
        assert len(G) % len(H) == 0
        index = len(G) // len(H)
        assert index in [1, 2]
        if index == 1:
            print("non-orientable")
        elif index == 2:
            print("orientable")
    
        g = f*e*v*e*f # colour symmetry
        H = Group.generate([e, v, g])
        print("|H| =", len(H))
        index = len(G)//len(H)
        print("|G:H| =", index)

        if index != 3:
            continue

        if len(G) > 1000:
            continue
        colours = G.left_cosets(H)
        B = get_adj(faces, colours)
        print(shortstr(B))
        print("colour:", B.sum(0), B.sum(1))

    #Gs = list(lookup.values())
    #for G1 in Gs:
    #  for G2 in Gs:
    for jdx in idxs:
      for idx in idxs:
        G1 = lookup[idx]
        G2 = lookup[jdx]
        if len(G2) < len(G1) and len(G1)%len(G2)==0:
            hom = close_hom(zip(G1.gen, G2.gen))
            if hom is not None:
                print("|G_%d| = %d"%(jdx, len(G2)), "subgroup of |G_%d| = %d"%(idx, len(G1)))


def unwrap_colour():
    print("unwrap_colour")
    key = (3, 8)
    idx0 = argv.get("idx0", 11)
    idx = argv.get("idx", idx0+1)
    N = len(lins_db.db[key])
    print("N =", N)

    G0 = Geometry(key, idx0, True).G
    print("|G_%d| = %d" % (idx0, len(G0)))

    idx -= 1
    while idx+1 < N:
        idx += 1
        geometry = Geometry(key, idx, True)
        G = geometry.G
        print("|G_%d| = %d" % (idx, len(G)))
        if len(G) % len(G0) != 0:
            continue
    
        hom = close_hom(zip(G.gen, G0.gen))
        if hom is not None:
            print("G_%d is subgroup of G_%d"%(idx0, idx))
            break



def make_genons():

    #solve.int_scalar = numpy.int64
    #print(solve.int_scalar)
    import qupy.ldpc.solve
    solve.int_scalar = qupy.ldpc.solve.int_scalar
    from qupy.ldpc.css import CSSCode

    key = (3, 6)
    idx = argv.get("idx", 12)
    geometry = Geometry(key, idx, True)
    #graph = geometry.build_graph(desc)
    G = geometry.G
    print("|G| = %d, idx = %d" % (len(G), idx))

    faces = geometry.get_cosets([0,1,1])
    edges = geometry.get_cosets([1,0,1])
    verts = geometry.get_cosets([1,1,0])
    print("faces=%d, edges=%d, verts=%d"%(len(faces), len(edges), len(verts)))

    A = get_adj(faces, verts)
    Az = get_adj(faces, edges)
    Ax = get_adj(verts, edges)
    print(A.shape, Az.shape, Ax.shape)
    print(rank(Az))
    print(rank(Ax))

    code = CSSCode(Hx=A, Hz=A)
    print(code)

    code = CSSCode(Hx=Ax, Hz=Az)
    print(code)
    print(shortstrx(Ax, Az.transpose()))
    #print(shortstr(code.Lx))
    #print()
    #print(shortstr(code.Lz))

    return

    Lx = code.Lx
    Hx = code.Hx
    k, n = Lx.shape
    d = n
    while d > 4:
        best = None
        for u in Lx:
            u = monte_carlo_css(Hx, u)
            d1 = u.sum()
            if d1 < d:
                best = u
                d = d1
    print("d =", d)
    print(best) # edges
    assert dot2(Az, best).sum() == 0

    print("Ax:", Ax.shape)

    print("best:")
    for i, x in enumerate(best):
        if x:
            print("\tedge:", i)

    v = zeros2(Ax.shape[0])
    #for i in range(mx):
    for vert, h in enumerate(Ax):
        hh = h*best
        if hh.sum() == 0:
            continue
        assert hh.sum() == 1
        edge = numpy.where(hh)[0][0]
        print("edge %d hits vertex %d" % (edge, vert))
        v[vert] = 1
    print(v)
    print(dot2(A, v))
    


def build_hyperbolic_4():
    key = (4, 3, 3, 5)
    idx = argv.get("idx", 5) # < 5 is trivial...
    geometry = Geometry(key, idx, False)
    graphs = []
    for desc in [
        #[0,1,1,1,1],
        #[1,0,1,1,1],
        [1,1,0,1,1],
        [1,1,1,0,1],
        #[1,1,1,1,0],
    ]:
        graph = geometry.build_graph(desc)
        print(desc, len(graph))
        graphs.append(graph)

    left, right = graphs[0], graphs[1]
    lhs, rhs = left.get_words(), right.get_words()
    #for l in lhs:
    #    right.follow_path(0, l)
    lookup = right.mulclose(left.hgens)
    print(len(lookup))
    lookup = left.mulclose(right.hgens)
    print(len(lookup))


def build_bicolour():

    key = (6,4)
    if 0:
        # search...
        idx = 0
        while 1:
            try:
                geometry = Geometry(key, idx)
            except IndexError:
                break
    
            G = geometry.G
        
            gens = G.gens
            a, b, c = gens
            i, j = ((a*b).order(), (b*c).order())
            if i==3:
                print("idx = %d, |G| = %d" % (idx, len(G)))
            idx += 1

    idx = 27
    geometry = Geometry(key, idx)
    G = geometry.G
    gens = G.gens
    a, b, c = gens
    i, j = ((a*b).order(), (b*c).order())

    print(i, j)



def build_group(idx=0, halve=False):
    # start with hyperbolic Coxeter reflection group: a--5--b--5--c
    # Then we add another generator "d" that halves these triangles.

    ngens = 3
    a, b, c, = range(ngens)
    rels_552 = [
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
    ]

    if halve:
        d = ngens
        ngens += 1
        rels_552 += [ (d,)*2, (c,d,a,d), (a,d)*4, (b,d)*2, (d,c)*4, ]
        # Now we get:
        # d--4--a--5--b or,
        # b--5--c--4--d
        # where a equals c.

    rels = [
        (0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1),
        (0, 1, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1,
            0, 2, 0, 1),
        (1, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1,
            0, 1, 0, 1, 0, 2),
        (0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 1, 2, 1, 2, 0, 1,
            0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 1, 0, 2, 0, 1,
            0, 2, 1, 2, 0, 2, 1, 2),
        (0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 1,
            0, 2, 1, 2, 1, 2, 1, 2),
        (1, 2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 1, 2,
            0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1),
        (1, 2, 1, 2, 0, 1, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2,
            1, 2, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 0, 1,
            0, 1, 0, 2, 1, 2, 0, 2, 1, 2, 0, 1, 0, 1, 0, 1),
    ]

    for rel in rels:
        assert len(rel)%2 == 0

    rel = rels[idx]
    graph = Schreier(ngens, rels_552 + [rel])
    graph.build()
    G = graph.get_group()
    return G

gap_54 = """
LoadPackage("recog");
LoadPackage("LINS");
F := FreeGroup("a", "b", "c");;
AssignGeneratorVariables(F);;
G := F/[a^2,b^2,c^2,(a*b)^5,(a*c)^2,(b*c)^4];;
gr:=LowIndexNormalSubgroupsSearchForAll(G, 600);
L := List(gr);
S := Grp(L[9]);
GeneratorsOfGroup(S);
"""


gap_344 = """
LoadPackage("recog");
LoadPackage("LINS");
F := FreeGroup("a", "b", "c", "d");;
AssignGeneratorVariables(F);;
G := F/[a^2,b^2,c^2,d^2,(a*b)^3,(a*c)^2,(a*d)^2,(b*c)^4, (b*d)^2, (c*d)^4];;
gr:=LowIndexNormalSubgroupsSearchForAll(G, 600);
L := List(gr);
S := Grp(L[35]);
GeneratorsOfGroup(S);
"""

def gap_parse(s, a, b, c, d=None):
    s = s.replace("^-1", "")
    s = s.replace("*", "+")
    s = s.replace("^", "*")
    s = s.replace(" ", "")
    return eval(s)


def build_group_524():
    ngens = 3
    a, b, c, = range(ngens)
    rels_524 = [
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*4, (a,c)*2,
    ]
    a, b, c = (a,), (b,), (c,)

    relss = []

    # order 160, non-checkerboard
    rels = [ 
        a+(b+a+c)*2+b+(c+a+b)*2+c, (a+c+b)*2+a+c+(b+c+a)*2+b, 
        (b+c+b+a)*2+(b+c+b+a)*2, (c+b+a+b)*2+(c+b+a+b)*2, 
        b+a+(b+a+c)*2+b+(c+a+b)*2+c+b, b+c+(b+a)*2+c+b+a+c+b+c+(a+b)*2+c,
    ]
    assert rels == gap_parse("""
    [ a*(b*a*c)^2*b*(c^-1*a^-1*b^-1)^2*c^-1, (a*c*b)^2*a*c*(b^-1*c^-1*a^-1)^2*b^-1,
    (b*c*b*a)^2*(b^-1*c^-1*b^-1*a^-1)^2, (c*b*a*b)^2*(c^-1*b^-1*a^-1*b^-1)^2,
    b*a*(b*a*c)^2*b*(c^-1*a^-1*b^-1)^2*c^-1*b^-1,
    b*c*(b*a)^2*c*b*a*c^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1 ]
    """, a, b, c)
    relss.append(rels)

    # order 240, checkerboard, [[30,8,3]] Bring's curve
    relss.append(gap_parse("""
    [ (b*a*c)^3*(b^-1*c^-1*a^-1)^3, (c*b*a)^3*c^-1*(b^-1*c^-1*a^-1)^2*b^-1*a^-1, 
    ((b*a)^2*c)^2*(b^-1*a^-1*b^-1*c^-1*a^-1)^2, 
    b*a*b*c*(b*a)^2*c*b*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*b^-1*a^-1, 
    (b*a*c*b*a)^2*(b^-1*c^-1*a^-1*b^-1*a^-1)^2, 
    b*c*(b*a)^2*c*b*a*b*c^-1*(a^-1*b^-1)^2*c^-1*(b^-1*a^-1)^2, 
    b*(c*b*a)^3*c^-1*(b^-1*c^-1*a^-1)^2*b^-1*a^-1*b^-1, 
    a*b*(c*b*a)^3*c^-1*(b^-1*c^-1*a^-1)^2*(b^-1*a^-1)^2 ]
    """, a, b, c))

    idx = argv.get("idx", 0)
    graph = Schreier(ngens, rels_524 + relss[idx])
    graph.build(maxsize=100000)
    n = len(graph)
    print("|G| =", n)
    G = graph.get_group()
    return G


def build_group_624():
    ngens = 3
    a, b, c, = range(ngens)
    rels_624 = [
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*6, (b,c)*4, (a,c)*2,
    ]
    a, b, c = (a,), (b,), (c,)

    relss = []

    rels = gap_parse("""
        [b*c*b*a*b*c*b^-1*a^-1*b^-1*c^-1*b^-1*a^-1, 
        c*b*a*b*c*b*a^-1*b^-1*c^-1*b^-1*a^-1*b^-1,
        a*c*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2, 
        b*a*c*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1,
        a*b*a*c*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^3,
        (b*a)^2*c*b*a*b*c*b^-1*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*a^-1,
        (c*b*a)^2*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1,
        (a*c*b)^2*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1,
        b*(c*b*a)^2*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*b^-1,
        (b*a*c)^2*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1*b^-1,
        a*(b*a*c)^2*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1*b^-1*a^-1,
        c*b*a*(b*a*c)^2*b*a*b*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*(a^-1*b^-1)^2 ]
    """, a, b, c)
    relss.append(rels) # [[30,11,3]]

    rels = gap_parse("""
[ ((b*a)^2*c)^2*(b^-1*a^-1*b^-1*c^-1*a^-1)^2, (b*a*c*b*a)^2*(b^-1*c^-1*a^-1*b^-1*a^-1)^2,
  (c*(b*a)^2)^2*c^-1*b^-1*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*a^-1,
  (b*a)^2*b*c*(b*a)^2*c*b^-1*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*b^-1*a^-1,
  b*a*b*c*(b*a)^2*c*b*a*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*(b^-1*a^-1)^2,
  (b*a*b*c)^2*b*a*c*b^-1*c^-1*(a^-1*b^-1*c^-1*b^-1)^2*a^-1,
  b*c*(b*a)^2*b*c*b*a*b*c^-1*(b^-1*c^-1*a^-1)^2*b^-1*c^-1*b^-1*a^-1,
  b*(c*(b*a)^2)^2*(c^-1*b^-1*a^-1*b^-1)^2*a^-1*b^-1*a^-1,
  (b*c*b*a)^2*b*a*b*c^-1*b^-1*c^-1*a^-1*b^-1*c^-1*b^-1*a^-1*b^-1*c^-1*a^-1,
  (b*c*b*a)^2*c*b*a*c^-1*(b^-1*c^-1*b^-1*a^-1)^2*b^-1*a^-1,
  b*(c*b*a)^2*b*c*b*a*(c^-1*b^-1)^2*(a^-1*b^-1)^2*c^-1*b^-1*a^-1,
  c*b*a*(b*a*b*c)^2*(b^-1*c^-1*a^-1)^2*b^-1*c^-1*b^-1*a^-1*b^-1,
  (c*b*a*b)^2*a*b*c*b^-1*c^-1*a^-1*b^-1*c^-1*b^-1*a^-1*b^-1*c^-1*a^-1*b^-1,
  (c*b*a*b)^2*c*b*a*(c^-1*b^-1*a^-1*b^-1)^2*c^-1*b^-1*a^-1,
  c*(b*a*c*b*a)^2*b^-1*c^-1*a^-1*b^-1*(a^-1*b^-1*c^-1)^2*a^-1,
  c*b*a*(c*b*a*b)^2*c^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*b^-1*a^-1*b^-1,
  a*c*b*a*(b*a*b*c)^2*(b^-1*c^-1*a^-1)^2*b^-1*c^-1*(b^-1*a^-1)^2,
  (b*a*c)^2*b*a*b*c*b*a*(b^-1*c^-1*a^-1)^2*b^-1*a^-1*b^-1*c^-1*b^-1*a^-1,
  b*c*b*a*(b*a*c)^2*b*a*b^-1*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*a^-1*b^-1*a^-1,
  c*b*a*(b*a*c)^2*b*a*b*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*(a^-1*b^-1)^2,
  (c*b*a)^3*b*a*c*b^-1*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1*b^-1,
  a*c*b*a*(b*a*c)^2*b*a*b*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*(a^-1*b^-1)^2*a^-1 ]
    """, a, b, c)
    relss.append(rels) # [[30,11,3]]

    idx = argv.get("idx", 0)
    graph = Schreier(ngens, rels_624 + relss[idx])
    graph.build(maxsize=100000)
    n = len(graph)
    print("|G| =", n)
    G = graph.get_group()
    return G


def test_prism():
    H = solve.parse("""
    111..1..
    1.11....
    11.11...
    .1..11..
    ..11..1.
    ...11.11
    ....11.1
    ..1..111
    """)
    code = find_code(H, False, True, True)
    print(code)
    n,k,d = code.get_params()
    print("[[%s,%s,%s]]"%(n,k,d))


def test_ldpc():
    n = argv.get("n", 10)
    m = n-1
    rw = argv.get("rw", 4) # row weight
    cw = argv.get("cw", 2) # minimum column weight
    distance = argv.get("distance", 2) # minimum distance
    yop = argv.get("yop", True)
    css = False

    while 1:
        while 1:
            H = zeros2(m, n)
            while 1:
                for i in range(m):
                    H[i] = 0
                    while H[i].sum() < rw:
                        j = randint(0, n-1)
                        H[i, j] = 1
                cols = H.sum(0)
                if cols.min() >= cw:
                    break
            code = find_code(H, css, yop)
            if code is not None:
                break
            print(".", end="", flush=True)
    
        n, k, d = code.get_params()
        print("[[%s,%s,%s]]"%(code.get_params()), end="", flush=True)
        if d >= distance:
            break
    print()
    print(code)


def test_surface():

    sides = argv.get("sides", 5)
    if sides==5:
        G = build_group_524()
    elif sides==6:
        G = build_group_624()
    else:
        assert 0

    if argv.dot:
        graph = Schreier(ngens, rels_524 + relss[idx])
        #graph.build([b+c+b+c]) # halve a vert
        graph.build([a+c]) # halve an edge
        make_dot(graph, "neato")
        return

    a, b, c = G.gens
    red, blue, green = G.gens

#    count = 0
#    for g in G:
#        if g.order() == 2:
#            count += 1
#    print("order 2 elements:", count)

    if 0:
        assert g.order() == 2
        H = [g, g*g]
        pairs = G.left_cosets(H)


    if argv.fold_vert:
        g = blue*green*blue*green # divide a vert in half
        assert g.order() == 2
    elif argv.fold_edge:
        g = red*green # divide an edge in half
        assert g.order() == 2
    else:
        g = G.identity

    H = Group.generate([g])
    #pairs = G.left_cosets(H)
    act = G.action_subgroup(H)
    faces = act.orbits(Group.generate([blue, red]))
    print("faces:", len(faces))
    print([len(f) for f in faces])
    edges = act.orbits(Group.generate([green, red]))
    print("edges:", len(edges))
    print([len(f) for f in edges])
    verts = act.orbits(Group.generate([green, blue]))
    print("verts:", len(verts))
    print([len(f) for f in verts])

#    def setkey(items):
#        items = list(items)
#        items.sort(key=str)
#        return str(items)
#
#    def canonical(items):
#        # XXX FAIL: items is a list of set's
#        items = list(items)
#        items.sort(key=setkey)
#        return items
#
#    faces = canonical(faces)
#    edges = canonical(edges)
#    verts = canonical(verts)

    A = get_adj(faces, edges)
    B = get_adj(verts, edges)

    AA = dot(A, A.transpose())
    BB = dot(B, B.transpose())

    print("face -- face")
    print(AA.shape)
    print(AA)

    AB = dot(A, B.transpose())
    print("face -- vert")
    print(AB.shape)
    print(AB)

    if argv.find_code:
        code = find_code(AB, argv.css, argv.yop)
        if code is None:
            return
        print(code)
        n,k,d = code.get_params()
        print("[[%s,%s,%s]]"%(n,k,d))
        if d is None:
            print(shortstr(code.get_logops()))
        return

    labels = string.ascii_lowercase

    def show_incident(A):
        print("graph")
        print("{")
        for (i,j) in numpy.ndindex(A.shape):
            if A[i,j] and i<j:
                print("\t", labels[i], "--", labels[j], ";")
        print("}")
    #print("vert -- edge")

#    print("face -- face")
#    show_incident(AA)
#    print("vert -- vert")
#    show_incident(BB)
#    print("face -- edge")
#    show_incident(A)

    if argv.autos:
        from bruhat.isomorph import Tanner, search
        U = Tanner.build2(A, B)
        V = Tanner.build2(A, B)
        perms = []
        for f in search(U, V):
            perms.append(f)
        print("autos:", len(perms))
    
    assert dot2(A, B.transpose()).sum() == 0
    Hx, Hz = A, B
    Hx = linear_independent(Hx)
    Hz = linear_independent(Hz)
    print(Hx.shape)
    print(shortstr(Hx))
    print()
    print(Hz.shape)
    print(shortstr(Hz))
    

    return

    # checkerboard subgroup
    H = mulclose([red, blue, blue*green*blue*green])
    is_checkerboard = 2*len(H)==len(G)
    print("is_checkerboard:", is_checkerboard)
    checkerboard = G.left_cosets(H)

    if 1:
        # orientation subgroup
        L = Group.generate([red*blue, blue*green])
        print("L:", len(L))
        orients = G.left_cosets(L)
        is_orientable = len(orients) == 2
        print("is_orientable:", is_orientable)
        orients = list(orients)

    # vert'ex subgroup
    H = mulclose([blue, green])
    verts = G.left_cosets(H)
    print("verts:", len(verts))
    n = len(verts) # qubits

    # face subgroup
    H = mulclose([red, blue])
    assert len(H) == 10
    faces = G.left_cosets(H)
    print("faces:", len(faces))

    # edge subgroup
    H = mulclose([green, red])
    assert len(H) == 4
    edges = G.left_cosets(H)
    print("edges:", len(edges))

    if not is_checkerboard:
        return

    mx = len(faces)//2 - 1
    mz = len(faces)//2 - 1
    print("k =", n-mx-mz)

#    for face in faces:
#        print(face)
#        for black in checkerboard:
#            print('\t', face.intersect(black))


    Hx = zeros2(mx, n)
    Hz = zeros2(mz, n)
    ix = iz = 0
    black, white = checkerboard
    black, white = [], []
    for face in faces:
        if face.intersect(checkerboard[0]):
            black.append(face)
        else:
            white.append(face)

    Hx = get_adj(black, verts)
    Hz = get_adj(white, verts)
    print(shortstr(Hx))
    print()
    print(shortstr(Hz))
    assert dot2(Hx, Hz.transpose()).sum() == 0



def build_group_344():

    ngens = 4
    a, b, c, d = range(ngens)
    a, b, c, d = (a,), (b,), (c,), (d,)
    rels_624 = [
        2*a, 2*b, 2*c, 2*d,
        (a+b)*3, (a+c)*2, (a+d)*2,
        (b+c)*4, (b+d)*2, (c+d)*4
    ]

    relss = []

    rels = gap_parse("""
    [ (c*b*d)^2*(c^-1*d^-1*b^-1)^2, (d*c*b)^2*d^-1*c^-1*d^-1*b^-1*c^-1*b^-1,
    a*(c*b*d)^2*(c^-1*d^-1*b^-1)^2*a^-1, a*(d*c*b)^2*d^-1*c^-1*d^-1*b^-1*c^-1*b^-1*a^-1,
    b*c*d*c*b*a*c*b^-1*c^-1*d^-1*c^-1*b^-1*c^-1*a^-1,
    c*d*c*b*a*c*b*c^-1*d^-1*c^-1*b^-1*c^-1*a^-1*b^-1,
    d*c*b*a*c*b*c*d^-1*c^-1*b^-1*c^-1*a^-1*b^-1*c^-1 ]
    """, a, b, c, d)

    graph = Schreier(ngens, rels_624 + rels)
    graph.build()
    n = len(graph)
    print("|G| =", n)
    G = graph.get_group()
    return G
    

def build_group_333():
    ngens = 3
    r, g, b = (0,), (1,), (2,)
    rels = [r*2, g*2, b*2, (r+b)*3, (r+g)*3, (b+g)*3]
    rels += [(b+r+g)*4]
    graph = Schreier(ngens, rels)
    graph.build(maxsize=100000)
    print(len(graph))


def make_hyperbolic():

    G = build_group(2)
    H = build_group(1)
    print("|G| =", len(G))
    print("|G.gens| =", len(G.gens))
    #a, b, c, d = G.gens

    count = 0
    for g in G:
        if g.order() != 2:
            continue
        for h in G:
            if g*h != h*g:
                break
        else:
            print(g)
            print()
            count += 1
    print(count)

    hom = dict(zip(G.gens, H.gens))
    hom = close_hom(hom)
    #assert is_hom(hom)

    kernel = []
    for g,h in hom.items():
        if h.is_identity():
            kernel.append(g)
    print(len(kernel))
    for g in kernel:
        if not g.is_identity():
            print(g)
        else:
            print("is_identity")


def make_kagome():

    N = 6
    lookup = {}
    
    cubes = [(2*i, 2*j, 2*k) for i in range(N) for j in range(N) for k in range(N)]
    for (x, y, z) in cubes:
        key = x+1, y, z
        lookup[key] = len(lookup)
        key = x, y+1, z
        lookup[key] = len(lookup)
        key = x, y, z+1
        lookup[key] = len(lookup)
    def get(x, y, z):
        return lookup[x%(2*N), y%(2*N), z%(2*N)]
    n = len(lookup)
    print("qubits:", n)
    Hx = []
    for (i, j, k) in cubes:
        op = [0]*n
        for idx in [
            get(i+1, j, k),
            get(i, j+1, k),
            get(i, j, k+1),
            get(i-1, j, k),
            get(i, j-1, k),
            get(i, j, k-1),
        ]:
            op[idx] = 1
        Hx.append(op)
    Hx = array2(Hx)
    #print("Hx:")
    #print(shortstr(Hx))
    #print(Hx.sum(0))
    Hx = linear_independent(Hx)
    print("Hx:", len(Hx))

    Hz = []

    # diamond faces
    for (i, j, k) in cubes:
        op = [0]*n
        for idx in [
            get(i+0, j+1, k+0),
            get(i+0, j+0, k+1),
            get(i+0, j+1, k+2),
            get(i+0, j+2, k+1),
        ]:
            op[idx] = 1
        Hz.append(op)

        op = [0]*n
        for idx in [
            get(i+1, j+0, k+0),
            get(i+0, j+0, k+1),
            get(i+1, j+0, k+2),
            get(i+2, j+0, k+1),
        ]:
            op[idx] = 1
        Hz.append(op)

        op = [0]*n
        for idx in [
            get(i+1, j+0, k+0),
            get(i+0, j+1, k+0),
            get(i+1, j+2, k+0),
            get(i+2, j+1, k+0),
        ]:
            op[idx] = 1
        Hz.append(op)

    for (i, j, k) in cubes:
        break
        op = [0]*n
        for idx in [
            get(i+0, j+2, k+1),
            get(i+0, j+1, k+0),
            get(i+1, j+0, k+0),
            get(i+2, j+0, k+1),
            get(i+2, j+1, k+2),
            get(i+1, j+2, k+2),
        ]:
            op[idx] = 1
        Hz.append(op)

        op = [0]*n
        for idx in [
            get(i+0, j+0, k+1),
            get(i+1, j+0, k+0),
            get(i+2, j+1, k+0),
            get(i+2, j+2, k+1),
            get(i+1, j+2, k+2),
            get(i+0, j+1, k+2),
        ]:
            op[idx] = 1
        Hz.append(op)

        op = [0]*n
        for idx in [
            get(i+0, j+1, k+0),
            get(i+0, j+0, k+1),
            get(i+1, j+0, k+2),
            get(i+2, j+1, k+2),
            get(i+2, j+2, k+1),
            get(i+1, j+2, k+0),
        ]:
            op[idx] = 1
        Hz.append(op)

        op = [0]*n
        for idx in [
            get(i+0, j+2, k+1),
            get(i+0, j+1, k+2),
            get(i+1, j+0, k+2),
            get(i+2, j+0, k+1),
            get(i+2, j+1, k+0),
            get(i+1, j+2, k+0),
        ]:
            op[idx] = 1
        Hz.append(op)
    Hz = array2(Hz)
    #print("Hz:", len(Hz))
    #print(shortstr(Hz))
    #print(Hz.sum(0))

    Hz = linear_independent(Hz)
    print("Hz:", len(Hz))

    assert dot2(Hz, Hx.transpose()).sum() == 0

    code = QCode.build_css(Hx, Hz)
    print(code.get_params())

#    L = code.get_logops()
#    #print("L:", L.sum(1))
#    #print(shortstr(code.get_logops()))
#    for l in code.get_logops():
#        items = []
#        for i in range(2*n):
#            if l[i] == 0:
#                continue
#            idx = i//2
#            for (key, value) in lookup.items():
#                if value == idx:
#                    break
#            else:
#                assert 0
#            items.append((idx, key, "XZ"[i%2]))
#        if len(items)<=4:
#            print(items)


def test_toric():
    Hx = parse("""
[[ 1  1  .  .  .  .  .  .  .  1  .  .  1  .  .  .  .  . ],
 [ .  1  1  .  .  .  .  .  .  .  1  .  .  1  .  .  .  . ],
 [ 1  .  1  .  .  .  .  .  .  .  .  1  .  .  1  .  .  . ],
 [ .  .  .  1  1  .  .  .  .  .  .  .  1  .  .  1  .  . ],
 [ .  .  .  .  1  1  .  .  .  .  .  .  .  1  .  .  1  . ],
 [ .  .  .  1  .  1  .  .  .  .  .  .  .  .  1  .  .  1 ],
 [ .  .  .  .  .  .  1  1  .  1  .  .  .  .  .  1  .  . ],
 [ .  .  .  .  .  .  .  1  1  .  1  .  .  .  .  .  1  . ],
 [ .  .  .  .  .  .  1  .  1  .  .  1  .  .  .  .  .  1 ]]
    """)
    Hz = parse("""
[[ 1  .  .  1  .  .  .  .  . ],
 [ .  1  .  .  1  .  .  .  . ],
 [ .  .  1  .  .  1  .  .  . ],
 [ .  .  .  1  .  .  1  .  . ],
 [ .  .  .  .  1  .  .  1  . ],
 [ .  .  .  .  .  1  .  .  1 ],
 [ 1  .  .  .  .  .  1  .  . ],
 [ .  1  .  .  .  .  .  1  . ],
 [ .  .  1  .  .  .  .  .  1 ],
 [ 1  1  .  .  .  .  .  .  . ],
 [ .  1  1  .  .  .  .  .  . ],
 [ 1  .  1  .  .  .  .  .  . ],
 [ .  .  .  1  1  .  .  .  . ],
 [ .  .  .  .  1  1  .  .  . ],
 [ .  .  .  1  .  1  .  .  . ],
 [ .  .  .  .  .  .  1  1  . ],
 [ .  .  .  .  .  .  .  1  1 ],
 [ .  .  .  .  .  .  1  .  1 ]]
    """).transpose()

    code = QCode.build_css(Hx, Hz)
    L = code.get_logops()
    print(L)
    print(code.get_params())
    

def test_color_24():

    n = 24

    R = parse("""
    ....11..11..11..........
    ......11..11..11........
    11..............11..11..
    ..11..............11..11
    """)
    G = parse("""
    11...11..11.............
    ..111..11..1............
    ............11...11..11.
    ..............111..11..1
    """)
    B = parse("""
    1..111..............1..1
    .11...11.............11.
    ........1..11..111......
    .........11..11...11....
    """)
    for A in [R, G, B]:
      for B in [R, G, B]:
        assert dot2(A, B.transpose()).sum() == 0

    H = numpy.concatenate((R, G, B), axis=0)
    assert H.shape[1] == n
    code = QCode.build_css(H, H)
    print(code.get_params())
    L = code.get_logops()
    print(shortstr(L))


def test_color():

    lattice = {}
    width = 6
    height = 6

    assert width % 3 == 0 
    assert height % 3 == 0 

    pt = lambda x, y : (x%width, y%height)
    def pt(x,y):
        i, j = (x%width, y%height)
        if (i,j) not in lattice:
            lattice[i,j] = len(lattice)
        idx = lattice[i, j]
        return idx

    hexagons = {}
    for x in range(width):
      for y in range(height):
        if x%3 != y%3:
            continue
        # clockwise:
        h = [
            pt(x+1,y),
            pt(x+1,y-1),
            pt(x,y-1),
            pt(x-1,y),
            pt(x-1,y+1),
            pt(x,y+1),
        ]
        assert len(set(h)) == 6
        hexagons[x,y] = h

    assert len(hexagons) == width*height // 3
    n = len(lattice)

    print("n =", n)

    if 0:
        H = []
        for h in hexagons.values():
            v = zeros2(n)
            v[h] = 1
            H.append(v)

    elif 0:
        # XZXZXZ version of the color code...
        H = []
        for h in hexagons.values():
            v = zeros2(n, 2)
            for i in h:
                v[i, i%2] = 1
            H.append(v)
            v = zeros2(n, 2)
            for i in h:
                v[i, 1-i%2] = 1
            H.append(v)
    
    elif 1:
        swap = [0, 1, 1, 1, 0, 0]
        H = []
        for key, h in hexagons.items():
            if key in [ (1,1), (3,0), (1,4), (3,3), (5,2), (5,5), ]:
                idxs = swap
            else:
                idxs = [0]*6
            v = zeros2(n, 2)
            for ii, i in enumerate(h):
                v[i, idxs[ii]] = 1
            H.append(v)
            v = zeros2(n, 2)
            for ii, i in enumerate(h):
                v[i, 1-idxs[ii]] = 1
            H.append(v)

    elif 0:
        swap00 = [0, 0, 0, 0, 0, 0]
        swap10 = [0, 1, 1, 1, 0, 0]
        swap01 = [0, 0, 1, 1, 1, 0]
        swap11 = [0, 1, 0, 0, 1, 0]
        H = []
        for key, h in hexagons.items():
            if key in [(1,4), (3,0), (5,2)]:
                idxs = swap10
            elif key in [(0,3), (2,5), (4,1)]:
                idxs = swap01
            elif key in [(1,1), (3,3), (5,5)]:
                idxs = swap11
            elif key in [(0,0), (2,2), (4,4)]:
                idxs = [0]*6
            else:
                assert 0, key
            v = zeros2(n, 2)
            for ii, i in enumerate(h):
                v[i, idxs[ii]] = 1
            H.append(v)
            v = zeros2(n, 2)
            for ii, i in enumerate(h):
                v[i, 1-idxs[ii]] = 1
            H.append(v)

    H = array2(H)
    #print(shortstr(H))
    code = QCode(H)
    #print(shortstr(code.H))
    print(code.get_params())

    L = code.get_logops()
    print(shortstr(L))

    print()
    ws = []
    for l in span(L):
        if l.sum() == 0:
            continue
        w = monte_carlo(code.H, l)
        ws.append(w)
    print(ws)
    print(ws.count(4))



if __name__ == "__main__":

    from time import time

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


