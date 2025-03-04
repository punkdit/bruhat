#!/usr/bin/env python3

"""
Algebraic groups: matrix groups over Z/pZ.

see also:
    quantale.py
    combinatorial.py
    orthogonal.py

"""


import sys, os
import random
from random import randint, choice
from functools import reduce
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import add, mul
from math import prod

import numpy

scalar = numpy.int64
#scalar = numpy.int8 # CAREFUL !! only works for p < 11

from bruhat.action import mulclose, mulclose_hom
from bruhat.spec import isprime
from bruhat.argv import argv
from bruhat.solve import parse, enum2, row_reduce, span, shortstr, rank, shortstrx, pseudo_inverse, intersect
from bruhat.solve import zeros2, identity2, solve
from bruhat.dev import geometry
from bruhat.util import cross, allperms, choose
from bruhat.smap import SMap
from bruhat.equ import Equ

EPSILON = 1e-8

DEFAULT_P = argv.get("p", 2)


def qchoose(n, m, p=DEFAULT_P):
    if m>n:
        return
    col = m
    row = n-m
    for A in geometry.get_cell(row, col, p):
        yield A

qchoose_2 = qchoose # backwards compat


def pdiv(i, j, p=DEFAULT_P):
    " i/j mod p "
    assert j!=0
    for k in range(1, p):
        if (j*k)%p == i:
            return k
    assert 0

def pinv(i, p=DEFAULT_P):
    return pdiv(1, i, p)


def swap_row(A, j, k):
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row

def swap_col(A, j, k):
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col

def row_reduce_p(A, p=DEFAULT_P, truncate=False, inplace=False, check=False, verbose=False):
    """ Remove zero rows if truncate==True
    """

    zero = 0
    one = 1

    assert len(A.shape)==2, A.shape
    m, n = A.shape
    if not inplace:
        A = A.copy()

    if m*n==0:
        if truncate and m:
            A = A[:0, :]
        return A

    if verbose:
        print("row_reduce")
        #print("%d rows, %d cols" % (m, n))

    i = 0
    j = 0
    while i < m and j < n:
        if verbose:
            print("i, j = %d, %d" % (i, j))
            print("A:")
            print(A)

        assert i<=j
        if i and check:
            assert (A[i:,:j]!=0).sum() == 0

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if A[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print("swap", i, i1)
            swap_row(A, i, i1)

        assert A[i, j] != zero
        for i1 in range(i+1, m):
            if A[i1, j]:
                r = p-pdiv(A[i1, j], A[i, j], p)
                if verbose:
                    print("add %d times row %s to %s" % (r, i, i1))
                A[i1, :] += r*A[i, :]
                A %= p
                assert A[i1, j] == zero

        i += 1
        j += 1

    if truncate:
        m = A.shape[0]-1
        #print("sum:", m, A[m, :], A[m, :].sum())
        while m>=0 and (A[m, :]!=0).sum()==0:
            m -= 1
        A = A[:m+1, :]

    if verbose:
        print()

    return A


# see orthogonal.py for numba version
def normal_form_p(A, p=DEFAULT_P, truncate=True):
    "reduced row-echelon form"
    #print("normal_form")
    #print(A)
    A = row_reduce_p(A, truncate=truncate)
    #print(A)
    m, n = A.shape
    j = 0
    for i in range(m):
        while j<n and A[i, j] == 0:
            j += 1
        if j==n:
            break
        r = A[i, j]
        if r != 1:
            A[i, :] = pinv(r, p) * A[i, :]
            A %= p
        assert A[i, j] == 1
        i0 = i-1
        while i0>=0:
            r = A[i0, j]
            if r!=0:
                A[i0, :] += (p-r)*A[i, :]
                A %= p
            i0 -= 1
        j += 1
    #print(A)
    return A

def test_row_reduce():
    p = 3
    m, n = 4, 4
    A = numpy.zeros((m, n))
    for idx in numpy.ndindex(A.shape):
        A[idx] = random.randint(0, p-1)
    B = row_reduce_p(A, p=p, verbose=False)
    C = normal_form_p(B, p)
    print(A)
    print(B)
    print(C)

_cache = {}
def normal_form(A, p=DEFAULT_P):
    "reduced row-echelon form"
    if p!=2:
        return normal_form_p(A, p)
    key = A.tobytes()
    if key in _cache:
        return _cache[key]
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
            r = A[i0, j]
            if r!=0:
                A[i0, :] += A[i, :]
                A %= p
            i0 -= 1
        j += 1
    #print(A)
    _cache[key] = A
    return A



def all_matrices(m, n, p=DEFAULT_P):
    shape = ((p,)*m*n)
    for idxs in numpy.ndindex(shape):
        M = numpy.array(idxs)
        M.shape = (m, n)
        yield M

def all_codes(m, n, p=DEFAULT_P):
    assert p==2
    for m1 in range(m+1):
        for M1 in geometry.all_codes(m1, n):
            M = numpy.zeros((m, n), dtype=scalar)
            M[:m1, :] = M1
            yield M


class Matrix(object):
    def __init__(self, A, p=DEFAULT_P, shape=None, name="?"):
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
        self.key = (self.p, self.A.tobytes())
        self._hash = hash(self.key)
        self.shape = A.shape
        self.name = name

    @classmethod
    def perm(cls, items, p=DEFAULT_P, name="?"):
        n = len(items)
        A = numpy.zeros((n, n), dtype=scalar)
        for i, ii in enumerate(items):
            A[ii, i] = 1
        return Matrix(A, p, name=name)

    @classmethod
    def identity(cls, n, p=DEFAULT_P):
        A = numpy.identity(n, dtype=scalar)
        return Matrix(A, p, name="I")

    def get_bitkey(self):
        assert self.p == 2
        A = numpy.packbits(self.A)
        return A.tobytes()

    def __str__(self):
        return str(self.A)

    def __repr__(self):
        return "Matrix(%s)"%str(self.A)

    def shortstr(self):
        return shortstr(self.A)

    def latex(self, zero=".", delim="$"):
        rows = ["&".join((str(i) if i else zero) for i in row) for row in self.A]
        rows = r"\\".join(rows)
        rows = r"\begin{bmatrix}%s\end{bmatrix}"%rows
        rows = delim+rows+delim
        return rows

    def __hash__(self):
        return self._hash

    def is_zero(self):
        return self.A.sum() == 0

    def is_upper_triangular(self):
        m, n = self.shape
        A = self.A
        for i in range(1,m):
          for j in range(i):
            if A[i,j]:
                return False
        return True

    def __len__(self):
        return len(self.A)

    def __eq__(self, other):
        assert self.p == other.p
        return self.key == other.key

    def __ne__(self, other):
        assert self.p == other.p
        return self.key != other.key

    def __lt__(self, other):
        assert self.p == other.p
        return self.key < other.key

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
            return Matrix(A, self.p, name=self.name+other.name)
        else:
            return NotImplemented

    def __rmul__(self, r):
        A = r*self.A
        return Matrix(A, self.p)

    def __getitem__(self, idx):
        A = self.A[idx]
        #print("__getitem__", idx, type(A))
        if type(A) is scalar:
            return A
        return Matrix(A, self.p)

    def transpose(self):
        return Matrix(self.A.transpose(), self.p)

    @property
    def t(self):
        return self.transpose()

    def sum(self):
        return self.A.sum()

    def mask(self, A):
        return Matrix(self.A * A, self.p) # pointwise multiply !

    def min(self, *arg, **kw):
        return self.A.min(*arg, **kw)

    def max(self, *arg, **kw):
        return self.A.max(*arg, **kw)

    def normal_form(self, cols=None):
        A = self.A
        m, n = A.shape
        if cols is not None:
            assert len(cols)==n
            A0 = A
            A = A[:, cols]
            inv = [cols.index(i) for i in range(n)]
            B = A[:, inv]
            #print()
            #print(A0)
            #print(A, cols)
            #print(B, inv)
            assert str(A0) == str(B)
        A = normal_form(A, self.p)
        if cols is not None:
            A = A[:, inv]
        return Matrix(A, self.p)

    @classmethod
    def all_codes(cls, m, n, p=DEFAULT_P):
        assert p==2
        for A in geometry.all_codes(m, n):
            yield cls(A, p)

    @classmethod
    def qchoose(cls, m, n, p=DEFAULT_P):
        for A in qchoose(m, n, p):
            yield cls(A, p)

    def solve(self, other):
        assert self.p == 2, "not implemented"
        A = solve(self.A, other.A)
        return Matrix(A, self.p)

    def get_pivots(M):
        A = M.A
        m, n = A.shape
        pivots = []
        col = 0
        for row in range(m):
            while col < n:
                if A[row, col]:
                    pivots.append(col)
                    col += 1
                    break
                col += 1
        return tuple(pivots)

    def inverse(self):
        assert self.p == 2
        B = pseudo_inverse(self.A)
        return Matrix(B, self.p)

    def span(self):
        V = []
        assert self.p == 2
        for v in span(self.A):
            v = Matrix(v, self.p)
            V.append(v)
        return V

    def order(self):
        n = len(self)
        I = Matrix.identity(n)
        count = 1
        g = self
        while g != I:
            g = self*g
            count += 1
        return count

    def direct_sum(self, other):
        A, B = self.A, other.A
        AB = numpy.zeros((A.shape[0]+B.shape[0], A.shape[1]+B.shape[1]), scalar)
        m,n = A.shape
        AB[:m,:n] = A
        AB[m:,n:] = B
        return Matrix(AB, self.p)


def test_matrix():
    M = Matrix([[1,1,0],[0,1,0]])
    M1 = Matrix([[1,0,0],[0,1,0]])
    assert M.normal_form() == M1

    M = Matrix([[1,2,0],[0,1,2]], p=3)
    M1 = Matrix([[1,0,2],[0,1,0]], p=3)
    assert M.normal_form() == M1


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


class Algebraic(object):
    def __init__(self, gen, order=None, p=DEFAULT_P, G=None, verbose=False, **kw):
        self.__dict__.update(kw)
        self.gens = self.gen = list(gen) # XXX rename as gens
        if G is not None:
            assert order is None or order==len(G)
            order = len(G)
        self.order = order
        self.G = G # elements
        self.p = p
        assert gen
        A = gen[0]
        self.n = len(A)
        assert p == A.p
        self.I = Matrix.identity(self.n, p)
        self.verbose = verbose

    def __str__(self):
        return "%s(group of order %d)"%(self.__class__.__name__, len(self))

    def get_elements(self):
        if self.G is None:
            I = self.I
            G = mulclose(self.gen, maxsize=self.order, verbose=self.verbose)
            G.remove(I)
            G.add(I)
            G = list(G)
            self.G = G
            self.order = len(self.G)
        return self.G

    def sample(self):
        gen = self.gen
        N = 10*len(gen)
        A = choice(gen)
        for i in range(N):
            A = A*choice(gen)
        return A

    def __len__(self):
        if self.order is None:
            self.get_elements()
        return self.order

    def __getitem__(self, idx):
        if self.G is None:
            self.get_elements()
        return self.G[idx]

    def __contains__(self, g):
        return g in self.get_elements()

    def __eq__(self, other):
        return set(self.get_elements()) == set(other.get_elements())

    def left_stabilizer(self, M):
        # find subgroup that stabilize the rowspace of M
        V = M.span()
        H = []
        for g in self:
            for v in V:
                if g*v not in V:
                    break
            else:
                H.append(g)
        return Algebraic(H)
    
    def right_stabilizer(self, M):
        # find subgroup that stabilize the rowspace of M
        V = M.span()
        H = []
        for g in self:
            for v in V:
                if v*g not in V:
                    break
            else:
                H.append(g)
        return Algebraic(H)
    
    def show(self):
        items = [M.A for M in self]
        M = numpy.array(items)
    
        smap = SMap()
        for (i, j) in numpy.ndindex(M.shape[1:]):
            if numpy.alltrue(M[:, i, j] == 0):
                smap[i,j] = "."
            elif numpy.alltrue(M[:, i, j] == 1):
                smap[i,j] = "1"
            else:
                smap[i,j] = "*"
        return str(smap)

    @classmethod
    def SL(cls, n, p=DEFAULT_P, **kw):
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
        return cls(gen, order, p=p, **kw)
    
    @classmethod
    def GL(cls, n, p=DEFAULT_P, **kw):
        "general linear group"
        assert int(n)==n
        assert int(p)==p
        assert n>0
        assert isprime(p)
    
        if n>1:
            H = cls.SL(n, p)
            gen = list(H.gen)
        else:
            gen = []
        if p>2:
            # find generator of GL(1,F_p)
            for a in range(1,p):
                items = set((a**i)%p for i in range(p))
                #print(items, len(items))
                if len(items) == p-1:
                    break
            else:
                assert 0
            for i in range(n):
                A = numpy.identity(n, scalar)
                A[i,i] = a # generator of GL(1,F_p)
                A = Matrix(A, p)
                gen.append(A)

        order = order_gl(n, p)
        return GL(gen, order, p=p, **kw)

    # See:
    # Pairs of Generators for Matrix _Algebraics. I
    # D. E. Taylor 2006
    # https://www.maths.usyd.edu.au/u/don/papers/genAC.pdf
    # Also:
    # http://doc.sagemath.org/html/en/reference/groups/sage/groups/matrix_gps/symplectic.html
    
    @classmethod
    def Sp_4_2(cls, F, **kw):
        A = numpy.array([[1,0,1,1],[1,0,0,1],[0,1,0,1],[1,1,1,1]], dtype=scalar)
        B = numpy.array([[0,0,1,0],[1,0,0,0],[0,0,0,1],[0,1,0,0]], dtype=scalar)
        gen = [Matrix(A, 2), Matrix(B, 2)]
        return Sp(gen, 720, p=2, invariant_form=F, **kw)

    @classmethod
    def Symplectic(cls, nn, p=DEFAULT_P, **kw):
        gen = []
        assert nn%2 == 0
        n = nn//2
        assert isprime(p)
        F = numpy.zeros((nn, nn), dtype=scalar)
        for i in range(n):
            F[i, n+i] = 1
            F[n+i, i] = p-1
        F = Matrix(F, p)
        I = numpy.identity(nn, dtype=scalar)
        I = Matrix(I, p)
        G = Symplectic([I], order_sp(nn, p), p=p, invariant_form=F, **kw)
        return G

    @classmethod
    def Sp(cls, n, p=DEFAULT_P, **kw):
        gen = []
        assert n%2 == 0

        assert isprime(p)
        F = numpy.zeros((n, n), dtype=scalar)
        for i in range(n//2):
            F[i, n-i-1] = 1
            F[i + n//2, n//2-i-1] = p-1
        F = Matrix(F, p)

        for i in range(1, p):
            vals = set((i**k)%p for k in range(p+1))
            if len(vals)==p-1:
                fgen = i # generates GL(1, p)
                break
        else:
            assert 0
        for i in range(1, p):
            if (i*fgen)%p == 1:
                ifgen = i
                break
        else:
            assert 0

        if n==2:
            G = Sp.SL(2, p, invariant_form=F, **kw)
            return G
        if n==4 and p==2:
            return Sp.Sp_4_2(F, **kw)

        if p==2:
            m = n//2
            A = numpy.zeros((n, n), dtype=scalar)
            B = numpy.zeros((n, n), dtype=scalar)
            for i in range(n):
                A[i, i] = 1
            A[0, m-1] = 1
            A[0, n-1] = 1
            A[m, n-1] = 1
            for i in range(m-1):
                B[i+1, i] = 1
                B[i+m, i+m+1] = 1
            B[0, m] = 1
            B[n-1, m-1] = 1
            gen = [Matrix(A, 2), Matrix(B, 2)]
        else:
            m = n//2
            A = numpy.zeros((n, n), dtype=scalar)
            B = numpy.zeros((n, n), dtype=scalar)
            for i in range(n):
                A[i, i] = 1
            A[0, 0] = fgen
            A[n-1, n-1] = ifgen

            for i in range(m-1):
                B[i+1, i] = 1
                B[m+i, m+i+1] = 1
            B[0, 0] = 1
            B[0, m] = 1
            B[n-2, m-1] = 1
            B[n-1, m-1] = ifgen

            gen = [Matrix(A, p), Matrix(B, p)]

        G = Sp(gen, order_sp(n, p), p=p, invariant_form=F, **kw)
        return G

    # rip this stuff out of sage:

    @classmethod
    def SO_9_2(cls, **kw):
        p = 2
        gens = [
            [[1,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[1,1,0,0,1,0,0,0,1],[0,0,0,0,0,1,0,0,1],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1]],
            [[1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],[0,1,0,0,0,0,0,0,0]]]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, order=47377612800, p=p,
            invariant_bilinear_form = Matrix(
                [[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],[0,1,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0]]),
            invariant_quadratic_form = Matrix(
                [[1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0]]))
        return G

    @classmethod
    def SO_7_3(cls, **kw):
        p = 3
        gens = [
[[2,0,0,0,0,0,0],[0,2,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1]],
[[2,1,2,0,0,0,0],[0,2,0,0,0,0,0],[0,1,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1]],
[[1,2,2,0,0,0,0],[2,2,1,2,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1],[0,2,1,0,0,0,0],[2,1,1,1,0,0,0]],
        ]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, order=9170703360, p=p,
            invariant_bilinear_form = Matrix(
[[0,1,0,0,0,0,0],[1,0,0,0,0,0,0],[0,0,2,0,0,0,0],[0,0,0,2,0,0,0],[0,0,0,0,2,0,0],[0,0,0,0,0,2,0],[0,0,0,0,0,0,2]], p),
            invariant_quadratic_form = Matrix(
[[0,1,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1]], p))
        return G

    @classmethod
    def SO_7_2(cls, **kw):
        p = 2
        gens = [
            [[1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,1,0,0,0,0],[1,1,0,1,0,0,1],
             [0,0,0,0,1,0,1],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1]],
            [[1,0,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,1,0,0],
             [0,0,0,0,0,1,0],[0,0,0,0,0,0,1],[0,1,0,0,0,0,0]]]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, order=1451520, p=p,
            invariant_bilinear_form = Matrix(
                [[0,0,0,0,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1],
                 [0,1,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0]], p),
            invariant_quadratic_form = Matrix(
                [[1,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],
                 [0,1,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0]], p))
        return G

    @classmethod
    def SO_5_2(cls, **kw):
        p = 2
        gens = [
            [[1,0,0,0,0],[1,0,1,0,1],[1,0,1,1,1],[0,1,0,0,1],[0,1,1,1,1]],
            [[1,0,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1],[0,1,0,0,0]]]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, p=p,
            invariant_bilinear_form = Matrix(
                [[0,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1],[0,1,0,0,0],[0,0,1,0,0]], p),
            invariant_quadratic_form = Matrix(
                [[1,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0]], p))
        return G

    @classmethod
    def SO_5_3(cls, **kw):
        p = 3
        gens = [
                [[2,0,0,0,0],[0,2,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]],
                [[2,1,2,0,0],[0,2,0,0,0],[0,1,1,0,0],[0,0,0,1,0],[0,0,0,0,1]],
                [[1,2,2,0,0],[2,2,1,2,0],[0,0,0,0,1],[0,2,1,0,0],[2,1,1,1,0]],
        ]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, p=p,
            invariant_bilinear_form = Matrix(
                [[0,1,0,0,0],[1,0,0,0,0],[0,0,2,0,0],[0,0,0,2,0],[0,0,0,0,2]], p),
            invariant_quadratic_form = Matrix(
                [[0,1,0,0,0],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], p))
        return G



    @classmethod
    def SO_3_5(cls, **kw):
        p = 5
        gens = [
            [[2,0,0],[0,3,0],[0,0,1]], [[3,2,3],[0,2,0],[0,3,1]], [[1,4,4],[4,0,0],[2,0,4]]]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, p=p,
            invariant_bilinear_form = Matrix([[0,1,0],[1,0,0],[0,0,2]], p),
            invariant_quadratic_form = Matrix([[0,1,0],[0,0,0],[0,0,1]], p))
        return G

    @classmethod
    def SO_3_3(cls, **kw):
        p = 3
        gens = [
            [[2,0,0],[0,2,0],[0,0,1]],
            [[2,1,2],[0,2,0],[0,1,1]],
            [[1,2,2],[2,0,0],[2,0,2]],
            ]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, p=p,
            invariant_bilinear_form = Matrix([[0,1,0],[1,0,0],[0,0,2]], 3),
            invariant_quadratic_form = Matrix([[0,1,0],[0,0,0],[0,0,1]], 3))
        return G

    @classmethod
    def SO_3_2(cls, **kw):
        p = 2
        gens = [[[1,0,0],[1,1,1],[0,0,1]],[[1,0,0],[0,0,1],[0,1,0]]]
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, p=p,
            invariant_bilinear_form = Matrix([[0,0,0],[0,0,1],[0,1,0]], p),
            invariant_quadratic_form = Matrix([[1,0,0],[0,0,0],[0,1,0]], p))
        return G

    @classmethod
    def make(cls, gens, invariant_bilinear_form, invariant_quadratic_form, order=None, p=DEFAULT_P, **kw):
        gens = [Matrix(A, p) for A in gens]
        G = Algebraic(gens, order, p=p,
            invariant_bilinear_form = Matrix(invariant_bilinear_form, p),
            invariant_quadratic_form = Matrix(invariant_quadratic_form, p))
        return G


    @classmethod
    def SO_2_2(cls, **kw):
        return cls.make(
            [[[1,0],[0,1]],
             [[0,1],[1,0]]],
            [[0,1],[1,0]],
            [[0,1],[0,0]],
            2 , 2 )

    @classmethod
    def SO_4_2(cls, **kw):
        return cls.make(
            [[[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],
             [[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]]],
            [[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]],
            [[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]],
            72 , 2 )

    @classmethod
    def SO_6_2(cls, **kw):
        return cls.make(
            [[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]],
             [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,1],[0,0,0,0,0,1],[0,0,1,0,1,0]],
             [[1,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[1,0,0,0,0,1]],
             [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]],
            [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]],
            [[0,1,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,0],[0,0,0,0,0,1],[0,0,0,0,0,0]],
            40320 , 2 )

    @classmethod
    def SO_8_2(cls, **kw):
        return cls.make(
            [[[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0]],
             [[0,1,0,0,0,0,0,0],[1,0,0,0,1,0,0,0],[0,0,1,0,0,0,0,0],[0,1,0,1,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,1,0,0,0,1,0,0],[0,0,1,0,1,0,0,0]],
             [[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0],[0,0,1,0,1,0,0,0]]],
            [[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0]],
            [[0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0]],
            348364800 , 2 )

    @classmethod
    def SO_10_2(cls, **kw):
        return cls.make(
            [[[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,1,0]],
             [[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0,0,0],[0,0,1,0,1,0,0,0,0,0]],
             [[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,1,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[1,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[0,1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]]],
            [[0,1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,1,0]],
            [[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,0]],
            46998591897600 , 2 )

    @classmethod
    def SO_2_3(cls, **kw):
        return cls.make(
            [[[2,0],[0,2]],
             [[2,0],[0,2]],
             [[1,0],[0,1]]],
            [[0,1],[1,0]],
            [[0,1],[0,0]],
            2 , 3 )

    @classmethod
    def SO_4_3(cls, **kw):
        return cls.make(
            [[[0,2,2,2],[0,1,1,2],[1,0,2,0],[1,2,2,0]],
             [[0,2,2,1],[0,2,1,1],[1,1,0,2],[2,0,0,1]],
             [[1,0,0,0],[2,1,0,2],[0,0,1,0],[2,0,0,1]]],
            [[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,2]],
            [[0,1,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]],
            576 , 3 )

    @classmethod
    def SO_6_3(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0],[0,2,0,0,0,0],[0,0,0,0,1,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,1]],
             [[2,1,0,2,0,0],[0,2,0,0,0,0],[0,0,0,0,0,1],[0,1,0,1,0,0],[0,0,1,0,0,0],[0,0,0,0,1,0]],
             [[1,2,0,2,0,0],[2,2,0,1,2,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,2,0,1,0,0],[2,1,0,1,1,0]]],
            [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,2,0,0,0],[0,0,0,2,0,0],[0,0,0,0,2,0],[0,0,0,0,0,2]],
            [[0,1,0,0,0,0],[0,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]],
            12130560 , 3 )

    @classmethod
    def SO_8_3(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0,0,0],[0,2,0,0,0,0,0,0],[0,0,0,1,1,0,0,0],[0,0,2,2,1,0,0,0],[0,0,1,2,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
             [[2,1,0,2,0,0,0,0],[0,2,0,0,0,0,0,0],[0,0,0,0,1,1,0,0],[0,1,0,1,0,0,0,0],[0,0,2,0,2,1,0,0],[0,0,1,0,2,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
             [[1,2,0,2,0,0,0,0],[2,2,0,1,2,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,2,0,1,0,0,0,0],[2,1,0,1,1,0,0,0]]],
            [[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,2,0,0,0,0],[0,0,0,0,2,0,0,0],[0,0,0,0,0,2,0,0],[0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,2]],
            [[0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,2,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
            19808719257600 , 3 )

    @classmethod
    def SO_10_3(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0,0,0,0,0],[0,2,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[2,1,0,2,0,0,0,0,0,0],[0,2,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,1,0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[1,2,0,2,0,0,0,0,0,0],[2,2,0,1,2,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[0,2,0,1,0,0,0,0,0,0],[2,1,0,1,1,0,0,0,0,0]]],
            [[0,1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,2,0,0,0,0,0,0,0],[0,0,0,2,0,0,0,0,0,0],[0,0,0,0,2,0,0,0,0,0],[0,0,0,0,0,2,0,0,0,0],[0,0,0,0,0,0,2,0,0,0],[0,0,0,0,0,0,0,2,0,0],[0,0,0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,0,0,2]],
            [[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
            None , 3 )

    # these have e == -1:

    @classmethod
    def SO_2_2_1(cls, **kw):
        return cls.make(
            [[[1,1],[1,0]],
             [[1,0],[1,1]]],
            [[0,1],[1,0]],
            [[1,1],[0,1]],
            6 , 2 )

    @classmethod
    def SO_4_2_1(cls, **kw):
        return cls.make(
            [[[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],
             [[0,1,0,0],[1,1,1,0],[0,1,0,1],[0,0,1,0]]],
            [[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]],
            [[0,1,0,0],[0,0,0,0],[0,0,1,1],[0,0,0,1]],
            120 , 2 )

    @classmethod
    def SO_6_2_1(cls, **kw):
        return cls.make(
            [[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]],
             [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,1],[0,0,0,0,0,1],[0,0,1,0,1,1]],
             [[1,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[1,0,0,0,0,1]],
             [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]],
            [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]],
            [[0,1,0,0,0,0],[0,0,0,0,0,0],[0,0,1,1,0,0],[0,0,0,1,0,0],[0,0,0,0,0,1],[0,0,0,0,0,0]],
            51840 , 2 )

    @classmethod
    def SO_8_2_1(cls, **kw):
        return cls.make(
            [[[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0]],
             [[0,1,0,0,0,0,0,0],[1,0,0,0,1,0,0,0],[0,0,1,0,0,0,0,0],[0,1,0,1,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,1,0,0,0,1,0,0],[0,1,1,0,1,1,0,0]],
             [[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0],[0,0,1,0,1,1,0,0]]],
            [[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0]],
            [[0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,1,1,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0]],
            394813440 , 2 )

    @classmethod
    def SO_10_2_1(cls, **kw):
        return cls.make(
            [[[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,1,0]],
             [[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0,0,0],[0,0,1,0,1,1,0,0,0,0]],
             [[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,1,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[1,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[0,1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]]],
            [[0,1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,1,0]],
            [[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,1,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,0]],
            50030759116800 , 2 )

    @classmethod
    def SO_12_2_1(cls, **kw):
        # sage code:
        # G = SO(12, GF(2), -1)
        # for g in G.gens():
        #     print(g)
        #     print()

        data = ("""
        [1 0 0 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 1]
        [0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0]""",
        """
        [0 1 0 0 0 0 0 0 0 0 0 0]
        [1 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0]
        [0 1 0 1 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 1]
        [0 1 0 0 0 1 0 0 0 0 0 0]
        [0 1 1 0 1 1 0 0 0 0 0 0]""",
        """
        [0 1 0 0 0 0 0 0 0 0 0 0]
        [1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 1]
        [0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 1 0 1 1 0 0 0 0 0 0]""")
        invariant_bilinear_form = parse("""
        [0 1 0 0 0 0 0 0 0 0 0 0]
        [1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 1]
        [0 0 0 0 0 0 0 0 0 0 1 0]
        """)
        invariant_quadratic_form = parse("""
        [0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 1 1 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 1]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        """)
        gen = [parse(d) for d in data]
        #for g in gen:
        #    print(g)
        return cls.make(gen, invariant_bilinear_form, invariant_quadratic_form, 
            103231467131240448000, 2)

    @classmethod
    def SO_2_3_1(cls, **kw):
        return cls.make(
            [[[2,2],[2,1]],
             [[1,1],[1,2]],
             [[1,0],[0,1]]],
            [[2,1],[1,1]],
            [[1,1],[0,2]],
            4 , 3 )

    @classmethod
    def SO_4_3_1(cls, **kw):
        return cls.make(
            [[[0,2,0,0],[2,1,0,1],[0,2,0,1],[0,0,1,0]],
             [[1,2,0,2],[1,1,1,1],[1,0,0,1],[1,2,1,2]],
             [[1,0,0,0],[1,1,2,1],[2,0,1,0],[1,0,0,1]]],
            [[0,1,0,0],[1,0,0,0],[0,0,2,0],[0,0,0,2]],
            [[0,1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]],
            720 , 3 )

    @classmethod
    def SO_6_3_1(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0],[0,2,0,0,0,0],[0,0,0,1,1,0],[0,0,2,2,1,0],[0,0,1,2,1,0],[0,0,0,0,0,1]],
             [[2,1,0,2,0,0],[0,2,0,0,0,0],[0,0,0,0,1,1],[0,1,0,1,0,0],[0,0,2,0,2,1],[0,0,1,0,2,1]],
             [[1,2,0,2,0,0],[2,2,0,1,2,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,2,0,1,0,0],[2,1,0,1,1,0]]],
            [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,2,0,0],[0,0,0,0,2,0],[0,0,0,0,0,2]],
            [[0,1,0,0,0,0],[0,0,0,0,0,0],[0,0,2,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]],
            13063680 , 3 )

    @classmethod
    def SO_8_3_1(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0,0,0],[0,2,0,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
             [[2,1,0,2,0,0,0,0],[0,2,0,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,1,0,1,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
             [[1,2,0,2,0,0,0,0],[2,2,0,1,2,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,2,0,1,0,0,0,0],[2,1,0,1,1,0,0,0]]],
            [[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,2,0,0,0,0,0],[0,0,0,2,0,0,0,0],[0,0,0,0,2,0,0,0],[0,0,0,0,0,2,0,0],[0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,2]],
            [[0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
            20303937239040 , 3 )

    @classmethod
    def SO_10_3_1(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0,0,0,0,0],[0,2,0,0,0,0,0,0,0,0],[0,0,0,1,1,0,0,0,0,0],[0,0,2,2,1,0,0,0,0,0],[0,0,1,2,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[2,1,0,2,0,0,0,0,0,0],[0,2,0,0,0,0,0,0,0,0],[0,0,0,0,1,1,0,0,0,0],[0,1,0,1,0,0,0,0,0,0],[0,0,2,0,2,1,0,0,0,0],[0,0,1,0,2,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[1,2,0,2,0,0,0,0,0,0],[2,2,0,1,2,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[0,2,0,1,0,0,0,0,0,0],[2,1,0,1,1,0,0,0,0,0]]],
            [[0,1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,2,0,0,0,0,0,0],[0,0,0,0,2,0,0,0,0,0],[0,0,0,0,0,2,0,0,0,0],[0,0,0,0,0,0,2,0,0,0],[0,0,0,0,0,0,0,2,0,0],[0,0,0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,0,0,2]],
            [[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,2,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
            None , 3 )

    @classmethod
    def SO_2_5(cls, **kw):
        return cls.make(
            [[[2,0],[0,3]],
             [[3,0],[0,2]],
             [[1,0],[0,1]]],
            [[0,1],[1,0]],
            [[0,1],[0,0]],
            4 , 5 )

    @classmethod
    def SO_4_5(cls, **kw):
        return cls.make(
            [[[0,1,0,0],[1,4,4,0],[0,0,0,4],[0,3,4,0]],
             [[0,4,0,0],[4,4,0,3],[0,4,0,4],[0,0,4,0]],
             [[4,0,0,0],[0,4,0,0],[0,0,1,0],[0,0,0,1]]],
            [[0,1,0,0],[1,0,0,0],[0,0,2,0],[0,0,0,2]],
            [[0,1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]],
            14400 , 5 )

    @classmethod
    def SO_6_5(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0],[0,3,0,0,0,0],[0,0,0,0,1,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,1]],
             [[3,2,0,3,0,0],[0,2,0,0,0,0],[0,0,0,0,0,1],[0,3,0,1,0,0],[0,0,1,0,0,0],[0,0,0,0,1,0]],
             [[1,4,0,4,0,0],[4,2,0,1,4,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,2,0,1,0,0],[2,3,0,3,1,0]]],
            [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,2,0,0,0],[0,0,0,2,0,0],[0,0,0,0,2,0],[0,0,0,0,0,2]],
            [[0,1,0,0,0,0],[0,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]],
            29016000000 , 5 )

    @classmethod
    def SO_8_5(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0,0,0],[0,3,0,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
             [[3,2,0,3,0,0,0,0],[0,2,0,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,3,0,1,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
             [[1,4,0,4,0,0,0,0],[4,2,0,1,4,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,2,0,1,0,0,0,0],[2,3,0,3,1,0,0,0]]],
            [[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,2,0,0,0,0,0],[0,0,0,2,0,0,0,0],[0,0,0,0,2,0,0,0],[0,0,0,0,0,2,0,0],[0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,2]],
            [[0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],
            None , 5 )

    @classmethod
    def SO_10_5(cls, **kw):
        return cls.make(
            [[[2,0,0,0,0,0,0,0,0,0],[0,3,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[3,2,0,3,0,0,0,0,0,0],[0,2,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,3,0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
             [[1,4,0,4,0,0,0,0,0,0],[4,2,0,1,4,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[0,2,0,1,0,0,0,0,0,0],[2,3,0,3,1,0,0,0,0,0]]],
            [[0,1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,2,0,0,0,0,0,0,0],[0,0,0,2,0,0,0,0,0,0],[0,0,0,0,2,0,0,0,0,0],[0,0,0,0,0,2,0,0,0,0],[0,0,0,0,0,0,2,0,0,0],[0,0,0,0,0,0,0,2,0,0],[0,0,0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,0,0,2]],
            [[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1]],
            None , 5 )

    @classmethod
    def SO(cls, n, p=DEFAULT_P, e=None, **kw):
        if e==None or e==1:
            attr = "SO_%d_%d"%(n, p)
        else:
            attr = "SO_%d_%d_1"%(n, p)
        method = getattr(cls, attr)
        if method:
            return method(**kw)
        assert 0, (n, p)


def _basep(n, p):
    e = n//p
    q = n%p
    if n == 0:
        return '0'
    elif e == 0:
        return str(q)
    else:
        return _basep(e, p) + str(q)


#def test_so():
#
#    n = argv.get("n", 3)
#    m = argv.get("m", 1)
#    p = argv.get("p", 2)
#    e = argv.get("e", 1)
#
#    G = Algebraic.SO(n, p, e)
#    if G.order is None:
#        G.get_elements()
#    print("|G| =", G.order)
#
#    B = G.invariant_bilinear_form
#    Q = G.invariant_quadratic_form
#
#    print("B =")
#    print(B)
#    print("Q =")
#    print(Q)
#
#    for g in G.gen:
#        assert g * B * g.transpose() == B
#        #assert g * Q * g.transpose() == Q # no ..
#
#    while 1:
#        M = [randint(0, p-1) for i in range(n*m)]
#        M = numpy.array(M)
#        M.shape = (m, n)
#        M = normal_form_p(M, p, truncate=True)
#        if len(M) < m:
#            continue
#        M = Matrix(M, p, shape=(m,n))
#        if argv.null:
#            w = M * M.transpose() # only works when p>2
#        elif p>2:
#            w = M * B * M.transpose() # only works when p>2
#        else:
#            w = M * Q * M.transpose() # use this one for p==2
#        if w.is_zero():
#            break
#    print("found: M =")
#    print(M)
#
#    orbit = set([M.key[1]])
#    bdy = set([M])
#    while bdy:
#        print("(%s)"%len(bdy), end="", flush=True)
#        _bdy = set()
#        for M in bdy:
#          for g in G.gen:
#            gM = M*g
#            gM = gM.normal_form()
#            s = gM.key[1]
#            if s in orbit:
#                continue
#            _bdy.add(gM)
#            orbit.add(s)
#        bdy = _bdy
#    print()
#    N = len(orbit)
#    print(N, "= [%s]_%d"%(_basep(N, p), p))


#def test_SO_generators():
#    # find generators
#    print("test_SO_generators()")
#    n = argv.get("n", 5)
#    p = argv.get("p", 2)
#    assert p==2
#
#    def A003053(n): 
#        return (1 << (n//2)**2)*prod((1 << i)-1 for i in range(2, 2*((n-1)//2)+1, 2))
#    order = A003053(n)
#    print("|O(%d,2)| = %d"%(n, order))
#
##    SO = Algebraic.SO(n=n, p=p)
##    if SO.order is None:
##        SO.get_elements()
##    N = len(SO)
##    print("|SO| =", N)
##    if n%2==0:
##        SOm = Algebraic.SO(n=n, p=p, e=-1)
##        Nm = len(SO)
##        print("|SO-| =", Nm)
#
#    gen = []
#    I = Matrix.identity(n, p)
#    while 1:
#        M = numpy.random.randint(0, p, (n, n))
#        M = Matrix(M, p, shape=(n,n))
#        if M * M.transpose() == I:
#            gen.append(M)
#            print("gen:", len(gen))
#            if len(gen)>1:
#                G = mulclose(gen, verbose=True)
#                N0 = len(G)
#                print("|G| =", N0)
#                if N0 == order:
#                    break
#    for g in gen:
#        print(g)

#    count = 0
#    for g in SO:
#        if g in G:
#            count += 1
#    print("intersection:", count)
    
#
#def test_SO():
#    n = argv.get("n", 6)
#    m = argv.get("m", 3)
#    p = 2
#
#    if n==6:
#        gen = [
#            [[1,1,1,0,1,1],
#             [0,1,0,1,1,0],
#             [1,0,1,1,1,1],
#             [0,1,0,1,0,1],
#             [0,1,1,1,0,0],
#             [1,1,0,1,0,0]],
#            [[1,1,1,1,1,0],
#             [0,1,1,0,0,1],
#             [0,0,1,1,0,1],
#             [0,0,1,0,1,1],
#             [1,1,0,1,1,1],
#             [1,0,1,0,0,1]]
#        ]
#        order = 23040
#    elif n==7:
#        gen = [
#            [[1,1,0,1,1,1,0],
#             [0,1,0,0,1,0,1],
#             [0,1,1,1,1,1,0],
#             [0,1,0,1,0,0,1],
#             [1,0,1,0,0,1,0],
#             [1,1,1,1,1,0,0],
#             [0,0,0,1,1,0,1]],
#            [[0,1,1,0,1,0,0],
#             [1,1,1,1,0,0,1],
#             [1,0,1,1,1,0,1],
#             [0,0,0,0,0,1,0],
#             [1,1,0,0,1,0,0],
#             [0,1,0,1,1,0,0],
#             [0,1,0,0,1,0,1]]
#        ]
#        order = 1451520
#    else:
#        assert 0
#
#    gen = [Matrix(g) for g in gen]
#    if argv.gen:
#        G = mulclose(gen, verbose=True)
#        assert len(G) == order
#
#    for count in range(1):
#        while 1:
#            M = numpy.random.randint(0, p, (m, n))
#            M = normal_form_p(M, p, truncate=True)
#            if len(M) < m:
#                continue
#            M = Matrix(M, p)
#            w = M * M.transpose()
#            if w.is_zero():
#                break
#        print(M)
#    
#        orbit = set([M.key[1]])
#        bdy = set([M])
#        while bdy:
#            #print("(%s)"%len(bdy), end="", flush=True)
#            _bdy = set()
#            for M in bdy:
#              for g in gen:
#                gM = M*g
#                #w = gM*gM.transpose()
#                #assert w.is_zero()
#                gM = gM.normal_form()
#                s = gM.key[1]
#                if s in orbit:
#                    continue
#                _bdy.add(gM)
#                orbit.add(s)
#                #print(gM)
#            bdy = _bdy
#        #print()
#        N = len(orbit)
#        print(N, "= [%s]_%d"%(_basep(N, p), p))
#

#def is_orthogonal(M):
#    if isinstance(M, numpy.ndarray):
#        M = Matrix(M)
#    assert isinstance(M, Matrix)
#    n = len(M)
#    return M*M.transpose() == Matrix.identity(n,2)
#
#def offdiag(n):
#    A = numpy.zeros((n,n),dtype=int)
#    A[:] = 1
#    A += numpy.identity(n,dtype=int)
#    A %= 2
#    return A
#
#def antidiag(n):
#    A = numpy.zeros((n,n),dtype=int)
#    for i in range(n):
#        A[i,n-i-1] = 1
#    return A
#
#
#def get_SO_gens(n):
#    assert n%2 == 1
#    assert is_orthogonal(offdiag(n-1))
#    assert is_orthogonal(antidiag(n))
#
#    gen = []
#    for i in range(n-1):
#        items = list(range(n))
#        items[i:i+2] = i+1, i
#        M = Matrix.perm(items)
#        gen.append(M)
#    
#    for k in range(1, n//2+1, 2):
#        A = antidiag(n)
#        A[k:,:-k] = offdiag(n-k)
#        M = Matrix(A)
#        #print("gen:")
#        #print(M)
#        assert is_orthogonal(M)
#        gen.append(M)
#
#    return gen
#
#
#def find_orbit(gen, M, verbose=False):
#    orbit = set([M.key[1]])
#    bdy = set([M])
#    yield M
#    while bdy:
#        if verbose:
#            print("(%s)"%len(bdy), end="", flush=True)
#        _bdy = set()
#        for M in bdy:
#          for g in gen:
#            gM = M*g
#            #w = gM*gM.transpose()
#            #assert w.is_zero()
#            gM = gM.normal_form()
#            s = gM.key[1]
#            if s in orbit:
#                continue
#            _bdy.add(gM)
#            orbit.add(s)
#            #print(gM)
#            yield gM
#        bdy = _bdy
#    if verbose:
#        print()
#    #return orbit
#
#
#def build_SO():
#    n = argv.get("n", 5)
#    m = argv.get("m", n//2)
#    p = 2
#
#    gen = get_SO_gens(n)
#    u = numpy.array([1]*n)
#
#    best_d = 1
#
#    #for m in range(1, n//2+1):
#    #m = n//2
#    #while m:
#    for m in [m]:
#        A = zeros2(m, n)
#        for i in range(m):
#            A[i,2*i] = 1
#            A[i,2*i+1] = 1
#        M = Matrix(A)
#        print("M =")
#        print(M)
#    
#        count = 0
#        for M in find_orbit(gen, M, verbose=True):
#            count += 1
#            if M.A.sum(0).min() == 0:
#                continue
#            d = n
#            for v in span(M):
#                d1 = ((v+u)%2).sum()
#                if d1<d:
#                    d = d1
#            if d > best_d:
#                best_d = d
#                print()
#                print(M, d)
#                print()
#        print(count, "= [%s]_%d"%(_basep(count, p), p))
#
#    if argv.mulclose:
#        G = mulclose(gen, verbose=True)
#        print(len(G))



class GL(Algebraic):
    def all_figs(self, dims, p=DEFAULT_P):
        dims = [self.n] + dims
        return list(Figure.qchoose(dims, p))

    def get_weyl(self):
        n = self.n
        items = list(range(n))
        gen = []
        for ii in range(n-1):
            perm = list(range(n))
            perm[ii:ii+2] = perm[ii+1], perm[ii]
            M = Matrix.perm(perm, self.p, name="w%d"%ii)
            gen.append(M)
            #print(M.name)
            #print(M.shortstr())
        return Algebraic(gen, p=self.p)


class Sp(Algebraic):

    def qchoose(self, m):
        F = self.invariant_form
        p = self.p
        n = self.n
        for M in qchoose(n, m, p):
            M = Matrix(M, p)
            A = M*F*M.transpose()
            if A.is_zero():
                yield M

    def all_figs(self, dims, p=DEFAULT_P):
        m = dims[0]
        items = [[M.A for M in self.qchoose(m)]]
        for m1 in dims[1:]:
            assert m1<=m
            items.append(list(qchoose(m, m1)))
            m = m1
        n = len(items)
        #print("Sp.all_figs")
        #print('\t', items)
        figs = []
        for select in cross(items):
            A = select[0]
            flag = [A]
            for i in range(n-1):
                A = numpy.dot(select[i+1], A) % p
                flag.append(A)
            flag = list(reversed(flag))
            flag = Figure(flag)
            figs.append(flag)
        return figs

    def slow_grassmanian(self, m):
        n = self.n//2
        p = self.p
        F = self.invariant_form
    
        bits = 2*n*m
        found = {}
    
        A = numpy.zeros((m, 2*n), dtype=int)
        pairs = self.get_pairs()
        cols = None
        #cols = reduce(add, pairs) # fail
        #cols = [pair[0] for pair in pairs] + [pair[1] for pair in pairs] # fail
        #print("slow_grassmanian:", cols)
        jdxs = list(numpy.ndindex(A.shape))
        assert len(jdxs) == bits
        for idxs in cross( [tuple(range(p))]*bits ):
            for (bit, jdx) in zip(idxs, jdxs):
                A[jdx] = bit
            M = Matrix(A)
            B = M*F*M.transpose()
            if not B.is_zero():
                continue
            M = M.normal_form(cols)
            if rank(M.A) < m:
                continue
            if M not in found:
                yield M
            found[M] = M

    def grassmanian(self, m):
        p = self.p
        if p != 2:
            #return self.qchoose(m)
            for M in self.qchoose(m):
                yield M
            return
        n = self.n//2
        F = self.invariant_form
        found = set()
        for left in all_codes(m, n):
          for right in all_matrices(m, n):
            M = numpy.concatenate((left, right), axis=1)
            assert M.shape == (m, 2*n), M.shape
            M = normal_form(M, p)
            if rank(M) < m:
                continue
            M = Matrix(M, p)
            A = M*F*M.transpose()
            if A.is_zero() and M not in found:
                found.add(M)
                yield M

    def get_pairs(self):
        n = self.n//2
        pairs = [(i, 2*n-i-1) for i in range(n)]
        return pairs

    def get_blocks(self, M):
        pairs = self.get_pairs()
        A = M.A[:, [pair[0] for pair in pairs]]
        B = M.A[:, [pair[1] for pair in pairs]]
        #A = Matrix(A, self.p)
        #B = Matrix(B, self.p)
        return A, B

    def from_blocks(self, A, B):
        n = self.n//2
        m = len(A)
        assert A.shape == (m, n)
        assert B.shape == (m, n)
        M = numpy.zeros((m, 2*n))
        pairs = self.get_pairs()
        M[:, [pair[0] for pair in pairs]] = A
        M[:, [pair[1] for pair in pairs]] = B
        #M = numpy.concatenate((A, B), axis=1)
        #assert M.shape == (m, 2*n)
        #print(M)
        M = Matrix(M, self.p)
        return M

    def is_isotropic(self, M):
        F = self.invariant_form
        MM = M*F*M.transpose()
        return MM.is_zero()

    def get_weyl(self):
        nn = self.n
        pairs = self.get_pairs()
        k = len(pairs)
        gen = []
        for ii in range(k-1):
            idxs = list(range(k))
            idxs[ii:ii+2] = idxs[ii+1], idxs[ii]
            qairs = [pairs[idx] for idx in idxs]
            src = reduce(add, pairs)
            tgt = reduce(add, qairs)
            #print(src, "-->", tgt)
            perm = list(range(nn))
            for i, j in zip(src, tgt):
                perm[i] = j
            #print('\t', perm)
            M = Matrix.perm(perm, self.p, name="w%d"%ii)
            gen.append(M)
            #print(M.name)
            #print(M.shortstr())
        perm = list(range(nn))
        a, b = pairs[-1]
        perm[a], perm[b] = perm[b], perm[a]
        M = Matrix.perm(perm, self.p, name="w%d"%(k-1))
        gen.append(M)
        #print(M.name)
        #print(M.shortstr())
        #return self.get_all_weyl()
        #print(len(gen))
        #print(len(mulclose(gen)))
        return Algebraic(gen, p=self.p)

    def get_all_weyl(self):
        nn = self.n
        pairs = self.get_pairs()
        k = len(pairs)
        I = Matrix.identity(nn, self.p)
        flips = []
        for (a, b) in pairs:
            perm = list(range(nn))
            perm[a], perm[b] = b, a
            M = Matrix.perm(perm)
            flips.append((I, M))
        W = []
        for qairs in allperms(pairs):
            perm = [None]*nn
            for src, tgt in zip(pairs, qairs):
                perm[src[0]] = tgt[0]
                perm[src[1]] = tgt[1]
            assert None not in perm
            M = Matrix.perm(perm, self.p)
            for items in cross(flips):
                N = reduce(mul, items)
                W.append(M*N)
        return W

    def get_borel(self):
        F = self.invariant_form
        pairs = self.get_pairs()
        gen = []
        lookup = {}
        for (i, j) in pairs:
            assert i<j
            A = numpy.identity(self.n, dtype=scalar)
            A[i, j] = 1
            M = Matrix(A, self.p, name="E%d%d"%(j,i))
            gen.append(M)
            lookup[j, i] = M
    
        k = len(pairs)
        for i in range(k):
          for j in range(i+1, k):
    
            a, b = pairs[i]
            c, d = pairs[j]
            assert d < b
            assert a < c
    
            A = numpy.identity(self.n, dtype=scalar)
            A[a, c] = 1
            A[d, b] = self.p-1
    
            M = Matrix(A, self.p, name="E%d%d"%(b,d))
            gen.append(M)
            lookup[b, d] = M # b>d
            assert (c, a) not in lookup
            lookup[c, a] = M

            # we don't need these gen's but we need the lookup
            c, d = d, c
            assert d < b
            assert a < c
            A = numpy.identity(self.n, dtype=scalar)
            A[a, c] = 1
            A[d, b] = 1
    
            M = Matrix(A, self.p, name="E%d%d"%(b,d))
            gen.append(M)
            lookup[b, d] = M
            assert (c, a) not in lookup
            lookup[c, a] = M
    
        for M in gen:
            assert M*F*M.transpose() == F
    
        B = Algebraic(gen, p=self.p)
        B.lookup = lookup
        return B


    U = None # cache
    def get_zip_uturn(self):
        if self.U is not None:
            return self.U
        nn = self.n
        n = nn//2
        U = numpy.zeros((nn,nn), dtype=scalar)
        cols = [2*i for i in range(n)] + list(reversed([2*i+1 for i in range(n)]))
        for i in range(nn):
            U[i, cols[i]] = 1
        U = Matrix(U, name="U")
        self.U = U
        return U

    def from_ziporder(self, M):
        nn = self.n
        assert M.shape == (nn, nn)
        U = self.get_zip_uturn()
        M1 = U*M*U.transpose()
        return M1

    def to_ziporder(self, M):
        nn = self.n
        assert M.shape == (nn, nn)
        U = self.get_zip_uturn()
        M1 = U.transpose()*M*U
        return M1




class Symplectic(Sp):
    def get_pairs(self):
        n = self.n//2
        pairs = [(i, i+n) for i in range(n)]
        return pairs


def test_grassmanian_fail():
    # here we are building bigger grassmanian's from smaller ones
    # _looking for the Sp(q)-deformed pascal triangle (and so far failing)

    p = argv.get("p", 2)
    n = argv.get("n", 3)
    m = argv.get("m", 2)

    print("test_grassmanian: n=%d, m=%d"%(n, m))

    G = Algebraic.Sp(2*n, p)
    G1 = Algebraic.Sp(2*(n+1), p)

    items = list(G.grassmanian(m))
    print(len(items))

    #gr1 = set(G1.grassmanian(m+1))
    #print(len(gr1))

    if 0:
        items1 = set()
        for M in items:
            assert G.is_isotropic(M)
            A, B = G.get_blocks(M)
            M1 = G.from_blocks(A, B)
            assert M1 == M
    
            # try adding a col
            A1 = numpy.zeros((m, n+1))
            B1 = numpy.zeros((m, n+1))
    
            bits = [tuple(range(p))]*m
            found = set()
            for left in cross(bits):
              A1[:, :n] = A
              A1[:, n] = left
              for right in cross(bits):
                B1[:, :n] = B
                B1[:, n] = right
        
                M1 = G1.from_blocks(A1, B1)
                if G1.is_isotropic(M1):
                    M1 = M1.normal_form()
                    found.add(M1)
            items1.update(found)
            print(len(found), end=" ", flush=True)
        print()
        print(len(items1))

    if 0:
        items1 = set()
        for M in items:
            assert G.is_isotropic(M)
            A, B = G.get_blocks(M)
            M1 = G.from_blocks(A, B)
            assert M1 == M
    
            # try adding a col & a row
            A1 = numpy.zeros((m+1, n+1))
            B1 = numpy.zeros((m+1, n+1))
    
            n1_bits = [tuple(range(p))]*(n+1)
            m_bits = [tuple(range(p))]*m
            found = set()
            for n1_left in cross(n1_bits):
             #for m_bit in cross(m_bits):
              A1[:m, :n] = A
              A1[m, :n+1] = n1_left
              #A1[:m, n] = m_bit # yes this finds all the bigger grassmanian's
              for n1_right in cross(n1_bits):
                B1[:m, :n] = B
                B1[m, :n+1] = n1_right
        
                M1 = G1.from_blocks(A1, B1)
                if not G1.is_isotropic(M1):
                    continue
                M1 = M1.normal_form()
                if rank(M1.A) == m+1:
                    #assert M1 in gr1
                    found.add(M1)
            items1.update(found)
            print(len(found), end=" ", flush=True)
            #break
        print()
        print(len(items1))
    
        
    F1 = G1.invariant_form
    ps = tuple(range(p))
    items1 = set()
    for M in items:
        assert G.is_isotropic(M)
        A, B = G.get_blocks(M)
        M1 = G.from_blocks(A, B)
        assert M1 == M
        print(M.shortstr())

        # try adding a col & a row
        A1 = numpy.zeros((m+1, n+1))
        B1 = numpy.zeros((m+1, n+1))

        found = set()
        A1[:m, :n] = A
        A1[m, n] = 1
        for bits in cross([ps]*(m+1)):
            B1[:m, :n] = B
            B1[:, n] = bits
    
            M1 = G1.from_blocks(A1, B1)
            if not G1.is_isotropic(M1):
                continue
            
            M1 = M1.normal_form()
            #assert M2==M1 # fail...
            assert( rank(M1.A) == m+1 )
            found.add(M1)
            print(M1.shortstr())
            print((M1 * F1 * M1.transpose()).shortstr())
            print()
        items1.update(found)
        #print(len(found), end=" ", flush=True)
        print("="*79)
        #break
    print()
    print(len(items1))


def show_cell(items):
    items = [M.A for M in items]
    M = numpy.array(items)

    smap = SMap()
    for (i, j) in numpy.ndindex(M.shape[1:]):
        if numpy.alltrue(M[:, i, j] == 0):
            smap[i,j] = "."
        elif numpy.alltrue(M[:, i, j] == 1):
            smap[i,j] = "1"
        else:
            smap[i,j] = "*"
    return str(smap)
    

def test_grassmanian():
    # here we are building bigger grassmanian's from smaller ones
    # _looking for the Sp(q)-deformed pascal triangle (and so far failing)

    p = argv.get("p", 2)
    n = argv.get("n", 3)
    m = argv.get("m", 2)

    print("test_grassmanian: n=%d, m=%d"%(n, m))

    G = Algebraic.Sp(2*n, p)
    #Ms = list(G.grassmanian(m))
    Ms = list(G.slow_grassmanian(m))

    # wrong symplectic form ... scrambles the bruhat cells
    #G = Algebraic.Symplectic(2*n, p)
    #Ms = list(G.grassmanian(m))

    print(len(Ms))

    bruhat = {}

    for M in Ms:
        #assert M == M.normal_form()
        key = M.get_pivots()
        items = bruhat.setdefault(key, set())
        items.add(M)

    print(len(bruhat))
    values = list(bruhat.values())
    values = [len(val) for val in values]
    values.sort()
    print(values)
    for key, values in bruhat.items():
        s = show_cell(values)
        print(key, len(values), "*", p**s.count("*")//len(values), "=", p**s.count("*"))
        print(s)
#        for M in values:
#            A, B = G.get_blocks(M)
#            print(shortstrx(A, B))
#            print()
#        print("="*79)


    
def test_1():
    n = argv.get("n", 3)
    m = argv.get("m", 2)
    p = argv.get("p", 2)

    G = Algebraic.Sp(2*n, p)
    F = G.invariant_form

    G1 = Algebraic.Sp(2*n + 2, p)
    F1 = G1.invariant_form

    items = set()
    for M in G.qchoose(m):
        A = M.A
        #left = A[:, :n]
        #right = A[:, n:]
        for u in all_matrices(m, 1):
          for v in all_matrices(m, 1):
            B = numpy.concatenate((u, A, v), axis=1)
            B = Matrix(B, p)
            #print(B.shape, F1.shape)
            C = B*F1*B.transpose()
            if not C.is_zero():
                continue
            B = B.normal_form()
            assert rank(B.A) == m
            items.add(B)
    print(len(items))
            

        
        
def get_subgroup(G, geom, check=False):
    if type(geom) is str:
        A = parse(geom)
    else:
        A = geom

    H = []
    for op in G:
        if op.mask(A) == op:
            H.append(op)

    if check:
        for a in H:
          for b in H:
            assert (a*b).mask(A) == (a*b)
    assert len(G) % len(H) == 0

    return H


def get_permrep(G):
    """
    permutation action of G on the non-zero vectors (in lexicographic order)
    """
    from bruhat import gset
    assert len(G)
    op = G[0]
    n, _ = op.shape
    p = op.p
    space = list(numpy.array(v) for v in enum2(n) if sum(v)!=0)
    lookup = dict((v.tobytes(), idx) for (idx, v) in enumerate(space))
    perms = []
    rep = {}
    for op in G:
        idxs = [lookup[(numpy.dot(op.A, v)%p).tobytes()] for v in space]
        perm = gset.Perm(idxs)
        rep[op] = perm
        perms.append(perm)
    G = gset.Group(perms)
    G.get_gens()
    G.rep = rep
    G.space = space
    return G


def test_dynkin():

    G = Algebraic.GL(4, 2)
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

    PonA = "1111 .111 .111 ..11"
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


class Building(object):
    def __init__(self, G):
        self.G = G
        self.nn = G.n
        self.W = G.get_weyl()
        self.B = G.get_borel()

    def decompose(self, g):
        nn = self.nn
        n = nn//2
        lookup = self.B.lookup
        b1 = b2 = self.G.I
        src = nn-1
        pivots = [] # cols
        while src >= 0:
            for col in range(nn):
                if g[src, col]:
                    break
            else:
                assert 0
            pivots.append((src, col))
            tgt = src-1
            while tgt >= 0:
                if g[tgt, col]:
                    g = lookup[src, tgt]*g
                    b1 = lookup[src, tgt]*b1
                tgt -= 1
            src -= 1

        for row, src in pivots:
            for tgt in range(src+1, nn):
                if g[row, tgt]==0:
                    continue
                b = lookup.get((tgt, src))
                if b is not None:
                    g = g*b
                    b2 = b2*b
        return b1, b2


def test_weyl():
    from bruhat.action import Perm
    from bruhat.action import Group

    n = argv.get("n", 3)
    letters = "abcdef"[:n]
    items = ["+"+c for c in letters] + ["-"+c for c in letters]
    pairs = [("+"+c, "-"+c) for c in letters]

    print(items)

    gen = []
    for i in range(n-1):
        perm = {c:c for c in items}
        a, b = pairs[i]
        c, d = pairs[i+1]
        perm[a] = c
        perm[c] = a
        perm[b] = d
        perm[d] = b
        gen.append(Perm(perm, items))
    perm = {c:c for c in items}
    a, b = pairs[-1]
    perm[a] = b
    perm[b] = a
    gen.append(Perm(perm, items))
    G = Group.generate(gen)
    print("|G| =", len(G))

    #for H in G.subgroups():
    #    print("\t", len(H))

    H = []
    a, b, c, aa, bb, cc = items
    #fix = [bb, aa]
    fix = {a, b}
    #print("fix:", fix)
    for g in G:
        send = {g[x] for x in fix}
        if fix == send:
            H.append(g)
            #print("send:", send)
    H = Group(H, items)
    print("|H| =", len(H))
    X = G.left_cosets(H)
    print("|X| =", len(X))

    nn = 2*n
    W0 = Algebraic.Sp(nn-2).get_weyl()
    W = Algebraic.Sp(nn).get_weyl()
    hom = mulclose_hom(gen, W.gen)
    for a in gen:
     for b in gen:
        assert hom[a*b] == hom[a]*hom[b]

    if argv.gom == "left":
        a, b = W0.gen
        c, d, e = W.gen
        gom = mulclose_hom([a, b], [c, d*e*d])
    else:
        gom = mulclose_hom(W0.gen, W.gen[1:])
    for a in gom:
     for b in gom:
        assert gom[a*b] == gom[a]*gom[b]

    #for h in H:
    #    g = hom[h]
    #    print(g.shortstr(), '\n')

    smap = SMap()
    for row, coset in enumerate(X):
        #print([(g[a],g[b]) for g in coset])
        for i, g in enumerate(coset):
            g = hom[g]
            g = g.transpose()
            smap[(i+2)*(nn+2), row*(nn+3)] = g.shortstr()

    cols = [0]*len(X)
    for g0 in W0:
        h = gom[g0]
        for row, coset in enumerate(X):
            for i, g in enumerate(coset):
                if h == hom[g]:
                    smap[cols[row]+1, row*(nn+3)+1] = g0.transpose().shortstr()
                    smap[(i+2)*(nn+2), (row+1)*(nn+3)-3] = "*"
                    cols[row] += nn+2

    print(smap)
    print(len(W0))
    print()


def test_building_0():

    B1 = Algebraic.Sp(2).get_weyl()
    B2 = Algebraic.Sp(4).get_weyl()
    B3 = Algebraic.Sp(6).get_weyl()

    if 1:
        hom = mulclose_hom(B2.gen, B3.gen[1:])
    else:
        a, b, c = B3.gen
        a, b = a, b*c*b
        assert b*b == B3.I
        assert a*b != B3.I
        assert a*b*a*b != B3.I
        assert a*b*a*b*a*b != B3.I
        assert a*b*a*b*a*b*a*b == B3.I
        hom = mulclose_hom(B2.gen, [a, b])
    assert len(hom) == len(B2)

    for g in B2:
        h = hom[g]

        smap = SMap()
        smap[1,1] = g.shortstr()

        smap[0,8] = h.shortstr()
        
        print(smap)
        print()

    return

def test_building():
    n = argv.get("n", 3)
    nn = 2*n
    p = argv.get("p", 2)
    G = Algebraic.Sp(nn, p)
    I = G.I
    F = G.invariant_form
    N = len(G)
    building = Building(G)

    print("|G| =", N)

    W = building.W
    print("|W| =", len(W))

    B = building.B
    print("|B| =", len(B))

    if argv.slow:
        found = {}
        for b1 in B:
         for w in W:
          b1w = b1*w
          for b2 in B:
            g = b1w*b2
            path = found.get(g)
            if path is None:
                found[g] = b1, w, b2
                continue
            b11, ww, b22 = path
            lhs = (len(b1.name), len(b2.name))
            rhs = (len(b11.name), len(b22.name))
            if lhs < rhs:
                found[g] = b1, w, b2
        assert len(found) == len(G)
    else:
        found = None
    
    lookup = B.lookup
    for trial in range(100):
        g = G.sample()
        print(g)
        b1, b2 = building.decompose(g)
        w = b1*g*b2
        assert b1 in B
        assert b2 in B
        assert w in W
        if found is not None:
            assert w == found[g][1]

    return

    W = building.W
    found = set()
    for w in W:
        g = w[4:, :]
        s = g.shortstr()
        if s not in found:
            found.add(s)
    
    found = list(found)
    found.sort(key=str)
    for s in found:
        print(s, '\n')


def test_stabilizer():

    n = argv.get("n", 3)
    nn = 2*n
    p = argv.get("p", 2)
    G = Algebraic.Sp(nn, p)

    I = G.I

    #POINT = Matrix([[0,0,0,1]])
    #LINE = Matrix([[0,0,1,0],[0,0,0,1]])

    print("|G| =", len(G))

    B = G.get_borel()
    #P = G.right_stabilizer(LINE)
    if n==3:
        LINE = parse("111111 111111 ..1111 ..1111  ..1111  ..1111")
    else:
        LINE = parse("1111 1111 ..11 ..11")
    P = get_subgroup(G, LINE)
    P = Algebraic(P)

    print("LINE:", len(G)//len(P))
    print(P.show())

    for b in B:
        assert b in P

    W = G.get_weyl()
    cosets = set()
    bruhat = {}
    for w in W:
        C = set(p1*w*p2 for p1 in P for p2 in B)
        C = list(C)
        C.sort()
        C = tuple(C)
        cosets.add(C)
        bruhat.setdefault(C, []).append(w)
        print("[%d]"%(len(cosets)), end="", flush=True)
    print()
    print("double cosets:", len(cosets))
    cosets = list(cosets)
    cosets.sort(key = len)
    size = 0
    for C in cosets:
        size += len(C)
        #for D in cosets:
        #    if C is D:
        #        continue
        #    for c in C:
        #        for d in D:
        #            assert c != d
    print("cosets:")
    for C in cosets:
        print("size:", len(C))
        smap = SMap()
        ws = bruhat[C]
        ws.sort()
        for i, w in enumerate(ws):
            smap[0, i*(nn+2)] = w.shortstr()
        print(smap)
        print()

    assert size==len(G)

    return

    P = G.right_stabilizer(LINE)
    print("LINE:", len(G)//len(P))
    print(P.show())


    return


    n = argv.get("n", 4)
    p = argv.get("p", 2)
    G = Algebraic.GL(n, p)

    I = G.I

    POINT = Matrix([[0,0,0,1]])
    LINE = Matrix([[0,0,1,0],[0,0,0,1]])
    PLANE = Matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1]])

    print("|G| =", len(G))

    P = G.right_stabilizer(POINT)
    print("POINT:", len(G)//len(P))
    print(P.show())

    P = G.right_stabilizer(LINE)
    print("LINE:", len(G)//len(P))
    print(P.show())

    P = G.right_stabilizer(PLANE)
    print("PLANE:", len(G)//len(P))
    print(P.show())

    return

    n = argv.get("n", 3)
    nn = 2*n
    p = argv.get("p", 2)
    G = Algebraic.Sp(nn, p)
    I = G.I
    F = G.invariant_form
    N = len(G)

    POINT = Matrix([[0,0,0,0,0,1]])
    LINE = Matrix([[0,0,0,0,1,0],[0,0,0,0,0,1]])
    PLANE = Matrix([[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])

    print("|G| =", len(G))

    H = G.right_stabilizer(POINT)
    print("POINT:", len(G)//len(H))
    print(H.show())

    H = G.right_stabilizer(LINE)
    print("LINE:", len(G)//len(H))
    print(H.show())

    H = G.right_stabilizer(PLANE)
    print("PLANE:", len(G)//len(H))
    print(H.show())



def test_partial_flag():
    n = argv.get("n", 3)
    nn = 2*n
    p = argv.get("p", 2)
    G = Algebraic.Sp(nn, p)
    I = G.I
    F = G.invariant_form
    N = len(G)

    POINT = parse("111111 .11111 .11111 .11111  .11111  .11111")
    LINE = parse("111111 111111 ..1111 ..1111  ..1111  ..1111")
    PLANE = parse("111111 111111 111111 ...111  ...111  ...111")

    H = get_subgroup(G, POINT)
    print("points:", N//len(H))

    H = get_subgroup(G, LINE)
    print("lines:", N//len(H))

    H = get_subgroup(G, PLANE)
    print("planes:", N//len(H))

    H = get_subgroup(G, POINT * LINE)
    print("point on line:", N//len(H))

    H = get_subgroup(G, POINT * PLANE)
    print("point on plane:", N//len(H))

    H = get_subgroup(G, LINE * PLANE)
    print("line on plane:", N//len(H))

    H = get_subgroup(G, POINT * LINE * PLANE)
    print("point on line on plane:", N//len(H))

    FLAG = "111111 .11111 ..1111 ...111  ....11  .....1"
    H = get_subgroup(G, FLAG)
    print("flags:", N//len(H))


CHECK = argv.get("check", False)

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
        self._str = b' '.join(A.tobytes() for A in items)
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

    def shortstr(self):
        items = ["[%s]"%shortstr(A).replace('\n', ' ') for A in self.items]
        return "[%s]"%(' '.join(items))

    def __hash__(self):
        return hash(self._str)

    def __eq__(self, other):
        return self._str == other._str

    def __ne__(self, other):
        return self._str != other._str

    def __lt__(self, other):
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
        p = self.p
        items = [(numpy.dot(B, A))%p for B in self.items]
        return Figure(items, p)
        
#    def __mul__(self, g):
#        assert isinstance(g, Matrix)
#        A = g.A
#        p = self.p
#        items = [(numpy.dot(A, B))%p for B in self.items]
#        return Figure(items, p)
#
#    def transpose(self):
#        items = [A.transpose() for A in self.items]
#        return Figure(items, self.p)
#
#    def is_zero(self):
#        for A in self.items:
#            if not A.is_zero():
#                return False
#        return True

    def intersect(lhs, rhs):
        "find intersect'ion numbers"
        rows = [[len(intersect(l, r))
            for r in rhs.items]
            for l in lhs.items]
        rows = numpy.array(rows)
        return rows
        
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
            items.append(list(qchoose(d0, d)))
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

    def __len__(self):
        return len(self.figures)

    def __getitem__(self, idx):
        return self.figures[idx]

    def __eq__(self, other):
        return self.figures == other.figures

    def __ne__(self, other):
        return self.figures != other.figures

    def __hash__(self):
        return hash(self.figures)



def test_orbit():

    figures = set()
    for fig in Figure.qchoose([4,3,1]):
        assert fig not in figures
        figures.add(fig)
    assert len(figures) == 105

    n = argv.get("n", 3)
    G = Algebraic.SL(n)
    print("|G| =", len(G))

    left = argv.get("left", [n,1]) 
    right = argv.get("right", [n,1]) 

    orbits = set()
    left = list(Figure.qchoose(left))
    right = list(Figure.qchoose(right))
    print("left:", len(left))
    print("right:", len(right))
    print("figures:", len(left)*len(right))
    found = set()
    for p in left:
      for q in right:
        fig = p+q
        if fig in found:
            continue

        # METHOD 1 ---------------------
#        orbit = set()
#        for g in G:
#            gfig = g*fig
#            orbit.add(g*fig)

        # METHOD 2 ---------------------
        orbit = set([fig])
        bdy = set(orbit)
        while bdy:
            _bdy = set()
            for fig in bdy:
              for g in G.gen:
                gfig = g*fig
                if gfig in orbit:
                    continue
                _bdy.add(gfig)
                orbit.add(gfig)
            bdy = _bdy
        # END METHODS ---------------------

        found.update(orbit)
        print("orbit:", len(orbit), end=", ", flush=True)
        orbit = Orbit(orbit)
        orbits.add(orbit)
    print()
    orbits = list(orbits)
    orbits.sort(key = len)
    print([len(orbit) for orbit in orbits])
    print([len(G)//len(orbit) for orbit in orbits])
    #n = len(orbits[0])
    #print([len(orbit)//n for orbit in orbits])
    print("orbits:", len(orbits))
    #for orbit in orbits:
    #    assert len(orbit)%n == 0


def test_intersect_gl():
    from bruhat import gset
    n = argv.get("n", 4)

    if n==4:
        point = Figure.get_atomic(1, n)
        line = Figure.get_atomic(2, n)
        plane = Figure.get_atomic(3, n)
        vol = Figure.get_atomic(4, n)
    
        figures = [
            point,
            line,
            plane,
            point + line,
            line + plane,
            point + plane,
            point + line + plane + vol,
        ]

    elif n==3:
        point = Figure.get_atomic(1, n)
        line = Figure.get_atomic(2, n)
        plane = Figure.get_atomic(3, n)
        figures = [point, line, point + line + plane]

    else:
        return

    G = Algebraic.GL(n)
    flag = figures[-1]

    bag = set()
    for g in G:
        other = g*plane
        M = flag.intersect(other)
        bag.add(str(M))
    for s in bag:
        print(s)
    print(len(bag))


def test_intersect_sp():

    n = 2
    nn = 2*n
    G = Algebraic.Sp(nn)
    
    items = []
    for m in range(1, nn+1):
        A = numpy.zeros((m, nn), dtype=int)
        I = numpy.identity(m, dtype=int)
        A[:, nn-m:] = I
        items.append(A)
    flag = Figure(items)
    
    if 0:
        bag = {}
        for m in G.grassmanian(2):
            fig = Figure([m.A])
            M = fig.intersect(flag)
            bag[str(M)] = m

    bag = {}
    for g in G:
        other = g*flag
        M = flag.intersect(other)
        bag[str(M)] = M
    keys = list(bag.keys())
    keys.sort(reverse=True)
    for s in keys:
        print(s)
        print()
    print(len(bag))

    from bruhat.poset import Poset
    pairs = []
    for l in keys:
     for r in keys:
        lhs, rhs = bag[l], bag[r]
        if numpy.alltrue(lhs <= rhs):
            pairs.append((l, r))
    poset = Poset(keys, pairs)
    poset.get_dot("order.dot")


def test_g2():
    # FAIL 
    n = 8

    # Wilson (4.27)
    desc = parse("""
    ....1234
    ..12..56
    .1.3.5.7
    1..4.67.
    .23.5..8
    2.4.6.8.
    34..78..
    5678....
    """)
    struct = numpy.zeros((n, n, n), dtype=int)
    for i in range(n):
      for j in range(n):
        k = desc[i, j]
        if k==0:
            continue
        k -= 1
        struct[i, j, k] = 1

    def mul(left, right):
        result = numpy.zeros((n,), dtype=int)
        for i in range(n):
         for j in range(n):
            result += left[i]*right[j]*struct[i,j,:]
        result = Matrix(result)
        return result

    basis = [[0]*n for i in range(n)]
    for i in range(n):
        basis[i][i] = 1
    basis = [Matrix(item) for item in basis]
    x1, x2, x3, x4, x5, x6, x7, x8 = basis
    zero = Matrix([0]*n)
    assert mul(x1, x4) == zero
    assert mul(x4, x1) == x1
    assert mul(mul(x1, x4), x7) == zero
    lhs, rhs = mul(x1, mul(x4, x7)), mul(x1, x7)
    assert lhs == rhs == x3

    e = x4 + x5
    for x in basis:
        assert mul(x, e) == x
        assert mul(e, x) == x

    def test_maufang(x, y, z):
        lhs = mul(mul(x, y), mul(z, x))
        rhs = mul(mul(x, mul(y, z)), x)
        assert lhs == rhs
        lhs = mul(x, mul(y, mul(x, z)))
        rhs = mul(mul(mul(x, y), x), z)
        assert lhs == rhs
        lhs = mul(mul(mul(y, x), z), x)
        rhs = mul(y, mul(x, mul(z, x)))
        assert lhs == rhs

    if 0:
        # Maufang laws work
        for x in basis:
         for y in basis:
          for z in basis:
            test_maufang(x, y, z)
    
        def mk_rand():
            x = [0]*n
            for i in range(n//2):
                x[randint(0, n-1)] = 1
            return Matrix(x)
    
        for trial in range(100):
            x = mk_rand()
            y = mk_rand()
            z = mk_rand()
            test_maufang(x, y, z)

    # See: Wilson 2009
    A = {7:[7,1], 8:[8,2]}
    B = {6:[6,1], 8:[8,3]}
    C = {4:[4,1], 5:[5,1], 6:[6,2], 7:[7,3], 8:[8,5,4,1]}
    D = {3:[3,1], 4:[4,2], 5:[5,2], 7:[7,4,5,2], 8:[8,6]}
    E = {3:[3,2], 7:[7,6]}
    F = {2:[2,1], 4:[4,3], 5:[5,3], 6:[6,4,5,3], 8:[8,7]}

    gens = []
    for op in A,B,C,D,E,F:
        a = numpy.zeros((n, n), dtype=int)
        for col in range(n):
            rows = op.get(col, [col])
            a[col, rows] = 1
        a = Matrix(a)
        gens.append(a)

    #for op in gens:
    #    print( op*e*op == e )

    A, B, C, D, E, F = gens
    I = Matrix.identity(n)

    def is_hom(op):
        #assert op*op == I
        for left in basis:
         for right in basis:
            lhs = mul(left*op, right*op)
            rhs = op*mul(left, right)
            if lhs != rhs:
                return False
        return True

    triangle = []
    for i in range(n):
     for j in range(i+1, n):
        triangle.append((i, j))
    print(triangle)
    for idxs in choose(triangle, 3):
        a = numpy.identity(n, dtype=int)
        for idx in idxs:
            a[idx] = 1
        op = Matrix(a)
        if is_hom(op):
            print("found!")
    print()

    # The unipotent subgroup (== Borel)
    U = Algebraic(gens)
    assert len(U) == 64

    r = Matrix.perm([i-1 for i in [1,3,2,4,5,7,6,8]])
    s = Matrix.perm([i-1 for i in [2,1,6,5,4,3,8,7]])

    # The Weyl group
    W = Algebraic([r, s])
    assert len(W) == 12

    G = mulclose( [B, C, D, E], verbose=True )
    print(len(G))

    return

    gens += [s]
    G = Algebraic(gens, verbose=True)
    print(len(G))
    

# https://mathoverflow.net/questions/270781/what-finite-simple-groups-we-can-obtain-using-octonions
def test_octonion():

    # FAIL 

    # Wilson (4.18)
    desc = """
    .361542
    3.40265
    64.5130
    105.624
    5216.03
    46320.1
    250431.
    """.strip().split()
    assert len(desc) == 7

    n = 8
    struct = numpy.zeros((n, n, n), dtype=int)
    for i in range(n-1):
      for j in range(n-1):
        k = desc[i][j]
        if k=='.':
            k = 0
        else:
            k = int(k)+1
        struct[i+1, j+1, k] = 1
    for i in range(n):
        struct[0, i, i] = 1
        struct[i, 0, i] = 1

    def mul(left, right):
        result = numpy.zeros((n,), dtype=int)
        for i in range(n):
         for j in range(n):
            result += left[i]*right[j]*struct[i,j,:]
        result = Matrix(result)
        return result

    basis = [[0]*n for i in range(n)]
    for i in range(n):
        basis[i][i] = 1
    basis = [Matrix(item) for item in basis]

    def lmul(left):
        rows = []
        for x in basis:
            rows.append(mul(left, x).A)
        A = numpy.array(rows)
        A = Matrix(A)
        return A

    e, x0, x1, x2, x3, x4, x5, x6 = basis
    zero = Matrix([0]*n)

    algebra = []
    for bits in cross([(0,1)]*n):
        v = numpy.array(bits)
        v = Matrix(v)
        algebra.append(v)

    assert mul(x0, x5) == x4
    for x in algebra:
        assert mul(e, x) == x
        assert mul(x, e) == x

    gen = []
    for x in algebra:
      for y in algebra:
        if mul(x, y) != mul(y, x):
            continue
        L = lmul(x)
        R = lmul(y).transpose()
        M = L*R*L*R 
        a = numpy.linalg.det(M.A)
        if a==1 and M*M != M:
            print(M, a)
            gen.append(M)

    gen = []
    for x in basis:
        L = lmul(x)
        R = lmul(x).transpose()
        gen.append(L)
        gen.append(R)
    G = mulclose(gen)
    print(len(gen))


def test_bruhat():
    from bruhat import gset
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

    G = Algebraic.GL(n)

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
            perms.append(gset.Perm(perm))
        tgt = gset.Group(perms)
        #send_perms = [tgt.lookup[perm] for perm in perms]
        #assert send_perms == list(range(len(send_perms)))
        send_perms = list(range(len(perms)))
        X = gset.GSet(repG, tgt, send_perms)
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
                K = gset.Group([g*h*~g for h in H])
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
    
def boost(g):
    A0 = g.A
    n = len(A0)+1
    A = zeros2(n,n)
    A[:n-1, :n-1] = A0
    A[n-1, n-1] = 1
    return Matrix(A)



def test_ASp4():
    n = 4
    G = Algebraic.Sp(n)
    assert len(G) == 720

    # now make affine Sp
    gen = []
    for g in G.gen:
        g1 = boost(g)
        gen.append(g1)

    for i in range(n):
        #I = Matrix.identity(n+1)
        A = identity2(n+1)
        A[i, n] = 1
        M = Matrix(A)
        gen.append(M)
    G = mulclose(gen)
    assert len(G) == 11520

    def get_weight(w):
        pts = []
        for idxs in choose(list(range(n)), w):
            v = zeros2(n+1, 1)
            v[n] = 1
            for idx in idxs:
                v[idx] = 1
            v = Matrix(v)
            pts.append(v)
        return pts
    pts = get_weight(1) + get_weight(3)
    assert len(pts) == 8
    pts = set(pts)
    H = []
    for g in G:
        qts = {g*v for v in pts}
        if qts == pts:
            H.append(g)
    #assert len(H) == 96 
    print(len(H))

    def get_key(pts):
        pts = [str(p) for p in pts]
        pts.sort()
        pts = ''.join(pts)
        return pts

    orbit = {get_key(pts):pts}
    bdy = [pts]
    while bdy:
        _bdy = []
        for g in gen:
          for pts in bdy:
            qts = [g*v for v in pts]
            key = get_key(qts)
            if key in orbit:
                continue
            orbit[key] = pts
            _bdy.append(qts)
        bdy = _bdy
        print(len(orbit), end=" ", flush=True)
    print()
    print("orbit:", len(orbit))
    #assert N % len(orbit) == 0

    # check even intersection
    ptss = list(orbit.values())
    vals = set()
    for p in ptss:
        for q in ptss:
            pq = [u for u in p if u in q]
            assert len(pq)%2 == 0
            vals.add(len(pq))

    return

    # permutation representation -------------------

    from bruhat.action import Perm, Group
    items = list(range(2**n))
    perm = lambda vs : Perm(vs, items)
    gen = [
        perm({0:0, 1:1, 2:14, 3:3, 4:4, 5:5, 6:11, 7:8, 8:7, 9:9, 10:13, 11:6, 12:12, 13:10, 14:2, 15:15}),
        perm({0:0, 1:1, 2:2, 3:10, 4:4, 5:11, 6:6, 7:9, 8:8, 9:7, 10:3, 11:5, 12:12, 13:13, 14:15, 15:14}),
        perm({0:0, 1:1, 2:3, 3:10, 4:11, 5:4, 6:6, 7:9, 8:8, 9:12, 10:2, 11:5, 12:7, 13:15, 14:13, 15:14}),
        perm({0:0, 1:1, 2:3, 3:2, 4:7, 5:9, 6:8, 7:4, 8:6, 9:5, 10:10, 11:12, 12:11, 13:15, 14:14, 15:13}),
        perm({0:0, 1:2, 2:1, 3:4, 4:3, 5:5, 6:14, 7:7, 8:15, 9:11, 10:12, 11:9, 12:10, 13:13, 14:6, 15:8}),
        perm({0:1, 1:0, 2:3, 3:2, 4:4, 5:6, 6:5, 7:7, 8:9, 9:8, 10:14, 11:11, 12:12, 13:15, 14:10, 15:13})]

    G = Group.generate(gen)
    assert len(G) == 11520

    # G is doubly transitive
    pairs = [(i, j) for i in items for j in items if i!=j]
    for trial in range(20):
        src = choice(pairs)
        tgt = choice(pairs)
        for g in G:
            if g(src[0]) == tgt[0] and g(src[1]) == tgt[1]:
                break
        else:
            assert 0, (src, tgt)
    print("OK")
            

    H = G.stabilizer(0)
    assert len(H) == 720 # Sp(4,2)
    
    H1 = G.stabilizer(0, 1)
    assert len(H1) == 48

    #for i in range(2, 2**n):
    #    print(len(G.stabilizer(0, 1, i)))
    #return

    A = parse("""
    [[0 0 0 0 0 1 1 1 0 0 1 0 1 0 1 0]
     [0 0 1 0 0 1 0 0 1 0 0 1 1 1 0 0]
     [0 0 0 1 1 1 0 1 1 0 0 0 0 0 0 1]
     [0 0 1 0 1 0 1 1 0 1 0 0 0 1 0 0]
     [0 1 0 0 0 1 0 0 0 1 1 0 0 1 0 1]
     [1 0 1 1 0 1 0 0 0 1 0 0 0 0 1 0]]
    """)

    def stabilizer(v):
        bits = [i for i in range(2**n) if v[i]]
        perms = []
        for g in G:
            _bits = [g[i] for i in bits]
            _bits.sort()
            if _bits == bits:
                perms.append(g)
        H = Group(perms, items)
        return H
    v = (A[0]+A[1])%2
    H = stabilizer(v)
    assert len(H) == 384

    for i in range(2**n):
        print(len(H.stabilizer(i)), end=" ")
    print()

    H = stabilizer(A[0])
    assert len(H) == 720
    #for g in H:
    #    print(g.fixed(), end=" ")
    #print()


def test_ASp6():
    # try to generalize test_ASp4, fail...

    n = 6
    G = Algebraic.Sp(n)
    N = len(G)
    assert N == 1451520

    # now make affine Sp
    gen = []
    for g in G.gen:
        g1 = boost(g)
        gen.append(g1)

    for i in range(n):
        #I = Matrix.identity(n+1)
        A = identity2(n+1)
        A[i, n] = 1
        M = Matrix(A)
        gen.append(M)
    #G = mulclose(gen, verbose=True)
    #assert len(G) == 11520

    N *= 2**n

    def get_weight(w):
        pts = []
        for idxs in choose(list(range(n)), w):
            v = zeros2(n+1, 1)
            v[n] = 1
            for idx in idxs:
                v[idx] = 1
            v = Matrix(v)
            pts.append(v)
        return pts

    # build a stabilizer
    #pts = get_weight(1) + get_weight(n-1)
    pts = get_weight(2) + get_weight(n-2)
    #pts = get_weight(2)
    #assert len(pts) == 2*n

    def get_key(pts):
        pts = [str(p) for p in pts]
        pts.sort()
        pts = ''.join(pts)
        return pts

    orbit = {get_key(pts):pts}
    bdy = [pts]
    while bdy:
        _bdy = []
        for g in gen:
          for pts in bdy:
            qts = [g*v for v in pts]
            key = get_key(qts)
            if key in orbit:
                continue
            orbit[key] = pts
            _bdy.append(qts)
        bdy = _bdy
        print(len(orbit), end=" ", flush=True)
    print()
    print("orbit:", len(orbit))
    assert N % len(orbit) == 0

    # look for even intersection & fail
    ptss = list(orbit.values())
    vals = set()
    #for p in ptss:
    #    for q in ptss:
    for trial in range(1000):
        p = choice(ptss)
        q = choice(ptss)
        pq = [u for u in p if u in q]
        vals.add(len(pq))
    print(vals) # these need to be even... fail

    return orbit




def main():

    n = argv.get("n", 2)
    p = argv.get("p", 2)

    print("|GL(%d, %d)| = %d"%(p, n, order_gl(n, p)))
    print("|SL(%d, %d)| = %d"%(p, n, order_sl(n, p)))
    if n%2 == 0:
        print("|Sp(%d, %d)| = %d"%(p, n, order_sp(n, p)))

    if argv.SL:
        G = Algebraic.SL(n, p)
        print("|G| =", len(G))

    if argv.GL:
        G = Algebraic.GL(n, p)
        print("|G| =", len(G))



if __name__ == "__main__":
    from time import time
    start_time = time()
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))

    

    




