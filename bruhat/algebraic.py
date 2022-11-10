#!/usr/bin/env python3

"""
Algebraic groups: matrix groups over Z/pZ.

"""


import sys, os
import random
from random import randint

import numpy

#scalar = numpy.int64
scalar = numpy.int8 # CAREFUL !!

#from bruhat.gset import Group, Perm, GSet
from bruhat import gset
from bruhat.action import mulclose
from bruhat.spec import isprime
from bruhat.argv import argv
from bruhat.solve import parse, enum2, row_reduce, span, shortstr, rank, shortstrx
from bruhat.dev import geometry
from bruhat.util import cross

EPSILON = 1e-8

DEFAULT_P = argv.get("p", 2)


def qchoose_2(n, m, p=DEFAULT_P):
    assert m<=n
    col = m
    row = n-m
    for A in geometry.get_cell(row, col, p):
        yield A


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


def normal_form(A, p=DEFAULT_P):
    "reduced row-echelon form"
    if p!=2:
        return normal_form_p(A, p)
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
        self.key = (self.p, self.A.tobytes())
        self._hash = hash(self.key)
        self.shape = A.shape

    def __str__(self):
        return str(self.A)

    def shortstr(self):
        return shortstr(self.A)

    def __hash__(self):
        return self._hash

    def is_zero(self):
        return self.A.sum() == 0

    def __len__(self):
        return len(self.A)

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

    def __getitem__(self, idx):
        A = self.A[idx]
        return Matrix(A, self.p)

    def transpose(self):
        A = self.A
        return Matrix(A.transpose(), self.p)

    def mask(self, A):
        return Matrix(self.A * A, self.p)

    def normal_form(self):
        A = normal_form(self.A, self.p)
        return Matrix(A, self.p)

    @classmethod
    def all_codes(cls, m, n, p=DEFAULT_P):
        assert p==2
        for A in geometry.all_codes(m, n):
            yield cls(A, p)


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


class Group(object):
    def __init__(self, gen, order=None, p=DEFAULT_P, **kw):
        self.__dict__.update(kw)
        self.gen = list(gen)
        self.order = order
        self.G = None
        self.p = p
        assert gen
        A = gen[0]
        self.n = len(A)
        assert p == A.p

    def get_elements(self):
        if self.G is None:
            self.G = mulclose(self.gen, maxsize=self.order)
            self.order = len(self.G)
        return self.G

    def __len__(self):
        if self.order is None:
            self.get_elements()
        return self.order

    def __getitem__(self, idx):
        if self.G is None:
            self.get_elements()
        return self.G[idx]
    
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
    
        H = cls.SL(n, p)
        gen = list(H.gen)
        for i in range(2, p):
            A = Matrix(i*numpy.identity(n, scalar), p)
            gen.append(A)
        order = order_gl(n, p)
        return cls(gen, order, p=p, **kw)

    # See:
    # Pairs of Generators for Matrix Groups. I
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
        G = Group(gens, order=47377612800, p=p,
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
        G = Group(gens, order=9170703360, p=p,
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
        G = Group(gens, order=1451520, p=p,
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
        G = Group(gens, p=p,
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
        G = Group(gens, p=p,
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
        G = Group(gens, p=p,
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
        G = Group(gens, p=p,
            invariant_bilinear_form = Matrix([[0,1,0],[1,0,0],[0,0,2]], 3),
            invariant_quadratic_form = Matrix([[0,1,0],[0,0,0],[0,0,1]], 3))
        return G

    @classmethod
    def SO_3_2(cls, **kw):
        p = 2
        gens = [[[1,0,0],[1,1,1],[0,0,1]],[[1,0,0],[0,0,1],[0,1,0]]]
        gens = [Matrix(A, p) for A in gens]
        G = Group(gens, p=p,
            invariant_bilinear_form = Matrix([[0,0,0],[0,0,1],[0,1,0]], p),
            invariant_quadratic_form = Matrix([[1,0,0],[0,0,0],[0,1,0]], p))
        return G

    @classmethod
    def make(cls, gens, invariant_bilinear_form, invariant_quadratic_form, order=None, p=DEFAULT_P, **kw):
        gens = [Matrix(A, p) for A in gens]
        G = Group(gens, order, p=p,
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


def basep(n, p):
    e = n//p
    q = n%p
    if n == 0:
        return '0'
    elif e == 0:
        return str(q)
    else:
        return basep(e, p) + str(q)

def test_so():

    n = argv.get("n", 3)
    m = argv.get("m", 1)
    p = argv.get("p", 2)
    e = argv.get("e", 1)

    G = Group.SO(n, p, e)
    print("|G| =", G.order)

    B = G.invariant_bilinear_form
    Q = G.invariant_quadratic_form

    print(B)

    for g in G.gen:
        assert g * B * g.transpose()  == B

    if 0:
        items = []
        #for u in Matrix.all_codes(m, n, p):
        for M in qchoose_2(n, m, p):
            M = Matrix(M, p, shape=(m,n))
            if p>2:  # .. ?.?.?
                w = M * B * M.transpose() # only works when p>2
            else:
                w = M * Q * M.transpose() # use this one for p==2 ... ?
            #print(M, v.is_zero(), w.is_zero())
            if not w.is_zero():
                continue
            items.append(M)
        print(len(items))
        orbit = set([items[0]])

    while 1:
        M = [randint(0, p-1) for i in range(n*m)]
        M = numpy.array(M)
        M.shape = (m, n)
        M = normal_form_p(M, p, truncate=True)
        if len(M) < m:
            continue
        M = Matrix(M, p, shape=(m,n))
        if p>2:  # .. ?.?.?
            w = M * B * M.transpose() # only works when p>2
        else:
            w = M * Q * M.transpose() # use this one for p==2 ... ?
        if w.is_zero():
            break
    print("found")
    print(M)

    orbit = set([M.key[1]])
    bdy = set([M])
    while bdy:
        print("(%s)"%len(bdy), end="", flush=True)
        _bdy = set()
        for M in bdy:
          for g in G.gen:
            gM = M*g
            gM = gM.normal_form()
            s = gM.key[1]
            if s in orbit:
                continue
            _bdy.add(gM)
            orbit.add(s)
        bdy = _bdy
    print()
    N = len(orbit)
    print(N, "[%s]_%d"%(basep(N, p), p))

    if 0:
        for M in orbit:
            w = M*B*M.transpose()
            assert w.is_zero()

    return
    
    ops = build_hecke(G, items, items)
    print("hecke:", len(ops))


class Sp(Group):

    def qchoose(self, m):
        F = self.invariant_form
        p = self.p
        n = self.n
        for M in qchoose_2(n, m, p):
            M = Matrix(M, p)
            A = M*F*M.transpose()
            if A.is_zero():
                yield M

    def all_flags(self, dims, p=DEFAULT_P):
        m = dims[0]
        items = [[M.A for M in self.qchoose(m)]]
        for m1 in dims[1:]:
            assert m1<=m
            items.append(list(qchoose_2(m, m1)))
            m = m1
        n = len(items)
        for select in cross(items):
            A = select[0]
            flag = [A]
            for i in range(n-1):
                A = numpy.dot(select[i+1], A) % p
                flag.append(A)
            flag = list(reversed(flag))
            flag = Figure(flag)
            yield flag

    def test_slow_grassmanian():
    
        p = argv.get("p", 2)
        n = argv.get("n", 3)
        m = argv.get("m", 2)
    
        G = Group.Sp(2*n, p)

    def slow_grassmanian(self, m):
        n = self.n//2
        p = self.p
        F = self.invariant_form
    
        bits = 2*n*m
        found = {}
    
        A = numpy.zeros((m, 2*n), dtype=int)
        jdxs = list(numpy.ndindex(A.shape))
        assert len(jdxs) == bits
        for idxs in cross( [tuple(range(p))]*bits ):
            for (bit, jdx) in zip(idxs, jdxs):
                A[jdx] = bit
            M = Matrix(A)
            B = M*F*M.transpose()
            if not B.is_zero():
                continue
            M = M.normal_form()
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

    def is_symplectic(self, M):
        F = self.invariant_form
        MM = M*F*M.transpose()
        return MM.is_zero()


def test_grassmanian():
    # here we are building bigger grassmanian's from smaller ones
    # looking for the Sp(q)-deformed pascal triangle (and so far failing)

    p = argv.get("p", 2)
    n = argv.get("n", 3)
    m = argv.get("m", 2)

    G = Group.Sp(2*n, p)
    G1 = Group.Sp(2*(n+1), p)

    items = list(G.grassmanian(m))
    print(len(items))

    #gr1 = set(G1.grassmanian(m+1))
    #print(len(gr1))

    if 0:
        items1 = set()
        for M in items:
            assert G.is_symplectic(M)
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
                if G1.is_symplectic(M1):
                    M1 = M1.normal_form()
                    found.add(M1)
            items1.update(found)
            print(len(found), end=" ", flush=True)
        print()
        print(len(items1))

    items1 = set()
    for M in items:
        assert G.is_symplectic(M)
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
            if not G1.is_symplectic(M1):
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

        


def test_symplectic():

    n = argv.get("n", 4)
    assert n%2==0
    m = argv.get("m", 1)
    assert m<=n

    p = argv.get("p", 2)

    G = Group.Sp(n, p)
    F = G.invariant_form
    for g in G.gen:
        assert g * F * g.transpose()  == F
        #print(g)

    #gs = list(G.get_elements())
    #print(len(gs))
    #for g in gs:
    #    A = g.A
    #    if numpy.alltrue(A.sum(0) == 1):
    #        print(A)
    #        A = numpy.array(g.A, dtype=float)
    #        a = numpy.linalg.det(A)
    #        print("det:", a)
    #return

    lookup = {}
    items = []
    for idx, v in enumerate(enum2(n)):
        v = numpy.array(v)
        v.shape = (n, 1)
        v = Matrix(v)
        lookup[v] = idx
        items.append(v)
    for g in G.gen:
        for idx, v in enumerate(items):
            u = g*v
            jdx = lookup[u]
            print("%s:%s"%(idx+1,jdx+1), end=" ")
        print()

    #for flag in G.all_flags([2, 1]):
    #    print(flag)
    print(len(G))

    if n==4:
        left = list(G.all_flags([2, 1]))
        print(len(left))
        right = list(G.all_flags([2, 1]))
    elif n==6:
        left = list(G.all_flags([3, 2, 1]))
        right = list(G.all_flags([3, 2, 1]))

    ops = build_hecke(G, left, right)
    print(len(ops))

    return

    items = list(G.qchoose(m))

    items.sort(key = str)
    if argv.show:
        for M in items:
            print(M)
    print(len(items))




    
def test_1():
    n = argv.get("n", 3)
    m = argv.get("m", 2)
    p = argv.get("p", 2)

    G = Group.Sp(2*n, p)
    F = G.invariant_form

    G1 = Group.Sp(2*n + 2, p)
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
    lookup = dict((v.tobytes(), idx) for (idx, v) in enumerate(space))
    perms = []
    for op in G:
        idxs = [lookup[(numpy.dot(op.A, v)%p).tobytes()] for v in space]
        perms.append(gset.Perm(idxs))
    G = gset.Group(perms)
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
            items.append(list(qchoose_2(d0, d)))
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


def test_hecke():

    n = argv.get("n", 3)
    G = Group.SL(n)
    print("|G| =", len(G))

    left = argv.get("left", [n,1]) 
    right = argv.get("right", left)

    left = list(Figure.qchoose(left))
    right = list(Figure.qchoose(right))

    ops = build_hecke(G, left, right)
    print("Hecke operators:", len(ops))

    if argv.eigvals:
      for J in ops:

        print(J.shape, int(round(J.sum())))
        vals = numpy.linalg.eigvals(J)
        #print(vals)
        ss = []
        for x in vals:
            if abs(x.imag)>EPSILON and abs(x.real)>EPSILON:
                ss.append(str(x))
            elif abs(x.imag)>EPSILON:
                ss.append("%.4f"%(x.imag)+"j")
            elif abs(x.real - int(round(x.real))) < EPSILON:
                ss.append("%.0f"%x.real)
            else:
                ss.append("%.4f"%x.real)
        ss = set(ss)
        print(' '.join(ss), '\n')



def build_hecke(G, left, right, verbose=argv.get("verbose", False)):

    if verbose:
        print("left:", len(left))
        print("right:", len(right))
        print("figures:", len(left)*len(right))

    llookup = dict((i, fig) for (fig, i) in enumerate(left))
    rlookup = dict((i, fig) for (fig, i) in enumerate(right))
    m = len(left)
    n = len(right)

    #H = numpy.zeros((m, n))
    remain = set(numpy.ndindex((m, n)))

    ops = []
    while remain:
        if verbose:
            print("[%s:%s]"%(len(ops),len(remain)), end="", flush=True)
        i, j = iter(remain).__next__()
        #assert H[i, j] == 0

        #fig = left[i] + right[j]
        J = numpy.zeros((m, n))

        bdy = set([(i, j)])
        while bdy:
            _bdy = set()
            #print(bdy)
            for (i, j) in bdy:
              l, r = left[i], right[j]
              for g in G.gen:
                i = llookup[g*l]
                j = llookup[g*r]
                if J[i, j]:
                    continue
                _bdy.add((i, j))
                J[i, j] = 1
                remain.remove((i, j))
            bdy = _bdy
        #print(shortstr(J))
        #print()
        ops.append(J)

    ops.sort(key = lambda J : J.sum())
    return ops


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
    G = Group.SL(n)
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
    



def main():

    n = argv.get("n", 2)
    p = argv.get("p", 2)

    print("|GL(%d, %d)| = %d"%(p, n, order_gl(n, p)))
    print("|SL(%d, %d)| = %d"%(p, n, order_sl(n, p)))
    if n%2 == 0:
        print("|Sp(%d, %d)| = %d"%(p, n, order_sp(n, p)))

    if argv.SL:
        G = Group.SL(n, p)
        print("|G| =", len(G))

    if argv.GL:
        G = Group.GL(n, p)
        print("|G| =", len(G))



if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK\n")

    

    




