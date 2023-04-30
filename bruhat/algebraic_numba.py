#!/usr/bin/env python


import sys, os
from time import time
start_time = time()
import random
from random import randint, choice
from functools import reduce
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import add, mul
from math import prod

import numpy
import numba


from bruhat.action import mulclose, mulclose_hom
from bruhat.spec import isprime
from bruhat.argv import argv
from bruhat.solve import parse, enum2, row_reduce, span, shortstr, rank, shortstrx, pseudo_inverse, intersect
from bruhat.solve import zeros2
from bruhat.dev import geometry
from bruhat.util import cross, allperms, choose
from bruhat.smap import SMap

from bruhat import algebraic
algebraic.scalar = numpy.int8 # ???
from bruhat.algebraic import Matrix


@numba.njit
def swap_row(A, j, k):
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row

@numba.njit
def swap_col(A, j, k):
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col

@numba.njit
def row_reduce(H, truncate=True, inplace=False, check=False):
    """Remove zero rows if truncate=True
    """

    #assert len(H.shape)==2, H.shape
    m, n = H.shape
    orig = H
    if not inplace:
        H = H.copy()

    if m*n==0:
        return H

    i = 0
    j = 0
    while i < m and j < n:
        assert i<=j
        if i and check:
            assert H[i:,:j].max() == 0 # XX rm

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if H[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            swap_row(H, i, i1)

        assert H[i, j]
        for i1 in range(i+1, m):
            if H[i1, j]:
                H[i1, :] += H[i, :]
                H[i1, :] %= 2


        #assert 0<=H.max()<=1, orig

        i += 1
        j += 1
    if truncate:
        m = H.shape[0]-1
        #print "sum:", m, H[m, :], H[m, :].sum()
        while m>=0 and H[m, :].sum()==0:
            m -= 1
        H = H[:m+1, :]

    return H


@numba.njit
def normal_form(A, p=2):
    "reduced row-echelon form"
    assert p==2
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



def _basep(n, p):
    e = n//p
    q = n%p
    if n == 0:
        return '0'
    elif e == 0:
        return str(q)
    else:
        return _basep(e, p) + str(q)

def is_orthogonal(M):
    if isinstance(M, numpy.ndarray):
        M = Matrix(M)
    assert isinstance(M, Matrix)
    n = len(M)
    return M*M.transpose() == Matrix.identity(n,2)

def offdiag(n):
    A = numpy.zeros((n,n),dtype=int)
    A[:] = 1
    A += numpy.identity(n,dtype=int)
    A %= 2
    return A

def antidiag(n):
    A = numpy.zeros((n,n),dtype=int)
    for i in range(n):
        A[i,n-i-1] = 1
    return A


def get_SO_gens(n):
    assert n%2 == 1
    assert is_orthogonal(offdiag(n-1))
    assert is_orthogonal(antidiag(n))

    gen = []
    for i in range(n-1):
        items = list(range(n))
        items[i:i+2] = i+1, i
        M = Matrix.perm(items)
        gen.append(M)

    for k in range(1, n//2+1, 2):
        A = antidiag(n)
        A[k:,:-k] = offdiag(n-k)
        M = Matrix(A)
        #print("gen:")
        #print(M)
        assert is_orthogonal(M)
        gen.append(M)

    return gen


def find_orbit(gen, M, verbose=False):
    orbit = set([M.key[1]])
    bdy = set([M])
    yield M
    while bdy:
        if verbose:
            print("(%s)"%len(bdy), end="", flush=True)
        _bdy = set()
        for M in bdy:
          for g in gen:
            gM = M*g
            #w = gM*gM.transpose()
            #assert w.is_zero()
            gM = gM.normal_form()
            s = gM.key[1]
            if s in orbit:
                continue
            _bdy.add(gM)
            orbit.add(s)
            #print(gM)
            yield gM
        bdy = _bdy
    if verbose:
        print()
    #return orbit



def build_SO():
    n = argv.get("n", 5)
    m = argv.get("m", n//2)
    p = 2

    gen = get_SO_gens(n)
    u = numpy.array([1]*n)

    best_d = 1

    #for m in range(1, n//2+1):
    #m = n//2
    #while m:
    for m in [m]:
        A = zeros2(m, n)
        for i in range(m):
            A[i,2*i] = 1
            A[i,2*i+1] = 1
        M = Matrix(A)
        print("M =")
        print(M)

        count = 0
        for M in find_orbit(gen, M, verbose=True):
            count += 1
            if M.A.sum(0).min() == 0:
                continue
            d = n
            for v in span(M):
                d1 = ((v+u)%2).sum()
                if d1<d:
                    d = d1
            if d > best_d:
                best_d = d
                print()
                print(M, d)
                print()
        print(count, "= [%s]_%d"%(_basep(count, p), p))

    if argv.mulclose:
        G = mulclose(gen, verbose=True)
        print(len(G))



if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))




