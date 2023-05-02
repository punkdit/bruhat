#!/usr/bin/env python

"""
orthogonal subspaces of Z/2 vector spaces
"""


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


from bruhat.action import mulclose
from bruhat.argv import argv
from bruhat.solve import parse, enum2, span, shortstr, rank, shortstrx, zeros2, dot2, array2
from bruhat.util import cross, allperms, choose_rev
from bruhat.smap import SMap
from bruhat.qcode import QCode
from bruhat.oeqc import gap_matrix

from bruhat import algebraic
#scalar = numpy.uint8  # numba doesn't like it
scalar = numpy.int8  # <---------  Careful with this, 
                     # <---------  matrices should have every dimension < 128.
algebraic.scalar = scalar
from bruhat.algebraic import Matrix


def find_kernel(A, inplace=False, check=False, verbose=False):
    """return a list of vectors that span the nullspace of A
    """

    if check:
        A0 = A.copy() # save

#    L, U = lu_decompose(A)
#    assert eq2(numpy.dot(L, U), A)

    U = row_reduce(A, inplace=inplace)

    # We are looking for a basis for the nullspace of A

    m, n = U.shape

    if verbose:
        print("find_kernel: shape", m, n)
        #print shortstr(U, deco=True)
        print()

    items = []
    for row in range(m):
        cols = numpy.where(U[row, :])[0]
        if not len(cols):
            break
        col = cols[0]
        items.append((row, col))

    #items.sort(key = lambda item : item[1])
    #print items
    #rows = [row for (row, col) in items]
    #U = U[rows]
    leading = [col for (row, col) in items]
    degeneracy = m - len(leading)

    if verbose:
        print("leading:", leading)
        print("degeneracy:", degeneracy)

    # Look for the free variables
    vars = []
    row = 0
    col = 0
    while row < m and col < n:
        #print row, col
        if U[row:, col].max() == 0: # XXX optimize this XXX
            #print "*"
            assert U[row:, col].max() == 0, U[row:, col]
            vars.append(col)
        else:
            #print U[row:, col]
            while row<m and U[row:, col].max():
                row += 1
                #print "row", row
                #if row<m:
                #    print U[row:, col]
        col += 1
    for k in range(col, n):
        vars.append(k)

    if verbose:
        print("found %d free vars:" % len(vars), vars)

    basis = []
    for var in vars:

        #print "var:", var
        v = numpy.zeros((n,), dtype=scalar)
        v[var] = 1
        row = min(var-1, m-1)
        while row>=0:
            u = numpy.dot(U[row], v)
            if u.sum()%2:
                col = leading[row]
                #print "\trow", row, "leading:", col
                v[col] = 1
                #print '\t', shortstr(v)
            assert numpy.dot(U[row], v).sum()%2==0, row
            row -= 1
        #print '\t', shortstr(v)
        if check:
            assert numpy.dot(A0, v).sum()%2 == 0, shortstr(v)
        basis.append(v)

    K = numpy.array(basis, dtype=scalar)
    if not basis:
        K.shape = (0, A.shape[1])
    else:
        assert K.shape[1] == A.shape[1]

    return K


def intersect(W1, W2):
    "find maximal subspace within rowspace W1 & rowspace W2"
    W = numpy.concatenate((W1, W2))
    K = find_kernel(W.transpose())
    W = dot2(K[:, :len(W1)], W1)
    W = row_reduce(W)
    return W


def complement(H):
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


def get_perms(n):
    gen = []
    for i in range(n-1):
        items = list(range(n))
        items[i:i+2] = i+1, i
        M = Matrix.perm(items)
        gen.append(M)
    return gen


def get_SO_gens_odd(n):
    assert n%2 == 1
    assert is_orthogonal(offdiag(n-1))
    assert is_orthogonal(antidiag(n))

    gen = get_perms(n)
    for k in range(1, n//2+1, 2):
        A = antidiag(n)
        A[k:,:-k] = offdiag(n-k)
        M = Matrix(A)
        #print("gen:")
        #print(M)
        assert is_orthogonal(M)
        gen.append(M)

    return gen


def get_SO_gens(n):
    if n%2:
        return get_SO_gens_odd(n)
    
    gen = get_perms(n)
    A = parse("""
    [[0 1 1 0 0 1]
     [1 0 1 1 1 1]
     [1 1 0 0 0 1]
     [1 1 1 1 1 0]
     [0 1 0 0 1 1]
     [0 1 0 1 0 1]]
    """)
    B = numpy.identity(n, dtype=int)
    n0 = (n-6)//2
    B[n0:n0+6, n0:n0+6] = A
    assert is_orthogonal(B)
    M = Matrix(B)
    gen.append(M)
    if argv.gap:
        gap_code(gen)
    return gen


def test_gap():
    n = argv.get("n", 5)
    gen = get_SO_gens(n)
    gap_code(gen)
    G = mulclose(gen)
    print(len(G))

def gap_code(gen):
    vs = []
    for i, g in enumerate(gen):
        print("m%d := %s;;"%(i, gap_matrix(g)))
        vs.append("m%d"%i)
    print("G := Group(%s);" % (','.join(vs)))


def find_gens_SO6():
    n = 6
    gen = get_perms(n)
    gap_code(gen)
    G = mulclose(gen, verbose=True)
    #assert len(G) == 23040

    I = numpy.identity(n, dtype=scalar)
    while 1:
        items = [randint(0,1) for i in range(n*n)]
        A = array2(items)
        A.shape = (n,n)
        #print(shortstr(A))
        B = dot2(A, A.transpose())
        if numpy.alltrue(I==B):
            M = Matrix(A)
            gen.append(M)
            G = mulclose(gen, verbose=True)
            if len(G) == 23040:
                break
    print()
    print(M)
    gap_code(gen)
        

def test_dickson():
    n = argv.get("n", 5)
    I = Matrix.identity(n)
    gen = get_SO_gens(n)
    G = mulclose(gen, verbose=True)
    print(len(G))
    SO = []
    for g in G:
        A = g-I
        #print(g)
        #print(A)
        #print()
        if rank(A.A)%2==0:
            SO.append(g)
    print(len(SO))
    SO = mulclose(SO, verbose=True)
    print(len(SO))


def find_orbit_fast(gen, M, verbose=False):
    print("find_orbit")
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
            #gM = gM.normal_form()
            gM = normal_form(gM.A)
            gM = Matrix(gM)
            key = gM.key[1]
            #key = gM.get_bitkey()
            if key in orbit:
                continue
            _bdy.add(gM)
            orbit.add(key)
            #print(gM)
            yield gM
        bdy = _bdy
    if verbose:
        print()
    #return orbit

def find_orbit(gen, M, verbose=False):
    print("find_orbit")
    #orbit = set([M.key[1]])
    orbit = set([M.get_bitkey()])
    bdy = set([M])
    yield M
    while bdy:
        if verbose:
            print("(%s:%s)"%(len(bdy),len(orbit)), end="", flush=True)
        _bdy = set()
        for M in bdy:
          for g in gen:
            gM = M*g
            #w = gM*gM.transpose()
            #assert w.is_zero()
            #gM = gM.normal_form()
            gM = normal_form(gM.A)
            gM = Matrix(gM)
            if gM in bdy or gM in _bdy: # hmm does not seem to help much..
                continue
            #key = gM.key[1]
            key = gM.get_bitkey() # slower
            #print(key)
            #print("sizeof=%d"%sys.getsizeof(key), len(key))
            #print(gM.key[1])
            #print("sizeof=%d"%sys.getsizeof(gM.key[1]), len(gM.key[1]))
            #assert 0
            if key in orbit:
                continue
            _bdy.add(gM)
            orbit.add(key)
            #print(gM)
            yield gM
        bdy = _bdy
    if verbose:
        print()
    #return orbit

def find_random(gen, M, verbose=False):
    print("find_random")
    gen = mulclose(gen, maxsize=1000)
    print("find_random", len(gen))
    gen = list(gen)
    while 1:
        g = choice(gen)
        M = M*g
        yield M


def test_12_2_3():
    A = parse("""
    [[1 1 1 1 0 0 0 0 0 0 0 0]
     [1 1 0 0 1 1 0 0 0 0 0 0]
     [0 0 0 0 0 0 1 1 1 1 0 0]
     [0 0 0 0 0 0 1 1 0 0 1 1]
     [1 0 1 0 1 0 1 0 1 0 1 0]]
    """)
    gen = get_SO_gens(12)
    M = Matrix(A)
    count = 0
    for M in find_orbit(gen, M, verbose=True):
        count += 1
    print(count)


def get_logops(H):
    assert isinstance(H, numpy.ndarray)
    m, n = H.shape
    W = complement(H)
    K = find_kernel(H, check=True)
    #print()
    #print("H =")
    #print(H)
    #print("K =")
    #print(K)
    HKt = numpy.dot(H, K.transpose()) % 2
    assert HKt.sum() == 0
    L = intersect(W, K)
    k = len(L)
    assert k+2*m == n
    HL = numpy.concatenate((H, L))
    assert rank(HL) == m+k
    HLt = numpy.dot(H, L.transpose()) % 2
    assert HLt.sum() == 0
    return L


def search_k1(n):
    m = n//2
    p = 2

    find = find_orbit
    if argv.find_orbit_fast:
        find = find_orbit_fast
    if argv.find_random:
        find = find_random

    gen = get_SO_gens(n)
    u = numpy.array([1]*n, dtype=scalar)

    best_d = 1

    assert m<=20, "um..?"
    B = numpy.array(list(numpy.ndindex((2,)*m)), dtype=scalar)

    A = zeros2(m, n)
    for i in range(m):
        A[i,2*i] = 1
        A[i,2*i+1] = 1
    M = Matrix(A)
    print("M =")
    print(M)

    count = 0
    for M in find(gen, M, verbose=True):
        count += 1
        if M.A.sum(0).min() == 0:
            continue # distance == 1
        #d0 = n
        #for v in span(M):
        #    d1 = ((v+u)%2).sum()
        #    if d1<d0:
        #        d0 = d1
        C = numpy.dot(B, M.A)
        C = (C+u)%2
        C = C.sum(1)
        d = numpy.min(C)
        #assert d == d0
        if d > best_d:
            best_d = d
            print()
            print(M, d)
            print()
    print(count, "= [%s]_%d"%(_basep(count, p), p))


def search(n, m):
    k = 2*n-m
    assert k>1
    p = 2

    gen = get_SO_gens(n)
    u = numpy.array([1]*n, dtype=scalar)

    best_d = 1

    assert m<=20, "um..?"
    B = numpy.array(list(numpy.ndindex((2,)*m)), dtype=scalar)

    A = zeros2(m, n)
    for i in range(m):
        A[i,2*i] = 1
        A[i,2*i+1] = 1
    
    A = parse("""
    [[1 1 1 1 0 0 0 0 0 0 0 0]
     [1 1 0 0 1 1 0 0 0 0 0 0]
     [0 0 0 0 0 0 1 1 1 1 0 0]
     [0 0 0 0 0 0 1 1 0 0 1 1]
     [1 0 1 0 1 0 1 0 1 0 1 0]]
    """)
    M = Matrix(A)
    print("M =")
    print(M)

    find = find_orbit
    if argv.find_orbit_fast:
        find = find_orbit_fast # about twice as fast
    if argv.find_random:
        find = find_random

    count = 0
    for M in find(gen, M, verbose=True):
        count += 1
        #if count > 100000:
        #    break
        if M.A.sum(0).min() == 0:
            continue # distance == 1
        #d0 = n
        #for v in span(M):
        #    d1 = ((v+u)%2).sum()
        #    if d1<d0:
        #        d0 = d1
        C = numpy.dot(B, M.A)
        C = (C+u)%2
        C = C.sum(1)
        d = numpy.min(C)
        #assert d == d0
        if d <= best_d:
            continue
        H = M.A
        L = get_logops(H)

        d = n
        for u in span(L):
            if u.sum() == 0:
                continue
            for v in span(H):
                uv = (u+v)%2
                d = min(d, uv.sum())
        if d > best_d:
            best_d = d
            print()
            print(M, d)
            print()
            code = QCode.build_css(H, H)
            params = code.get_params()
            print(params)
            assert params[2] is None or params[2] == d
    
#
#        C = numpy.dot(B, H)
#        for u in span(L):
#            if u.sum() == 0:
#                continue
#            C1 = (C+u)%2
#            d = numpy.min(C1.sum(1))
#            if d > best_d:
#                best_d = d
#                print()
#                print(M, d)
#                print(C1)
#                print(C1.sum(1))
#                print()
#                code = QCode.build_css(H, H)
#                params = code.get_params()
#                print(params)
#                assert params[2] == d
    print(count, "= [%s]_%d"%(_basep(count, p), p))



def main_stabs():

    H = parse("""
[[1 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0]
 [0 1 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0]
 [0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 1]
 [0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1]
 [0 0 0 0 1 0 0 1 0 0 1 0 1 1 1 1 1]
 [0 0 0 0 0 1 0 1 0 0 0 1 1 1 1 1 1]
 [0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0]
 [0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0]]
    """)
    min_stabs(H)


def min_stabs(H):
    print(H, H.sum(1))
    m, n = H.shape
    assert rank(H) == m
    V = list(span(H))
    ws = [(v.sum(),v) for v in V if v.sum()]
    ws.sort(key = lambda w:w[0])
    print([w[0] for w in ws])
    
    for basis in choose_rev(ws, m):
        vs = [b[1] for b in basis]
        A = array2(vs)
        m1 = rank(A)
        assert 0<m1<=m
        if m1==m:
            break
    print(A, A.sum(1))


def main():
    n = argv.get("n", 5)
    k = argv.get("k")
    m = argv.get("m", n//2 if k is None else (n-k)//2)
    k = n-2*m
    assert k>=0
    print("[[%d,%d]]"%(n,k))
    if k==1:
        search_k1(n)
    else:
        search(n, m)


if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))




