#!/usr/bin/env python

"""
m by n isotropic matrices correspond to m by 2n+1 biorthogonal matrices ... ??!?!?
"""

from functools import cache
from time import time
start_time = time()

import numpy
from numpy import dot

from bruhat import sp_pascal 
from bruhat import css 
from bruhat.css import qchoose_2
from bruhat.solve import dot2, span, zeros2
from bruhat.argv import argv
from bruhat.smap import SMap

from bruhat.poly import Poly
from bruhat.element import Fraction, Q, Ring, FiniteField


def selfdual_css(n, k, p=2):

    m = (n-k)//2
    #print("selfdual_css", n, m)
    count = 0
    found = 0
    items = []
    for Hx in qchoose_2(n, m, p):
        if (dot(Hx, Hx.transpose())%p).sum() != 0:
            continue
        piv = tuple(list(row).index(1) for row in Hx)
        items.append((piv, Hx))
        count += 1
    return items


def dump(pivs):
    smap = SMap()
    dist = {piv:0 for piv in pivs}
    N = len(pivs[0])
    shape = [0]*N
    for piv in pivs:
        dist[piv] += 1
        for i, p in enumerate(piv):
            shape[i] = max(shape[i], p+1)
    A = numpy.zeros(shape, dtype=int)
    for k,v in dist.items():
        A[k] = v
    print(A)
    print()


def test_pivs(n, m):
    pivs = [tuple(piv) for piv, H in sp_pascal.grassmannian(n, m)]
    print("sp_pascal:")
    dump(pivs)

    k = n-m
    n, k = 2*n+1, k+1
    #print(n, k, "self-dual css")
    pivs = [tuple(piv) for piv, H in selfdual_css(n, m)]
    print("selfdual_css:")
    dump(pivs)


def get_dist(H):
    m, n = H.shape
    counts = [0 for i in range(n+1)]
    for v in span(H):
        counts[v.sum()] += 1
    return tuple(counts)

def get_dist_sy(H):
    m, nn = H.shape
    assert nn%2 == 0
    n = nn//2
    counts = [0 for i in range(n+1)]
    for v in span(H):
        d = 0
        for i in range(n):
            if v[i] or v[n+i]:
                d += 1
        counts[d] += 1
    return tuple(counts)
    

def show_dist(Hs, cb):
    items = [get_dist(H) for H in Hs]
    dist = {}
    for H in Hs:
        k = cb(H)
        dist[k] = dist.get(k, 0) + 1
    keys = list(dist.keys())
    #keys.sort(key = lambda k:dist[k])
    keys.sort() # ?
    print(len(keys))
    count = 0
    for k in keys:
        count += dist[k]
        print('\t', k, dist[k])
    print("\t =", count)
    

def test_Hs(n, m):
    k = n-m
    print("[[n=%d,m=%d,k=%d]]"%(n, m, k))
    Hs = [H for piv, H in sp_pascal.grassmannian(n, m)]
    show_dist(Hs, cb=get_dist_sy)

def test_Hs_selfdual(n, m):
    k = n - 2*m
    print("[[n=%d,m=%d,k=%d]]"%(n, m, k), "self-dual css")
    Hs = [H for piv, H in selfdual_css(n, k)]
    for H in Hs:
        assert H.shape == (m,n), (H.shape)
    show_dist(Hs, cb=get_dist)

def test():
    #test_pivs(3, 2)
    #test_pivs(4, 2)
    # hmm....

    test_Hs(3, 3)
    test_Hs_selfdual(7,3)
    print()

    test_Hs(3, 2)
    test_Hs_selfdual(7,2)
    print()

    test_Hs(3, 1)
    test_Hs_selfdual(7,1)
    print()

    test_Hs(4, 2)
    test_Hs_selfdual(9,2)
    print()

    test_Hs(4, 3)
    test_Hs_selfdual(9,3)
    print()


def omega(n, p=2):
    nn = 2*n
    F = zeros2(nn,nn)
    for i in range(n):
        F[i, n+i] = 1
        F[n+i, i] = p-1
    return F

def test_p3():
    p = 3
    for n in range(5):
     nn = 2*n
     F = omega(n, p)
     for m in range(n+1):
        count = 0
        for H in qchoose_2(nn, m, p):
            HFHt = dot(dot(H, F), H.transpose()) % p
            if HFHt.sum() == 0:
                count += 1
        print(count, " ", end="", flush=True)
     print()

def test_css_p3():
    p = 3
    for n in range(10):
      print("n=%s"%n, end=" ")
      for m in range(n):
        count = 0
        for H in qchoose_2(n, m, p):
            if (dot(H, H.transpose())%p).sum() != 0:
                continue
            count += 1
        print(count, " ", end="", flush=True)
      print()


def algebraic():
    # yes these agree...
    for n in range(1, 5):
     for m in range(1, n+1):
        get_dim_iso(n, m)
    print()
    for n in range(1, 5):
     for m in range(1, n+1):
        get_dim_biorth(n, m)

def get_dim_biorth(n, m):
    # biorthogonal 
    nn = 2*n + 1
    H = numpy.empty((m, nn), dtype=object)
    ring = Q
    for i in range(m):
     for j in range(nn):
        v = 'h(%d,%d)'%(i,j)
        v = Poly(v, ring)
        H[i,j] = v
    #print(H)
    lhs = dot(H, H.transpose())
    eqs = []
    for i in range(m):
     for j in range(m):
        eq = lhs[i, j]
        if eq==0:
            continue
        if eq in eqs:
            continue
        if -eq in eqs:
            continue
        #print(eq, "=", 0)
        eqs.append(eq)
    print("n=%d, m=%d"%(n,m), end=" ")
    print("vars:", m*nn, end=" ")
    print("eqs:", len(eqs), end=" ")
    print("dim =", m*nn - len(eqs))


def get_dim_iso(n, m):
    # Isotropic 
    nn = 2*n
    H = numpy.empty((m, nn), dtype=object)
    ring = Q
    for i in range(m):
     for j in range(nn):
        v = 'h(%d,%d)'%(i,j)
        v = Poly(v, ring)
        H[i,j] = v
    #print(H)
    F = numpy.zeros((nn, nn))
    for i in range(n):
        F[i, n+i] = 1
        F[n+i, i] = -1
    lhs = dot(dot(H, F), H.transpose())
    eqs = []
    for i in range(m):
     for j in range(m):
        eq = lhs[i, j]
        if eq==0:
            continue
        if eq in eqs:
            continue
        if -eq in eqs:
            continue
        #print(eq, "=", 0)
        eqs.append(eq)
    print("n=%d, m=%d"%(n,m), end=" ")
    print("vars:", m*nn, end=" ")
    print("eqs:", len(eqs), end=" ")
    print("dim =", m*nn - len(eqs))


if __name__ == "__main__":
    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))




