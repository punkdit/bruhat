#!/usr/bin/env python

"""
build dual containing codes (self-dual CSS) from hecke operators
"""


from time import time
start_time = time()
from random import choice
from functools import reduce
from operator import mul

import numpy

from bruhat.gset import Group, Perm
from bruhat.argv import argv
from bruhat.solve import zeros2, dot2, span, shortstr, linear_independent
from bruhat.util import choose
from bruhat.orthogonal import get_logops, row_reduce
from bruhat.todd_coxeter import Schreier



def make_hecke(X, Y):
    #print("make")
    G = X.src
    assert Y.src is G
    N = len(G)
    m, n = X.tgt.rank, Y.tgt.rank
    #print(X.send_perms)
    #print(X.tgt)
    remain = set((i,j) for i in range(m) for j in range(n))
    while remain:
        A = zeros2(m, n)
        ij = iter(remain).__next__()
        i, j = ij
        for idx in range(N):
            x = X.tgt[X.send_perms[idx]]
            y = Y.tgt[Y.send_perms[idx]]
            ij = x[i], y[j]
            if A[ij] == 0:
                A[ij] = 1
                remain.remove(ij)
        yield A


def get_lower_distance(H, L, max_d=4):
    m, n = H.shape
    idxs = list(range(n))
    for w in range(2, max_d):
      for ii in choose(idxs, w):
        v = zeros2(n, 1)
        for i in ii:
            v[i, 0] = 1
        if dot2(H, v).sum() == 0 and dot2(L, v).sum():
            return w
    return None



def get_distance_xz(H, L, min_d=2):

    m, n = H.shape
    if len(L) > 24:
        d = min(L.sum(1))
        return None, d

    if len(L) + len(H) > 24:
        d = n
        for u in span(L):
            d1 = u.sum()
            if d1 == 0:
                continue
            d = min(d, d1)
            if d<=min_d:
                return None, d
        return None, d

    vs = list(span(H))

    d = n
    for u in span(L):
        if u.sum() == 0:
            continue
        for v in vs:
            uv = (u+v)%2
            d = min(d, uv.sum())
        if d<=min_d:
            return d,d
    return d,d


def get_distance(H, L=None, min_d=2):
    # return lower & upper bound on distance
    m, n = H.shape
    if L is None:
        L = get_logops(H)

    max_d = 4 if n>100 else 5
    ld = get_lower_distance(H, L, max_d)
    if ld is not None:
        return ld, ld

    return get_distance_xz(H, L, min_d)


def get_distance_css(code):
    Hx, Hz, Lx, Lz = code.Hx, code.Hz, code.Lx, code.Lz
    dx = get_lower_distance(Hx, Lz)
    dz = get_lower_distance(Hz, Lx)
    if dx is None:
        _, dx = get_distance_xz(Hx, Lx)
    if dz is None:
        _, dz = get_distance_xz(Hz, Lz)

    return dx, dz

    
def get_group(name):
    if "*" in name:
        Gs = [get_group(name) for name in name.split("*")]
        G = reduce(mul, Gs)
        return G

    graph = None
    stem = name[0]
    n = int(name[1:]) if name[1:] else None
    if stem == "A":
        G = Group.alternating(n)
    elif stem == "S":
        G = Group.symmetric(n)
    elif stem == "I":
        G = Group.dihedral(n)
    elif stem == "C":
        G = Group.cyclic(n)
    elif stem == "K":
        ngens = 2
        a, ai, b, bi = range(2*ngens)
    
        # Klein Quartic
        rels = [ (ai, a), (bi, b), (a,)*3, (b,)*7, (a,b)*2 ]
        rels += [ (a,bi,bi)*4 ]
        graph = Schreier(2*ngens, rels)
        graph.build()
        assert len(graph) == 168, len(graph)

    elif stem == "D":
        graph = Schreier.make_D(n) # D4 has order 192

    elif stem == "B":
        #graph = Schreier.make_B(3) # nada
        graph = Schreier.make_B(n) # B4 has order 384

    elif stem == "Z":
        ngens = 3
        r, g, b = (0,), (1,), (2,)
        rels = [r*2, g*2, b*2, (r+b)*3, (r+g)*3, (b+g)*3]
        rels += [(b+r+g)*4]
        graph = Schreier(ngens, rels)
        graph.build(maxsize=100000)

    elif stem == "G":
        from bruhat.qcode import Geometry
        key = argv.get("key", (4, 4))
        idx = argv.get("idx", 24)
        geometry = Geometry(key, idx)
        graph = geometry.build_graph()

    else:
        assert 0, name

    if graph is not None:
        print("graph:", len(graph))
        G = graph.get_group()
        N = len(G)
        perms = []
        for g in G:
            perm = [g[i] for i in range(N)]
            perm = Perm(perm)
            perms.append(perm)
        G = Group(perms)

    return G


def main():

    n = argv.get("n", 5)
    name = argv.next()
    G = get_group(name)

    print("|G| =", len(G))
    Hs = list(G.subgroups())
    Hs = [H for H in Hs if len(H)<len(G)]

    print("proper subgroups:", len(Hs))
    Xs = [G.action_subgroup(H) for H in Hs]

    i = 0
    while i < len(Xs):
        j = i+1
        Xi = Xs[i]
        while j < len(Xs):
            Xj = Xs[j]
            if Xi.rank==Xj.rank and Xs[i].isomorphic(Xs[j]):
                Xs.pop(j)
                print(".", end="", flush=True)
            else:
                j += 1
        i += 1
    print()
    print("|Xs| =", len(Xs))

    print("make_hecke...")
    ops = {}
    N = len(Xs)
    for i in range(N):
      for j in range(N):
        ops[i,j] = []
    for i, X in enumerate(Xs):
      for j, Y in enumerate(Xs):
        for A in make_hecke(X, Y):
            A = linear_independent(A)
            m, n = A.shape
            if m<3 or n<5:
                continue
            ops[i,j].append(A)
    print("ops:", len(ops))

    min_d = argv.get("min_d", 3)

    if argv.css:
        print("CSSCode's")
        from qupy.ldpc.css import CSSCode
        for i in range(N):
         for j in range(N):
          for k in range(j, N):
            Hxs = ops[j,i]
            Hzs = ops[k,i]
            for Hx in Hxs:
             for Hz in Hzs:
                assert Hx.shape[1] == Hz.shape[1]
                mx, n = Hx.shape
                mz, n = Hz.shape
                k = n-mx-mz
                if k==0 or k>=n//2:
                    continue
                if Hx.tobytes() == Hz.tobytes(): # self-dual
                    continue
                if dot2(Hx, Hz.transpose()).sum() == 0:
                    code = CSSCode(Hx=Hx, Hz=Hz)
                    #dx, dz = code.distance(2)
                    dx, dz = get_distance_css(code)
                    if (dx is not None and dx<min_d) or (dz is not None and dz < min_d):
                        print(".", end="", flush=True)
                        continue
                    print()
                    print(type(dx) is int, dx<min_d)
                    print("Hx=")
                    print(shortstr(Hx))
                    print("Hz=")
                    print(shortstr(Hz))
                    print("[[%d, %d, %s %s]"%(n, k, dx, dz), min_d)
        return


    found = set()
    for key in ops:
      for A in ops[key]:
        #print(A, A.shape)
        H = row_reduce(A)
        m, n = H.shape
        k = n-2*m
        if k<=0:
            continue
        ws = H.sum(0)
        if ws.min() == 0 or ws.max() == 1:
            continue
        HHt = dot2(H, H.transpose())
        if HHt.max():
            continue

        dmin, dmax = get_distance(H)
        if dmax<min_d:
            print(choice("/\\"), end='', flush=True)
            continue
        w = A.sum(1)[0]
        desc = "[[%d, %d, %s<=d<=%s]]_%d"%(n, k, dmin, dmax, w)
        if desc in found:
            continue
        found.add(desc)
        print()
        print("A.shape =", A.shape)
        A = linear_independent(A)
        print(shortstr(A))
        print(desc)
        print()




if __name__ == "__main__":
    fn = "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds.\n"%(time() - start_time))



