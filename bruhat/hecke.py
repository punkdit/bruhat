#!/usr/bin/env python

"""
build dual containing codes (self-dual CSS) from hecke operators
"""


from time import time
start_time = time()
from random import choice

import numpy

from bruhat.gset import Group, Perm
from bruhat.argv import argv
from bruhat.solve import zeros2, dot2, span, shortstr
from bruhat.orthogonal import get_logops, row_reduce
from bruhat.todd_coxeter import Schreier



def make(G, X, Y):
    #print("make")
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


def get_distance(H, min_d=2):
    # return lower, upper bound on distance
    m, n = H.shape
    L = get_logops(H)

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

def get_group():

    if argv.klein:
        ngens = 2
        a, ai, b, bi = range(2*ngens)
    
        # Klein Quartic
        rels = [ (ai, a), (bi, b), (a,)*3, (b,)*7, (a,b)*2 ]
        rels += [ (a,bi,bi)*4 ]
        graph = Schreier(2*ngens, rels)
        graph.build()
        assert len(graph) == 168, len(graph)

    elif argv.D4:
        #graph = Schreier.make_B(3) # nada
        graph = Schreier.make_D(4) # order 192

    elif argv.B4:
        graph = Schreier.make_B(4) # order 384

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
    if argv.alternating:
        G = Group.alternating(n)
    elif argv.symmetric:
        G = Group.symmetric(n)
    elif argv.dihedral:
        G = Group.dihedral(n)
    elif argv.cyclic:
        G = Group.cyclic(n)
    else:
        G = get_group()
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

    print("make...")
    found = set()
    for X in Xs:
      for Y in Xs:
        for A in make(G, X, Y):
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
            if dmax<3:
                print(choice("/\\"), end='', flush=True)
                continue
            desc = "[[%d, %d, %s<=d<=%s]]"%(n, k, dmin, dmax)
            if desc in found:
                continue
            found.add(desc)
            print()
            print(shortstr(H))
            print(H.shape, desc)
            print()




if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))



