#!/usr/bin/env python

"""
build dual containing codes (self-dual CSS) from hecke operators
"""


from time import time
start_time = time()
from functools import reduce
from operator import add

import numpy

from bruhat.gset import Group
from bruhat.argv import argv
from bruhat.solve import zeros2, row_reduce, dot2, span
from bruhat.orthogonal import get_logops



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
        return
    print("|G| =", len(G))
    Hs = list(G.subgroups())
    Hs = [H for H in Hs if len(H)>1 and len(H)<len(G)]
    print("proper subgroups:", len(Hs))
    Xs = [G.action_subgroup(H) for H in Hs]

    print("make...")
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

            L = get_logops(H)
            vs = list(span(H))
    
            d = n
            for u in span(L):
                if u.sum() == 0:
                    continue
                for v in vs:
                    uv = (u+v)%2
                    d = min(d, uv.sum())
            if d<3:
                continue

            print(H, H.shape, "[[%d, %d, %d]]"%(n, k, d))




if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))



