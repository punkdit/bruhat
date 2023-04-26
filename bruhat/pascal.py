#!/usr/bin/env python
"""
build pascal triangles for various code families.
"""

from functools import cache

from bruhat.dev.geometry import all_codes
from bruhat.codes import is_morthogonal
from bruhat.util import choose
from bruhat.argv import argv

from time import time
start_time = time()



def is_morthogonal(G, m):
    k = len(G)
    if m==1:
        for v in G:
            if v.sum()%2 != 0:
                return False
        return True
    if m>1 and not is_morthogonal(G, m-1):
        return False
    items = list(range(k))
    for idxs in choose(items, m):
        v = G[idxs[0]]
        for idx in idxs[1:]:
            v = v * G[idx]
        if v.sum()%2 != 0:
            return False
    return True

def is_strict_morthogonal(G, level):
    assert level>=1
    if level==1:
        for v in G:
            if v.sum()%2 != 0:
                return False
        return True
    k = len(G)
    if level==2:
     for i0 in range(k):
      for i1 in range(i0+1, k):
        v = G[i0]*G[i1]
        if v.sum()%2 != 0:
            return False
     return True
    if level==3:
     for i0 in range(k):
      for i1 in range(i0+1, k):
       for i2 in range(i1+1, k):
        v = G[i0]*G[i1]*G[i2]
        if v.sum()%2 != 0:
            return False
     return True
    assert 0

if argv.get("jit", True):
    import numba
    is_strict_morthogonal = numba.njit(is_strict_morthogonal)


def qchoose_2(n, m):
    assert m<=n
    col = m
    row = n-m
    for A in geometry.get_cell(row, col, 2):
        yield A


if argv.get("cache", True):
    def get_cache(m, n):
        return list(all_codes(m,n))
    get_cache = cache(get_cache)

    def even_codes(m, n, level=0):
        if level == 0:
            return get_cache(m, n)
        items = []
        for H in even_codes(m, n, level-1):
            if is_strict_morthogonal(H, level):
                items.append(H)
        return items

else:
    get_cache = all_codes

    def even_codes(m, n, level=0):
        if level == 0:
            for H in get_cache(m, n):
                yield H
            return
        for H in even_codes(m, n, level-1):
            if is_strict_morthogonal(H, level):
                yield H


def main():
    N = argv.get("N", 10)

    for level in range(4):
      print("level =", level)
      for n in range(N):
        nn = n+1
        if level:
            nn = n or 1
        if level > 1:
            nn = max(1, n-1)
        for m in range(nn):
            count = 0
            for H in even_codes(m, n, level):
                count += 1
            print("%8s"%count, end=" ", flush=True)
        print()
      print()


if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()



