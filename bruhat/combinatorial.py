#!/usr/bin/env python


"""
Weyl groups...

"""


import sys, os
from time import time
start_time = time()
from functools import reduce
from operator import add, mul
from string import ascii_lowercase
letters = ascii_lowercase

from bruhat.action import Perm, Group, Coset, mulclose, Action
from bruhat.smap import SMap
from bruhat.argv import argv


def build_A(n):
    letters = ascii_lowercase[:n]
    items = list(letters)

    gen = []
    for i in range(n-1):
        perm = {c:c for c in items}
        perm[items[i]] = items[i+1]
        perm[items[i+1]] = items[i]
        gen.append(Perm(perm, items))
    G = Group.generate(gen, items=items)
    return G


def build_B(n):
    letters = ascii_lowercase[:n]
    items = ["+"+c for c in letters] + ["-"+c for c in letters]
    pairs = [("+"+c, "-"+c) for c in letters]

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
    if n:
        perm = {c:c for c in items}
        a, b = pairs[-1]
        perm[a] = b
        perm[b] = a
        gen.append(Perm(perm, items))
    G = Group.generate(gen, items=items)
    return G


def test_weyl():

    N = 6

    Gs = {} # Group's
    Hs = {} # stabilizer Coset's
    Xs = {} # Action's
    for n in range(N):
        G = build_A(n)
        Gs[n] = G
        items = G.items
        #print(len(G))
        for m in range(n+1):
            H = []
            fix = {items[j] for j in range(m)}
            for g in G:
                send = {g[x] for x in fix}
                if fix==send:
                    H.append(g)
            Hs[n, m] = H
            H = Coset(H, G.items)
            X = G.action_subgroup(H)
            X.basepoint = H
            Xs[n, m] = X

    def get_pos(n, m):
        K = 4
        return n, m*K + int(round((3-n/2)*K))

    smap = SMap()
    for n in range(N):
        for m in range(n+1):
            G = Gs[n]
            H = Hs[n, m]
            smap[get_pos(n, m)] = str(len(G)//len(H))
    print(smap)
    print()

    def check_group_hom(G0, G1, hom):
        assert len(G0) == len(hom)
        for g in G0:
          for h in G0:
            assert hom[g*h] == hom[g]*hom[h]
        assert len(set(hom.values())) == len(hom), "not injective"

    smap = SMap()
    lhoms = {}
    rhoms = {}
    for n in range(N-1):
        # construct two Group hom's, G0 --> G1
        G0 = Gs[n]
        G1 = Gs[n+1]
        items = G1.items

        # add new item on the right
        hom = {}
        for g in G0:
            perm = dict(g.perm)
            perm[items[-1]] = items[-1]
            h = Perm(perm, items)
            assert h in G1
            hom[g] = h
        #print("hom:", len(hom))
        check_group_hom(G0, G1, hom)
        rhoms[n] = hom

        # add new item on the left
        lookup = {src : letters[letters.index(src)+1] for src in G0.items}
        hom = {}
        for g in G0:
            perm = {lookup[src] : lookup[g.perm[src]] for src in G0.items}
            #print(g.perm, perm)
            perm[items[0]] = items[0]
            #print(g.perm, perm)
            h = Perm(perm, items)
            assert h in G1
            hom[g] = h
        check_group_hom(G0, G1, hom)
        #print("hom:", len(hom))
        lhoms[n] = hom

    lfuncs = {}
    rfuncs = {}
    for n in range(N-1):
      for m in range(n+1):
        print("\n%d choose %d" % (n, m))

        # construct a right map
        src = Xs[n, m]
        tgt = Xs[n+1, m+1]
        print(len(src.items), ">-->", len(tgt.items))
        func = {} # map src --> tgt
        for x in src.items:
            # find a coset representative
            for g in src.G:
                if src[g](src.basepoint) == x:
                    break
            else:
                assert 0
            h = lhoms[n][g]
            y = tgt[h](tgt.basepoint)
            func[x] = y
        assert len(func) == len(src.items)
        assert len(set(func.values())) == len(func)
        rfuncs[n, m] = func

        # construct a left map
        src = Xs[n, m]
        tgt = Xs[n+1, m]
        func = {} # map src --> tgt
        for x in src.items:
            # find a coset representative
            for g in src.G:
                if src[g](src.basepoint) == x:
                    break
            else:
                assert 0
            h = lhoms[n][g]
            y = tgt[h](tgt.basepoint)
            func[x] = y
        assert len(func) == len(src.items)
        assert len(set(func.values())) == len(func)
        lfuncs[n, m] = func

    n, m = 4, 2
    rmap = rfuncs[n-1, m-1]
    lmap = lfuncs[n-1, m]

    tgt = Xs[n, m]
    print(len(tgt.items))
    found = set(tgt.items)
    image = set(rmap.values())
    print(len(set(rmap.keys())), "-->", len(image))
    for item in rmap.values():
        print(item in found)
        #assert item in found
        #found.remove(item)
    for item in lmap.values():
        assert item in found
        found.remove(item)

    return

    def show(coset):
        items = []
        for g in coset.perms:
            keys = list(g.perm.keys())
            keys.sort()
            #s = " ".join("%s:%s"%(key,g[key]) for key in keys)
            s = " ".join("%s"%(g[key]) for key in keys)
            items.append("{%s}"%s)
        return "[%s]"%(", ".join(items))

    def dump(n, m):
        X = Xs[n, m]
        print("|%d choose %d| = %d"%(n, m, len(X.items)))
        for x in X.items:
            print('\t', show(x))
    dump(3, 2)
    dump(4, 2)


if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))






