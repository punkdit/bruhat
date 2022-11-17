#!/usr/bin/env python


"""
Weyl groups...

"""


import sys, os
from time import time
start_time = time()
import random
from random import randint, choice
from functools import reduce
from operator import add, mul
from string import ascii_lowercase

from bruhat.action import Perm, Group, Coset, mulclose
from bruhat.action import Group
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

    Hs = {}
    Gs = {}
    for n in [0, 1, 2, 3, 4]:
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
            print(len(G) // len(H), end=" ")
        print()

    print()

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
        G = Gs[n]
        H = Hs[n, m]
        X = G.left_cosets(H)
        print("|%d choose %d| = %d"%(n, m, len(X)))
        for x in X:
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






