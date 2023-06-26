#!/usr/bin/env python3

from time import time
start_time = time()

from bruhat.action import Group, Perm, Action
from bruhat.equ import Equ
from bruhat.argv import argv


#def test_0():
#
#    G = Group.symmetric(4)
#    print(G)
#
#    H = Group([g for g in G if g[3] == 3], G.items)
#    print(H)
#
#    #X = H.tautological_action()
#
#    H1 = Group.symmetric(3)
#    items = H1.items
#    send_perms = {}
#    for h0 in H:
#      for h1 in H1:
#        if [h0[i] for i in items] == [h1[i] for i in items]:
#            send_perms[h0] = h1
#    X = Action(H, send_perms, items)
#    print("X:", X)
#
#    Y = G.action_subgroup(H)
#    print("Y:", Y)
#
#    #cosets = G.left_cosets(H)
#    cosets = Y.items
#    assert Y.basepoint in cosets
#    assert H is Y.basepoint
#
#    lookup = {H : G.identity}
#    for g in G:
#        H1 = Y.send_perms[g][H]
#        lookup[H1] = g
#    assert len(lookup) == len(cosets)
#
#    items = [(a, b) for a in Y.items for b in X.items]
#    print(len(items))
#    send_perms = {}
#    for g in G:
#        f = Y.send_perms[g]
#        H1 = f[H]
#        g1 = lookup[H1]
#        #assert f == Y.send_perms[g1]
#        f1 = Y.send_perms[g1]
#        h = (~g1)*g
#        assert h in H
#        perm = {}
#        for (a, b) in items:
#            perm[a, b] = (f1[a], h[b])
#        perm = Perm(perm, items)
#        send_perms[g] = perm
#    induced = Action(G, send_perms, items, check=True) # FAIL


def test():

    G = Group.symmetric(4)
    print(G)

    H = Group([g for g in G if g[3] == 3], G.items)
    print(H)

    H1 = Group.symmetric(3)
    items = H1.items
    send_perms = {}
    for h0 in H:
      for h1 in H1:
        if [h0[i] for i in items] == [h1[i] for i in items]:
            send_perms[h0] = h1
    X = Action(H, send_perms, items)

    # construct induced GSet action as G*X modulo an equivalence relation
    GX = [(g, x) for g in G for x in X.items]
    lookup = {gx:Equ(gx) for gx in GX}
    for (g,x) in GX:
      for h in H:
        lhs = lookup[g*h, x]
        rhs = lookup[g, h[x]]
        lhs.merge(rhs)

    equs = set(equ.top for equ in lookup.values())
    equs = list(equs)
    items = [tuple(equ.items) for equ in equs]
    send_perms = {}
    for g in G:
        perm = {}
        for src in items:
            g1, i = src[0]
            g1 = g*g1
            for tgt in items:
                if (g1, i) in tgt:
                    break
            else:
                assert 0
            perm[src] = tgt
        perm = Perm(perm, items)
        send_perms[g] = perm
    induced = Action(G, send_perms, items, check=True)
    print(induced)

    # co-induction...
    src = G.cayley_action(H)
    tgt = X
    count = 0
    for f in all_homs(src, tgt):
        count += 1
    print("homs:", count)


def all_homs_atomic(src, tgt):
    assert len(src.orbits()) == 1
    assert len(tgt.orbits()) == 1


def all_homs(src, tgt):
    assert src.G is tgt.G
    G = src.G
    orbits = src.orbits()
    print("orbits:", orbits, len(orbits))
    #yield None

    remain = list(src.items)
    send = {}
    while remain:
        item = remain.pop()


def test_hom():

    G = Group.symmetric(4)
    H = Group([g for g in G if g[3] == 3], G.items)

    A = G.tautological_action()
    B = H.tautological_action()
    print(B.orbits())


if __name__ == "__main__":
    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))




