#!/usr/bin/env python3

"""
see also: biset.py, combinatorial.py

"""

from time import time
start_time = time()

from bruhat.action import Group, Perm, Morphism, Action, Hom
from bruhat.equ import Equ
from bruhat.argv import argv


def induced_action(G, H, X):
    # construct induced GSet action as G*X _modulo an equivalence relation
    assert X.G is H
    GX = [(g, x) for g in G for x in X.items]
    lookup = {gx:Equ(gx) for gx in GX}
    for (g,x) in GX:
      for h in H:
        lhs = lookup[g*h, x] # XXX should this be ~h ???
        rhs = lookup[g, X(h)[x]]
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
    action = Action(G, send_perms, items, check=True)
    return action


def coinduced_action(G, H, X, check=True):
    src = G.cayley_action(H)
    tgt = X
    items = [f for f in src.get_homs(tgt)]

    found = set(items)
    assert len(found) == len(items)
    #print("homs:", len(items))

    # now we find an action of G on items
    send_perms = {}
    for g in G:
        perm = {}
        for hom in items:
            # hom : G --> X
            send_items = {h:hom.send_items[h*g] for h in src}
            gom = Hom(src, tgt, send_items)
            assert gom in found
            perm[hom] = gom
        perm = Perm(perm, items)
        send_perms[g] = perm
    action = Action(G, send_perms, items, check=check)
    return action


def restricted_action(morphism, action, check=True):
    src = morphism.src
    tgt = morphism.tgt
    assert action.G is tgt
    send_perms = {g:action(morphism(g)) for g in src}
    action = Action(src, send_perms, action.items, check=check)
    return action


def test_restriction():
    n = argv.get("n", 3)
    G = Group.symmetric(n)
    #H = Group([g for g in G if g[n-1] == n-1], G.items)
    H = Group([g for g in G if g.sign()==1], G.items)

    morphism = Morphism(H, G, {g:g for g in H})
    
    X = G.tautological_action()
    Y = restricted_action(morphism, X)
    print(Y)


def test_coinduction_products():

    n = argv.get("n", 3)
    G = Group.symmetric(n)
    #H = Group([g for g in G if g[n-1] == n-1], G.items)
    H = Group([g for g in G if g.sign()==1], G.items)

    # right adjoints preserve products
    Xs = [H.action_subgroup(K) for K in H.subgroups()]
    for X1 in Xs:
      Y1 = coinduced_action(G, H, X1)
      for X2 in Xs:
        Y2 = coinduced_action(G, H, X2)
        X = X1 * X2
        Y = coinduced_action(G, H, X)
        #print(len(Y1), len(Y2), len(Y))
        assert len(Y1) * len(Y2) == len(Y)


def test_induction():

    G = Group.symmetric(4)

    H = Group([g for g in G if g[3] == 3], G.items)

    H1 = Group.symmetric(3)
    items = H1.items
    send_perms = {}
    for h0 in H:
      for h1 in H1:
        if [h0[i] for i in items] == [h1[i] for i in items]:
            send_perms[h0] = h1
    X = Action(H, send_perms, items, check=True)
    action = induced_action(G, H, X)
    assert len(action) == 12

    print(len(list(action.get_homs(G.tautological_action()))))

    return

    # left adjoints preserve coproducts
    Xs = [H.action_subgroup(K) for K in H.subgroups()]
    for X1 in Xs:
      Y1 = induced_action(G, H, X1)
      for X2 in Xs:
        Y2 = induced_action(G, H, X2)
        X = X1 + X2
        Y = induced_action(G, H, X)
        assert len(Y1) + len(Y2) == len(Y)


def test_coinduction_general():
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
    X = Action(H, send_perms, items, check=True)

    action = coinduced_action(G, H, X)
    print(action)


def test_hom():

    C3 = Group.cyclic(3)
    X = C3.tautological_action()
    assert len(list(X._get_homs_atomic(X))) == 3

    G = Group.symmetric(3)
    Xs = [G.action_subgroup(H) for H in G.subgroups()]
    for X in Xs:
      for Y in Xs:
        homs = list(X._get_homs_atomic(Y))
        #print("%2d"%len(homs), end=" ")
      #print()

    G = Group.symmetric(4)
    H = Group([g for g in G if g[3] == 3], G.items)

    A = G.tautological_action()
    B = H.tautological_action()
    assert len(B.orbits()) == 2

    G = Group.cyclic(2)
    X = G.tautological_action()
    X = X.coproduct(X)
    assert len(list(X.get_homs(X))) == 16
    #for hom in X.get_homs(X):
    #    print(hom)


if __name__ == "__main__":
    fn = argv.next() or "test"

    print("%s()"%fn)

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))




