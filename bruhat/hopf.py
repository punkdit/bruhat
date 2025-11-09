#!/usr/bin/env python3

import numpy

from bruhat.argv import argv
from bruhat import elim
from bruhat.action import Perm, Group, mulclose, mulclose_hom
if argv.fast:
    from bruhat import _element as element
else:
    from bruhat import element
from bruhat.chain import Space, Lin


def test_hom():

    ring = element.Q
    n = argv.get("n", 3)

    V = Space(ring, n)

    # build action of symmetric group on the space V
    items = list(range(n))
    gen1 = []
    gen2 = []
    for i in range(n-1):
        perm = dict((item, item) for item in items)
        perm[items[i]] = items[i+1]
        perm[items[i+1]] = items[i]
        g = Perm(perm, items)
        gen1.append(g)
        A = elim.zeros(ring, n, n)
        for k,v in perm.items():
            A[v,k] = ring.one
        lin = Lin(V, V, A)
        gen2.append(lin)

    perms = mulclose(gen1)
    G = Group(perms, items)

    #print(G)

    action = mulclose_hom(gen1, gen2)
    for g in G:
      for h in G:
        assert action[g*h] == action[g]*action[h] # check it's a group hom


def main(n=3):

    ring = element.Q
    zero = ring.zero
    one = ring.one

    n = argv.get("n", n)
    if argv.cyclic:
        G = Group.cyclic(n)
    else:
        G = Group.symmetric(n)

    comm = G.is_abelian()

    print(G)

    d = len(G)
    K = Space(ring, 1, name="K")
    V = Space(ring, d, name="V")
    VV = V@V

    scalar = K.identity()
    I = V.identity()
    swap = VV.get_swap()

    lunit = Lin(V, K@V, elim.identity(ring, d))
    runit = Lin(V, V@K, elim.identity(ring, d))

    cap = Lin(K, V@V) # tgt, src
    cup = Lin(V@V, K) # tgt, src
    for i in range(d):
        cup[i + d*i, 0] = one
        cap[0, i + d*i] = one
    
    # green spiders
    g_ = Lin(K, V) # uniform discard
    _g = Lin(V, K) # uniform create
    g_gg = Lin(VV, V) # copy (comul)
    gg_g = Lin(V, VV) # pointwise mul

    for i in range(d):
        g_[0, i] = one
        _g[i, 0] = one
        g_gg[i + d*i, i] = one
        gg_g[i, i + d*i] = one

    eq = lambda lhs, rhs : lhs.weak_eq(rhs)

    assert eq(g_gg >> (g_ @ I), I) # counit
    assert eq(g_gg >> (I @ g_), I) # counit
    assert eq(g_gg >> (g_gg @ I), g_gg >> (I@g_gg)) # coassoc

    assert eq(gg_g * (_g @ I), I) # unit
    assert eq(gg_g * (I @ _g), I) # unit
    assert eq(gg_g * (gg_g @ I), gg_g * (I@gg_g)) # assoc

    assert eq((g_gg @ I) >> (I @ gg_g), (I @ g_gg) >> (gg_g @ I)) # frobenius
    assert eq((g_gg @ I) >> (I @ gg_g), gg_g >> g_gg) # extended frobenius

    assert eq(_g >> g_, d*scalar)

    assert eq(gg_g >> g_, cap)
    assert eq(_g >> g_gg, cup)

    # red spiders
    r_ = Lin(K, V)    # discard unit
    _r = Lin(V, K)    # create unit
    r_rr = Lin(VV, V) # comul
    rr_r = Lin(V, VV) # mul

    # hopf involution
    inv = Lin(V, V)

    lookup = dict((v,k) for (k,v) in enumerate(G))
    for i in range(d):
        g = G[i]
        if g.is_identity():
            r_[0, i] = one
            _r[i, 0] = one
        inv[lookup[~g], i] = one
            
        for j in range(d):
            h = G[j]
            gh = g*h
            k = lookup[gh]
            rr_r[k, i+j*d] = one
            r_rr[i+j*d, k] = one

    assert eq(r_rr >> (r_ @ I), I) # unit
    assert eq(r_rr >> (I @ r_), I) # unit
    assert eq(r_rr >> (r_rr @ I), r_rr >> (I@r_rr)) # assoc

    assert eq(rr_r * (_r @ I), I) # unit
    assert eq(rr_r * (I @ _r), I) # unit
    assert eq(rr_r * (rr_r @ I), rr_r * (I@rr_r)) # assoc

    assert eq((r_rr @ I) >> (I @ rr_r), (I @ r_rr) >> (rr_r @ I)) # frobenius
    assert eq((r_rr @ I) >> (I @ rr_r), rr_r >> r_rr) # extended frobenius

    assert eq((_r >> r_) , scalar )

    print("gg_g")
    print(gg_g)
    print("g_gg")
    print(g_gg)
    print("rr_r")
    print(rr_r)
    print(r_)

    if n>2:
        assert not eq(rr_r >> r_, cap)
        assert not eq(_r >> r_rr, cup)
    else:
        print( eq(rr_r >> r_, cap) )
        print( eq(_r >> r_rr, cup) )

    # K[G] is a bialgebra
    assert eq( rr_r >> g_, g_ @ g_)
    assert eq( _r >> g_gg, _r @ _r )
    assert eq( _r >> g_, scalar )
    if not argv.skip:
        assert eq( rr_r >> g_gg , (g_gg @ g_gg) >> (I @ swap @ I) >> (rr_r @ rr_r) )
    print("K[G] is comm  ", eq(swap >> rr_r, rr_r))
    print("K[G] is cocomm", eq(g_gg >> swap, g_gg))

    # K[G] is hopf
    rhs = g_ >> _r
    assert eq( g_gg >> (I @ inv) >> rr_r, rhs)
    assert eq( g_gg >> (inv @ I) >> rr_r, rhs)

    # k^G is a bialgebra
    assert eq( gg_g >> r_, r_ @ r_)
    assert eq( _g >> r_rr, _g @ _g )
    assert eq( _g >> r_, scalar )
    if not argv.skip:
        assert eq( gg_g >> r_rr , (r_rr @ r_rr) >> (I @ swap @ I) >> (gg_g @ gg_g) )

    # k^G is hopf
    rhs = r_ >> _g
    assert eq( r_rr >> (I @ inv) >> gg_g, rhs)
    assert eq( r_rr >> (inv @ I) >> gg_g, rhs)
    print("k^G is comm   ", eq(swap >> gg_g , gg_g))
    print("k^G is cocomm ", eq(r_rr >> swap, r_rr))

    #print(rr_r)
    #print(r_rr)

    # unimodular
    r_cup = _r >> r_rr
    g_cap = gg_g >> g_
    assert eq(r_cup >> (I @ g_), _g)
    assert eq(r_cup >> (g_ @ I), _g)
    assert eq((I @ _r) >> g_cap, r_)
    assert eq((_r @ I) >> g_cap, r_)
    assert eq(inv, (I @ r_cup) >> (swap @ I) >> (I @ g_cap))
    assert eq(inv, (r_cup @ I) >> (I @ swap) >> (g_cap @ I))

    assert eq(r_rr >> rr_r, d*I)
    assert eq(g_gg >> gg_g, I) # special

    # complementary frobenius structures ?
    # Heunen & Vicary eq (6.4)
    lhs = (_r @ I) >> (r_rr @ g_gg) >> (I @ gg_g @ I) >> (I @ g_ @ I) >> (I @ lunit) >> rr_r
    rhs = g_ >> _r
    assert eq(lhs, rhs)

    lhs = (I @ _r) >> (r_rr @ g_gg) >> (I @ gg_g @ I) >> (I @ g_ @ I) >> (I @ lunit) >> rr_r
    #assert eq(lhs, rhs) # FAIL

    lhs = (_g @ I) >> (g_gg @ r_rr) >> (I @ rr_r @ I) >> (I @ r_ @ I) >> (I @ lunit) >> gg_g
    rhs = r_ >> _g
    assert eq(lhs, rhs)

    lhs = (I @ _g) >> (g_gg @ r_rr) >> (I @ rr_r @ I) >> (I @ r_ @ I) >> (I @ lunit) >> gg_g
    #assert eq(lhs, rhs) # FAIL

    # Heunen & Vicary eq (6.5)
    lhs = (_r @ I) >> (r_rr @ I) >> (I @ gg_g) >> (I @ g_)
    rhs = (I @ _r) >> (I @ r_rr) >> (gg_g @ I) >> (g_ @ I)
    assert eq(lhs, rhs)

    lhs = (_g @ I) >> (g_gg @ I) >> (I @ rr_r) >> (I @ r_)
    rhs = (I @ _g) >> (I @ g_gg) >> (rr_r @ I) >> (r_ @ I)
    assert eq(lhs, rhs)

    assert eq( r_rr, r_rr >> swap ) == G.is_abelian()
    assert eq( rr_r, swap >> rr_r ) == G.is_abelian()

    assert eq( _r >> r_rr, _r >> r_rr >> swap )
    assert eq( rr_r >> r_, swap >> rr_r >> r_ )

    #print(r_rr >> gg_g)
    #print(r_ >> _g)
    #print(g_gg >> rr_r)
    #print(rr_r >> r_)
    #print(_r >> r_rr)

    class NS:
        pass
    ns = NS()
    ns.__dict__.update(locals())
    return ns


if __name__ == "__main__":

    main()



