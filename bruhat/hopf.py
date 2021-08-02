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


def main():

    ring = element.Q
    zero = ring.zero
    one = ring.one

    n = argv.get("n", 3)
    G = Group.symmetric(n)

    print(G)

    d = len(G)
    K = Space(ring, 1, name="K")
    V = Space(ring, d, name="V")
    VV = V@V

    scalar = K.identity()
    I = V.identity()
    swap = VV.get_swap()
    
    # green spiders
    g_ = Lin(K, V) # uniform discard
    _g = Lin(V, K) # uniform create
    g_gg = Lin(VV, V) # copy
    gg_g = Lin(V, VV) # pointwise mul

    for i in range(d):
        g_[0, i] = one
        _g[i, 0] = one
        g_gg[i + d*i, i] = one
        gg_g[i, i + d*i] = one

    eq = lambda lhs, rhs : lhs.weak_eq(rhs)

    assert eq(g_gg >> (g_ @ I), I) # unit
    assert eq(g_gg >> (I @ g_), I) # unit
    assert eq(g_gg >> (g_gg @ I), g_gg >> (I@g_gg)) # assoc

    assert eq(gg_g * (_g @ I), I) # unit
    assert eq(gg_g * (I @ _g), I) # unit
    assert eq(gg_g * (gg_g @ I), gg_g * (I@gg_g)) # assoc

    assert eq((g_gg @ I) >> (I @ gg_g), (I @ g_gg) >> (gg_g @ I)) # frobenius
    assert eq((g_gg @ I) >> (I @ gg_g), gg_g >> g_gg) # extended frobenius

    assert eq(_g >> g_, d*scalar)

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

    # K[G] is a bialgebra
    assert eq( rr_r >> g_, g_ @ g_)
    assert eq( _r >> g_gg, _r @ _r )
    assert eq( _r >> g_, scalar )
    assert eq( rr_r >> g_gg , (g_gg @ g_gg) >> (I @ swap @ I) >> (rr_r @ rr_r) )

    # K[G] is hopf
    rhs = g_ >> _r
    assert eq( g_gg >> (I @ inv) >> rr_r, rhs)
    assert eq( g_gg >> (inv @ I) >> rr_r, rhs)

    # k^G is a bialgebra
    assert eq( gg_g >> r_, r_ @ r_)
    assert eq( _g >> r_rr, _g @ _g )
    assert eq( _g >> r_, scalar )
    assert eq( gg_g >> r_rr , (r_rr @ r_rr) >> (I @ swap @ I) >> (gg_g @ gg_g) )

    # k^G is hopf
    rhs = r_ >> _g
    assert eq( r_rr >> (I @ inv) >> gg_g, rhs)
    assert eq( r_rr >> (inv @ I) >> gg_g, rhs)

    


if __name__ == "__main__":

    main()



