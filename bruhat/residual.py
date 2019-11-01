#!/usr/bin/env python3

from bruhat.argv import argv
from bruhat.weyl import Weyl
from bruhat.action import Group, Coset, Perm, mulclose, burnside
from bruhat.util import all_subsets
from bruhat.poly import Poly
from bruhat import element



def show_qpoly(G):
    gen = G.gen
    roots = G.roots
    els = G.generate()
    e = G.identity
    G = Group(els, roots)
    print("order:", len(els))

    ring = element.Z
    value = zero = Poly({}, ring)
    q = Poly("q", ring)
    for g in els:
        #print(g.word)
        value = value + q**(len(g.word))
    print(value.qstr())

    #lookup = dict((g, g) for g in G) # remember canonical word

    n = len(gen)
    for i in range(n):
        gen1 = gen[:i] + gen[i+1:]
        H = mulclose_short([e]+gen1)
        eH = Coset(H, roots)
        H = Group(H, roots)
        #gHs = G.left_cosets(H)

        cosets = set([eH])
        bdy = set(cosets)
        while bdy:
            _bdy = set()
            for coset in bdy:
              for g in gen:
                gH = Coset([g*h for h in coset], roots)
                if gH not in cosets:
                    cosets.add(gH)
                    _bdy.add(gH)
            bdy = _bdy

        value = zero
        for gH in cosets:
            items = list(gH)
            items.sort(key = lambda g : len(g.word))
            #for g in items:
            #    print(g.word, end=" ")
            #print()
            g = items[0]
            value = value + q**len(g.word)
        #print(len(gH))
        print(value.qstr())




def main():

    desc = argv.next()
    assert desc

    print("constructing %s"%desc)
    assert "_" in desc, repr(desc)

    name, idx = desc.split("_")
    idx = int(idx)
    attr = getattr(Weyl, "build_%s"%name)
    G = attr(idx)

    print(G)

    e = G.identity
    gen = G.gen
    roots = G.roots
    els = G.generate()
    G = Group(els, roots)
    print("order:", len(els))

    ring = element.Z
    value = zero = Poly({}, ring)
    q = Poly("q", ring)
    for g in els:
        #print(g.word)
        value = value + q**(len(g.word))
    print(value.qstr())

    n = len(gen)
    Hs = []
    for idxs in all_subsets(n):
        print(idxs, end=" ")
        gen1 = [gen[i] for i in idxs] or [e]
        H = Group(mulclose(gen1), roots)
        Hs.append(H)
        gHs = G.left_cosets(H)
        value = zero
        for gH in gHs:
            items = list(gH)
            items.sort(key = lambda g : len(g.word))
            #for g in items:
            #    print(g.word, end=" ")
            #print()
            g = items[0]
            value = value + q**len(g.word)
        #print(len(gH))
        print(value.qstr())

    G = Group(els, roots)
    burnside(G, Hs)




if __name__ == "__main__":
    main()



