#!/usr/bin/env python


"""
Weyl groups...

"""


import sys, os
from time import time
start_time = time()
from functools import reduce, lru_cache
from operator import add, mul
from string import ascii_lowercase
letters = ascii_lowercase

from bruhat.action import Perm, Group, Coset, mulclose, Action
from bruhat.smap import SMap
from bruhat.argv import argv

cache = lru_cache(maxsize=None)

def show(coset):
    items = []
    for g in coset.perms:
        keys = list(g.perm.keys())
        keys.sort()
        #s = " ".join("%s:%s"%(key,g[key]) for key in keys)
        s = " ".join("%s"%(g[key]) for key in keys)
        items.append("{%s}"%s)
    return "[%s]"%(", ".join(items))

def check_group_hom(G0, G1, hom):
    assert len(G0) == len(hom)
    for g in G0:
      for h in G0:
        assert hom[g*h] == hom[g]*hom[h]
    assert len(set(hom.values())) == len(hom), "not injective"


def check_injection(src, tgt, func):
    assert set(func.keys()) == set(src)
    output = set(func.values())
    assert len(output) == len(func)
    assert output.issubset(tgt)



class Pascal(object):
    "Pascal's triangle for thin geometry Weyl group"

    def __init__(self):
        pass

    def str(self, N):
        def get_pos(n, m):
            K = 4
            return n, m*K + int(round((3-n/2)*K))
    
        smap = SMap()
        for n in range(N):
            for m in range(n+1):
                G = Gs[n]
                H = Hs[n, m]
                smap[get_pos(n, m)] = str(len(G)//len(H))
        return str(smap)

    def get_G(self, n):
        assert 0, "abstract base class"

    @cache
    def get_X(self, n, m):
        # fix 0,....,m-1
        G = self.get_G(n)
        items = G.items
        assert 0<=m<=n
        H = []
        fix = {items[j] for j in range(m)}
        for g in G:
            send = {g[x] for x in fix}
            if fix==send:
                H.append(g)
        #Hs[n, m] = H
        H = Coset(H, G.items)
        X = G.action_subgroup(H)
        assert X.basepoint == H
        #Xs[n, m] = X
        return X

    @cache
    def get_Y(self, n, m):
        # fix 1,...,m
        G = self.get_G(n)
        items = G.items
        assert 0<=m<n
        J = []
        fix = {items[j] for j in range(1,m+1)}
        for g in G:
            send = {g[x] for x in fix}
            if fix==send:
                J.append(g)
        #Js[n, m] = J
        J = Coset(J, G.items)
        Y = G.action_subgroup(J)
        assert Y.basepoint == J
        #Ys[n, m] = Y
        return Y

    @cache
    def get_iso(self, n, m):
        src = self.get_Y(n, m)
        tgt = self.get_X(n, m)
        #assert src.isomorphic(tgt) # yes

        G = self.get_G(n)
        items = G.items
        assert len(items) == n
        k = {items[j]:items[j-1] for j in range(1,n)}
        k[items[0]] = items[-1]
        k = Perm(k, items)
        assert k in G

        # Construct an isomorphism of G-sets X-->Y
        send = {}
        for y in src.items:
            g = src.repr[y]
            x = tgt[g*~k](tgt.basepoint)
            assert x in tgt.items
            send[y] = x
        src.check_isomorphism(tgt, send)
        #isos[n, m] = send

        # this also works
        #iso = iter(src.isomorphisms(tgt)).__next__()
        #assert iso is not None
        #isos[n, m] = iso

        return send

    @cache
    def get_rhom(self, n):
        # construct Group hom, G0 --> G1
        G0 = self.get_G(n)
        G1 = self.get_G(n+1)
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
        return hom

    @cache
    def get_lhom(self, n):
        # construct Group hom, G0 --> G1
        G0 = self.get_G(n)
        G1 = self.get_G(n+1)
        items = G1.items

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
        #lhoms[n] = hom
        return hom

    @cache
    def get_rfunc(self, n, m):
        assert 0<=n
        assert 0<=m<=n
        #print("\n%d choose %d" % (n, m))

        # construct a right map
        src = self.get_X(n, m)
        tgt = self.get_X(n+1, m+1)
        #print(len(src.items), ">-->", len(tgt.items))
        func = {} # map src --> tgt
        for x in src.items:
            # find a coset representative
            g = src.repr[x]
            assert src[g](src.basepoint) == x
            h = self.get_lhom(n)[g]
            y = tgt[h](tgt.basepoint)
            func[x] = y
        check_injection(src.items, tgt.items, func)
        #rfuncs[n, m] = func
        return func

    @cache
    def get_lfunc(self, n, m):
        assert 0<=n
        assert 0<=m<=n
        # construct a left map
        src = self.get_X(n, m)
        Y = self.get_Y(n+1, m)
        tgt = self.get_X(n+1, m)
        func = {} # map src --> tgt
        for x in src.items:
            # find a coset representative
            g = src.repr[x]
            assert src[g](src.basepoint) == x
            h = self.get_lhom(n)[g]
            y = Y[h](Y.basepoint)
            y = self.get_iso(n+1,m)[y] # send Y --> tgt
            func[x] = y
        check_injection(src.items, tgt.items, func)
        #lfuncs[n, m] = func
        return func


class PascalA(Pascal):

    @cache
    def get_G(self, n):
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
    
    
class PascalB(Pascal):

    @cache
    def get_G(self, n):
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

    print("test_weyl:")

    triangle = PascalA()
    assert triangle.get_G(3) is triangle.get_G(3)
    assert triangle.get_G(3) is not PascalA().get_G(3)

    N = 6
    n, m = 4, 2
    for n in range(2, N):
      for m in range(1, n):
        rfunc = triangle.get_rfunc(n-1, m-1)
        lfunc = triangle.get_lfunc(n-1, m)
    
        tgt = triangle.get_X(n, m)
        found = set(tgt.items)
        for item in rfunc.values():
            assert item in found
            found.remove(item)
    
        for item in lfunc.values():
            assert item in found
            found.remove(item)
        assert not found
        print("(%d %d)"%(n,m), end=" ")
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






