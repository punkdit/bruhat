#!/usr/bin/env python


"""
Weyl groups (thin geometry) and their Pascal triangles,
and bijective (combinatorial) correspondences.

"""


import sys, os
from time import time
start_time = time()
import random
from functools import reduce, lru_cache
from operator import add, mul
from string import ascii_lowercase

from bruhat.action import Perm, Group, Coset, mulclose, mulclose_hom, Action
from bruhat.algebraic import Matrix, Algebraic
from bruhat.poset import Poset
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


class Space(Action):
    "Homogeneous Space."
    "Make a choice as a subset that is fix'ed by a Group action."
    def __init__(self, G, fix):
        items = G.items
        fix = set(fix)
        assert fix.issubset(items)
        H = []
        for g in G:
            send = {g[x] for x in fix}
            if fix==send:
                H.append(g)
        H = Coset(H, G.items)
        X = G.action_subgroup(H)
        assert X.basepoint == H
        self.__dict__.update(X.__dict__) # doh..
        self.H = H
        self.fix = fix
        #self.remain = set(G.items).difference(fix)

    def get_iso(src, tgt):
        # Construct an isomorphism of G-sets
        assert isinstance(tgt, Space)
        assert tgt.G is src.G
        assert len(tgt.fix) == len(src.fix)
        G = src.G

        for g in G:
            for x in src.fix:
                if g(x) not in tgt.fix:
                    break
            else:
                break
        else:
            assert 0

        send = {}
        for y in src.items:
            h = src.repr[y]
            x = tgt(h*~g)(tgt.basepoint)
            assert x in tgt.items
            send[y] = x
        src.check_isomorphism(tgt, send)

        # this also works
        #iso = iter(src.isomorphisms(tgt)).__next__()
        #assert iso is not None

        return send



class Pascal(object):
    "Pascal's triangle for thin geometry Weyl group"

    def __init__(self):
        self.choices = {} # map fix -> Space

    def str(self, N):
        def get_pos(n, m):
            K = 4
            return n, m*K + int(round((3-n/2)*K))
    
        smap = SMap()
        for n in range(N):
            for m in range(n+1):
                X = self.get_X(n, m)
                smap[get_pos(n, m)] = str(len(X.items))
        return str(smap)

    def get_G(self, n):
        assert 0, "abstract base class"

    def get_fix(self, G, fix):
        key = list(fix)
        key.sort()
        key = G, tuple(key)
        if key in self.choices:
            return self.choices[key]
        choice = Space(G, fix)
        self.choices[key] = choice
        return choice

    @cache
    def get_X(self, n, m):
        # fix 0,....,m-1
        assert 0<=m<=n
        G = self.get_G(n)
        fix = {G.items[j] for j in range(m)}
        return self.get_fix(G, fix)

    @cache
    def get_Y(self, n, m):
        # fix 1,...,m
        assert 0<=m<n
        G = self.get_G(n)
        fix = {G.items[j] for j in range(1,m+1)}
        return self.get_fix(G, fix)

    @cache
    def get_iso(self, n, m):
        src = self.get_Y(n, m)
        tgt = self.get_X(n, m)

        send = src.get_iso(tgt)
        return send

    @cache
    def get_hom(self, n):
        # construct Group hom, G0 --> G1
        G0 = self.get_G(n)
        G1 = self.get_G(n+1)
        items = G1.items

        # add new item on the left
        lookup = {src : items[items.index(src)+1] for src in G0.items}
        hom = {}
        for g in G0:
            perm = {lookup[src] : lookup[g.perm[src]] for src in G0.items}
            #print(g.perm, perm)
            for item in items:
                if item not in perm:
                    perm[item] = item
            #print(g.perm, perm)
            h = Perm(perm, items)
            assert h in G1
            hom[g] = h
        check_group_hom(G0, G1, hom)
        #print("hom:", len(hom))
        return hom

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
            g = src.repr[x] # a coset representative for x
            assert src(g)(src.basepoint) == x
            h = self.get_hom(n)[g]
            y = Y(h)(Y.basepoint)
            y = self.get_iso(n+1,m)[y] # send Y --> tgt
            func[x] = y
        check_injection(src.items, tgt.items, func)
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
            g = src.repr[x] # a coset representative for x
            assert src(g)(src.basepoint) == x
            h = self.get_hom(n)[g]
            y = tgt(h)(tgt.basepoint)
            func[x] = y
        check_injection(src.items, tgt.items, func)
        return func

    
class _PascalB(Pascal):

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

    @cache
    def get_Z(self, n, m):
        G = self.get_G(n)
        fix = {G.items[j] for j in range(1, m)}
        fix.add(G.items[n])
        Z = self.get_fix(G, fix)
        return Z

    @cache
    def get_rfunc(self, n, m, idx):
        assert 0<=n
        assert 0<=m<=n
        assert 0<=idx<=1

        src = self.get_X(n, m)
        tgt = self.get_X(n+1, m+1)

        if idx==0:
            # construct a right map
            #print(len(src.items), ">-->", len(tgt.items))
            func = {} # map src --> tgt
            for x in src.items:
                g = src.repr[x] # a coset representative for x
                assert src[g](src.basepoint) == x
                h = self.get_hom(n)[g]
                y = tgt[h](tgt.basepoint)
                func[x] = y
            check_injection(src.items, tgt.items, func)
            return func

        # construct a another right map
        Z = self.get_Z(n+1, m+1)
        iso = Z.get_iso(tgt)
        func = {} # map src --> tgt
        for x in src.items:
            g = src.repr[x] # a coset representative for x
            assert src[g](src.basepoint) == x
            h = self.get_hom(n)[g]
            y = Z[h](Z.basepoint)
            func[x] = iso[y]
        check_injection(src.items, tgt.items, func)
        return func


class PascalB(Pascal):

    @cache
    def get_G(self, n):
        letters = ascii_lowercase[:n]
        items = ["+"+c for c in letters] + list(reversed(["-"+c for c in letters]))
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

    @cache
    def get_hom(self, n):
        # construct Group hom, G0 --> G1
        G0 = self.get_G(n)
        G1 = self.get_G(n+1)
        items = G1.items

        # add new item on the left
        lookup = {}
        for src in G0.items:
            if src[0]=="+":
                tgt = items[items.index(src)+1]
            else:
                tgt = items[items.index(src)-1]
            lookup[src] = tgt
        hom = {}
        for g in G0:
            perm = {lookup[src] : lookup[g.perm[src]] for src in G0.items}
            #print(g.perm, perm)
            for item in items:
                if item not in perm:
                    perm[item] = item
            #print(g.perm, perm)
            h = Perm(perm, items)
            assert h in G1
            hom[g] = h
        check_group_hom(G0, G1, hom)
        #print("hom:", len(hom))
        return hom

    @cache
    def get_Z(self, n, m):
        G = self.get_G(n)
        fix = {G.items[j] for j in range(1, m)}
        fix.add(G.items[-1])
        Z = self.get_fix(G, fix)
        return Z

    @cache
    def get_rfunc(self, n, m, idx):
        assert 0<=n
        assert 0<=m<=n

        src = self.get_X(n, m)
        tgt = self.get_X(n+1, m+1)

        if idx==0:
            # construct a right map
            #print(len(src.items), ">-->", len(tgt.items))
            func = {} # map src --> tgt
            for x in src.items:
                g = src.repr[x] # a coset representative for x
                assert src(g)(src.basepoint) == x
                h = self.get_hom(n)[g]
                y = tgt(h)(tgt.basepoint)
                func[x] = y
            check_injection(src.items, tgt.items, func)
            return func

        # construct a another right map
        Z = self.get_Z(n+1, m+1)
        if idx==1:
            iso = Z.get_iso(tgt)
            func = {} # map src --> tgt
            for x in src.items:
                g = src.repr[x] # a coset representative for x
                assert src(g)(src.basepoint) == x
                h = self.get_hom(n)[g]
                y = Z(h)(Z.basepoint)
                func[x] = iso[y]
            check_injection(src.items, tgt.items, func)
            return func

        print("get_rfunc")
        isos = list(Z.isomorphisms(tgt))
        print("isos:", len(isos))
        #iso = random.choice(isos) # there's two of these, but only one works 

        funcs = []
        for iso in isos:
            func = {} # map src --> tgt
            for x in src.items:
                g = src.repr[x] # a coset representative for x
                assert src[g](src.basepoint) == x
                h = self.get_hom(n)[g]
                y = Z[h](Z.basepoint)
                func[x] = iso[y]
            check_injection(src.items, tgt.items, func)
            funcs.append(func)
        return funcs


def test_weyl_A():

    print("test_weyl_A:")

    triangle = PascalA()
    assert triangle.get_G(3) is triangle.get_G(3)
    assert triangle.get_G(3) is not PascalA().get_G(3)

    N = 6
    print(triangle.str(N))
    print()

    for n in range(N):
        G = triangle.get_G(n)
        items = G.items
        for m in range(1, n):
            src = Space(G, items[:m])
            tgt = Space(G, items[1:m+1])
            iso = src.get_iso(tgt)

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
    print()



def test_weyl_B():

    print("test_weyl_B:")

    triangle = PascalB()

    N = 5
    print(triangle.str(N))
    print()

    for n in range(N):
        G = triangle.get_G(n)
        items = G.items
        for m in range(1, n):
            src = Space(G, items[:m])
            tgt = Space(G, items[1:m+1])
            iso = src.get_iso(tgt)

    for n in range(N):
        hom = triangle.get_hom(n)
        for m in range(n):
            Y = triangle.get_Y(n, m)
            iso = triangle.get_iso(n, m)

    for n in range(2, N):
      for m in range(1, n+1):
        tgt = triangle.get_X(n, m)
        found = set(tgt.items)
        print("(%d %d)=%d"%(n,m,len(found)), end=" ", flush=True)

        rfunc = triangle.get_rfunc(n-1, m-1, 0)
        for item in rfunc.values():
            assert item in found
            found.remove(item)

        rfunc = triangle.get_rfunc(n-1, m-1, 1)
        for item in rfunc.values():
            assert item in found
            found.remove(item)
    
        if m<n:
            lfunc = triangle.get_lfunc(n-1, m)
            for item in lfunc.values():
                assert item in found
                found.remove(item)

        assert not found

      print()
    print()
    print()


def test_32B():
    # show what these functions look like on the Matrix repr of the Weyl's

    get_perm = lambda g : Matrix.perm(g.get_idxs())

    triangle = PascalB()
    B2 = triangle.get_G(2)
    B3 = triangle.get_G(3)
    print("|B2| =", len(B2))
    print(B2.items)
    #B3 = triangle.get_G(3)

    #for g in B3.gens:
    #    M = get_perm(g)
    #    print(M.shortstr())
    #    print()

    def show_func(src, tgt, func):
        for x in src.items:
            smap = SMap()
            col = 0
            g = src.get_repr(x)
            H = src.basepoint
            for h in H:
                M = get_perm(g*h)
                smap[0, col] = M.transpose().shortstr()
                col += len(M)+2
            smap[3, col] = "=>"
            col += 4
            y = func[x]
            
            g = tgt.get_repr(y)
            H = tgt.basepoint
            for h in H:
                M = get_perm(g*h)
                smap[0, col] = M.transpose().shortstr()
                col += len(M)+2
    
            print(smap)
            print()

    X21 = triangle.get_X(2, 1)
    H21 = X21.basepoint
    X32 = triangle.get_X(3, 2)
    H32 = X32.basepoint

    func = triangle.get_rfunc(2, 1, 0) # X21 --> X32
    show_func(X21, X32, func)

    for func in triangle.get_rfunc(2, 1, None): # X21 --> X32
        print()
        show_func(X21, X32, func)


def test_bruhat_order():

    # See:
    # Combinatorics of Coxeter Groups
    # Anders Bjorner
    # Francesco Brenti
    # Chapter 2.

    n = argv.get("n", 3)
    nn = 2*n

    if argv.GL:
        G = Algebraic.GL(n)
    else:
        G = Algebraic.Sp(nn)
    W = G.get_weyl()
    length = {}
    for w in W:
        l = w.name.count('w')
        length[w] = l
    assert length[G.I] == 0

    items = list(W)
    items.sort( key = lambda M : length[M] )
    for M in items:
        print(length[M], end=" ")
    print()
    #for M in items:
    #    print("name:", M.name)
    #    print(M.shortstr())

    print(len(W))
    refls = set()
    for g in W.gen:
        for w in W:
            h = w*g*w.transpose()
            refls.add(h)
    print(len(refls))

    #for i in range(len(refls)+1):
    #    print(list(length.values()).count(i))

    print("gen:", len(W.gen))
    if n==3:
        a, b, c = W.gen
    else:
        a, b = W.gen
    gen = [a]
    H = Algebraic(gen)
    print("|H| =", len(H))

    # the Poset is graded by word lengths
    pairs = set()
    for u in H:
     for v in H:
        if u.transpose()*v in refls and length[u]<length[v]:
            pairs.add((u.shortstr(), v.shortstr())) # u-->v
    hom = Poset.generate(pairs, check=True)
    poset = hom.tgt
    poset.show()


def test_thick():

    n = argv.get("n", 3)

    triangle = PascalB()
    G = triangle.get_G(n)
    items = G.items

    print(items)

    print("|G| =", len(G))

    #for H in G.subgroups():
    #    print("\t", len(H))

    H = []
    a, b, c, aa, bb, cc = items
    #fix = [bb, aa]
    fix = {a, b}
    #print("fix:", fix)
    for g in G:
        send = {g[x] for x in fix}
        if fix == send:
            H.append(g)
            #print("send:", send)
    H = Group(H, items)
    print("|H| =", len(H))
    X = G.left_cosets(H)
    print("|X| =", len(X))

    nn = 2*n
    W0 = Algebraic.Sp(nn-2).get_weyl()
    W = Algebraic.Sp(nn).get_weyl()
    hom = mulclose_hom(G.gens, W.gen)
    check_group_hom(G, W, hom)

    gom = mulclose_hom(W0.gen, W.gen[1:])
    check_group_hom(W0, W, gom)

    smap = SMap()
    for row, coset in enumerate(X):
        #print([(g[a],g[b]) for g in coset])
        for i, g in enumerate(coset):
            g = hom[g]
            g = g.transpose()
            smap[(i+2)*(nn+2), row*(nn+3)] = g.shortstr()

    cols = [0]*len(X)
    for g0 in W0:
        h = gom[g0]
        for row, coset in enumerate(X):
            for i, g in enumerate(coset):
                if h == hom[g]:
                    smap[cols[row]+1, row*(nn+3)+1] = g0.transpose().shortstr()
                    smap[(i+2)*(nn+2), (row+1)*(nn+3)-3] = "*"
                    cols[row] += nn+2

    print(smap)
    print(len(W0))
    print()



def test():
    test_weyl_A()
    test_weyl_B()



if __name__ == "__main__":
    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))






