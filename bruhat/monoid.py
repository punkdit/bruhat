#!/usr/bin/env python

import sys, os

from bruhat.util import choose, cross
from bruhat.argv import argv


def mulclose(els, verbose=False, maxsize=None):
    els = set(els)
    changed = True
    while changed:
        if verbose:
            print("mulclose:", len(els))
        changed = False
        _els = list(els)
        for A in _els:
            for B in _els:
                for C in [A*B, B*A]:
                  if C not in els:
                    els.add(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
                    changed = True
    return els


class Func(object):
    def __init__(self, tgt, src, send, name=''):
        self.send = dict(send)
        self.tgt = list(tgt)
        self.src = list(src)
        self.tgt.sort()
        self.src.sort()
        self.name = name

    def __mul__(self, other):
        "multiply right to left"
        assert self.src == other.tgt
        send = {}
        for k, v in other.send.items():
            v = self.send[v]
            send[k] = v
        return Func(self.tgt, other.src, send, self.name+other.name)

    def __repr__(self):
        send = self.send
        keys = list(send.keys())
        keys.sort()
        return "Func(%s, %s, {%s})"%(
            self.tgt, self.src,
            ', '.join("%s:%s"%(k, send[k]) for k in keys))

    def __str__(self):
        return self.name

    def __eq__(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send == other.send

    def __ne__(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send != other.send

    def __hash__(self):
        return hash(self.__repr__())

    def __lt__(self, other):
        keys = self.src
        lhs = [self.send[k] for k in keys]
        rhs = [other.send[k] for k in keys]
        return lhs < rhs

    @classmethod
    def all_funcs(cls, tgt, src=None):
        if src is None:
            src = tgt
        n = len(src)
        #tgt = tuple(tgt)
        for items in cross([tgt]*n):
            send = dict(zip(src, items))
            yield cls(tgt, src, send)

    @classmethod
    def identity(cls, items):
        return cls(items, items, dict(zip(items, items)))



class Monoid(object):
    "_concrete monoid"
    def __init__(self, funcs):
        funcs = list(funcs)
        funcs.sort()
        self.funcs = tuple(funcs)

    def __str__(self):
        return "Monoid(%s)"%(self.funcs,)

    def __len__(self):
        return len(self.funcs)

    def __getitem__(self, i):
        return self.funcs[i]

    def __eq__(self, other):
        return self.funcs == other.funcs

    def __hash__(self):
        return hash(self.funcs)

    def check(self):
        funcs = set(self.funcs)
        assert funcs
        for f in funcs:
          for g in funcs:
            assert f*g in funcs
        items = f.src
        I = Func.identity(items) 
        assert I in funcs

    @classmethod
    def generate(cls, funcs, *args, **kw):
        funcs = list(mulclose(funcs, *args, **kw))
        return cls(funcs)

    def find(self, f):
        for x in self.funcs:
            if x==f:
                return x


def test_cayley():

    n = 3

    items = list(range(1, n+1))
    tgt = src = items

    I = Func(tgt, src, dict((i, i) for i in items), "I")
    a = Func(tgt, src, dict((i, min(i+1, n)) for i in items), "a")
    b = Func(tgt, src, dict((i, max(i-1, 1)) for i in items), "b")

    M = Monoid.generate([I, a, b])
    print(n, len(M))

    # left cayley
    f = open("left_monoid.dot", "w")
    print("digraph {", file=f)
    for x in M:
        ax = M.find(a*x)
        bx = M.find(b*x)
        print("%s -> %s [label=a];" % (x, ax), file=f)
        print("%s -> %s [label=b];" % (x, bx), file=f)
    print("}", file=f)
    f.close()

    # right cayley
    f = open("right_monoid.dot", "w")
    print("digraph {", file=f)
    for x in M:
        xa = M.find(x*a)
        xb = M.find(x*b)
        print("%s -> %s [label=a];" % (x, xa), file=f)
        print("%s -> %s [label=b];" % (x, xb), file=f)
    print("}", file=f)
    f.close()





def test_sequence():

    # this is OEIS A081489
    for n in range(1, 7):
        items = list(range(1, n+1))
        tgt = src = items
    
        I = Func(tgt, src, dict((i, i) for i in items), "")
        up = Func(tgt, src, dict((i, min(i+1, n)) for i in items), "u")
        dn = Func(tgt, src, dict((i, max(i-1, 1)) for i in items), "d")
    
        M = Monoid.generate([I, up, dn], maxsize=1200)
        print(n, len(M))



def test():

    n = argv.get("n", 3)
    items = list(range(1, n+1))
    tgt = src = items

    I = Func(tgt, src, dict((i, i) for i in items), "")
    up = Func(tgt, src, dict((i, min(i+1, n)) for i in items), "u")
    dn = Func(tgt, src, dict((i, max(i-1, 1)) for i in items), "d")

    M = Monoid.generate([I, up, dn], maxsize=12)
    assert len(M) == 8

    # ---------------------------------------------------

    X = [0, 1, 2]

    M = [
        Func(X, X, {0:0, 1:1, 2:2}, 'I'),
        Func(X, X, {0:0, 1:0, 2:1}, 'a'),
        Func(X, X, {0:0, 1:0, 2:2}, 'b'),
        Func(X, X, {0:0, 1:1, 2:1}, 'c'),
        Func(X, X, {0:1, 1:1, 2:2}, 'd'),
        Func(X, X, {0:0, 1:2, 2:2}, 'e'),
        Func(X, X, {0:1, 1:2, 2:2}, 'f'),
        Func(X, X, {0:0, 1:0, 2:0}, 'g'),
        Func(X, X, {0:1, 1:1, 2:1}, 'h'),
        Func(X, X, {0:2, 1:2, 2:2}, 'i'),
    ]

    N = Monoid.generate(M, maxsize=100)
    assert len(N) == len(M)

    I, a, b, c, d, e, f, g, h, i = M
    assert repr(b*a) == "Func([0, 1, 2], [0, 1, 2], {0:0, 1:0, 2:0})"
    assert b*a == g
    assert a*b == a

    for gen in choose(M[1:], 3):
        N = Monoid.generate(gen)
        if len(N) == len(M)-1:
            print([str(g) for g in gen], end=" ")
            print(len(N))
    
    gens = [a, b, f]
    show_table(M, M)


def show_table(M, gens):
    lookup = dict((f, f) for f in M)
    print(" |", end="")
    for b in M:
        print(b, end=" ")
    print()
    print("-+" + "--"*len(M))
    for b in M:
      print(b, end="|")
      for a in M:
        c = lookup[a*b]
        if a in gens and b in gens:
            print(c, end=" ")
        else:
            print(" ", end=" ")
      print()



def all_submonoids(M):

    funcs = list(M)
    f = funcs[0]
    items = f.src

    I = Func.identity(items)
    funcs.remove(I)

    #print(len(funcs))

    M = Monoid([I])
    yield M
    bdy = [M]
    found = set(bdy)
    while bdy:
        _bdy = []
        for M0 in bdy:
            avoid = set(M0)
            gen = list(M0)
            for f in funcs:
                if f in avoid:
                    continue
                M = Monoid.generate(gen+[f])
                if M in found:
                    continue
                yield M
                found.add(M)
                _bdy.append(M)
                #print(len(M), end=' ', flush=True)
        bdy = _bdy
        #print("%d:%d"%(len(found), len(bdy)))
    #print("found:", len(found))




def test_submonoids():

    # https://oeis.org/A343140
    # Number of submonoids of the monoid of maps from an n-element set to itself.

    for n in [2,3]:
        items = list(range(n))
    
        funcs = []
        for f in Func.all_funcs(items, items):
            funcs.append(f)

        M = Monoid(funcs)

        count = 0
        for M1 in all_submonoids(M):
            count += 1
        assert count == [1,6,699][n-1]


def test_endo():

    # The image of M in End(M)
    # is characterised as the functions M->M
    # that commute with all right translations.

    n = 3
    items = list(range(n))

    funcs = []
    for f in Func.all_funcs(items, items):
        funcs.append(f)

    M1 = Monoid(funcs)

    for M in all_submonoids(M1):
        M.check()
        if len(M) == 6:
            print(M)
            test_regular(M)


def test_regular(M):
    X = M.funcs

    left = []
    right = []
    for m in M:
        send = dict((i, m*i) for i in X)
        f = Func(X, X, send)
        left.append(f)
        send = dict((i, i*m) for i in X)
        f = Func(X, X, send)
        right.append(f)
    L = Monoid(left)
    L.check()
    assert len(L) == len(M)
    R = Monoid(right)
    R.check()
    assert len(R) == len(M)

    #print(L == R)

    found = []
    for f in Func.all_funcs(X):
        for g in right:
            if f*g != g*f:
                break
        else:
            found.append(f)
    #print(len(found))

    M = Monoid(found)
    M.check()
    assert M == L
    
        



if __name__ == "__main__":
    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "test"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()






