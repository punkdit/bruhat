#!/usr/bin/env python

import sys, os

from bruhat.util import choose
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


class Monoid(object):
    "concrete monoid"
    def __init__(self, funcs):
        self.funcs = list(funcs)
        for func in funcs:
            assert isinstance(func, Func)

    def __len__(self):
        return len(self.funcs)

    def __getitem__(self, i):
        return self.funcs[i]

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






if __name__ == "__main__":

    #test_sequence()
    #test()
    test_cayley()

    print("OK\n")




