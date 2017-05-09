#!/usr/bin/env python

import sys, os

from bruhat.argv import argv


def mulclose(els, verbose=False, maxsize=None):
    els = set(els)
    changed = True
    while changed:
        if verbose:
            print "mulclose:", len(els)
        changed = False
        _els = list(els)
        for A in _els:
            for B in _els:
                C = A*B
                if C not in els:
                    els.add(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
                    changed = True
    return els


class Func(object):
    def __init__(self, func, word=''):
        self.func = dict(func)
        self.word = word

    def __mul__(self, other):
        func = {}
        keys = self.func.keys()
        for k, v in self.func.items():
            v = other.func[v]
            func[k] = v
        return Func(func, self.word+other.word)

    def __str__(self):
        func = self.func
        keys = func.keys()
        keys.sort()
        return "Func({%s})"%(
            ', '.join("%s:%s"%(k, func[k]) for k in keys))

    def __eq__(self, other):
        return self.func == other.func

    def __ne__(self, other):
        return self.func != other.func

    def __hash__(self):
        func = self.func
        keys = func.keys()
        keys.sort()
        items = tuple((k, func[k]) for k in keys)
        return hash(items)


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
    def generate(cls, funcs):
        funcs = list(mulclose(funcs))
        return cls(funcs)



def test():

    n = argv.get("n", 3)
    items = range(1, n+1)

    I = Func(dict((i, i) for i in items), "")
    up = Func(dict((i, min(i+1, n)) for i in items), "u")
    dn = Func(dict((i, max(i-1, 1)) for i in items), "d")

    M = Monoid.generate([I, up, dn])

    print len(M)

    for x in M:
        print x.word,
    print


if __name__ == "__main__":

    test()




