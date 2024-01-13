#!/usr/bin/env python3

"""

"""


from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import add, mul

from bruhat.equ import Equ
from bruhat.util import factorial, cross, is_prime
from bruhat.action import mulclose
from bruhat.gset import Simplicial, Cone
from bruhat.argv import argv


class Set(object):
    def __init__(self, items):
        items = list(items)
        #items.sort()
        self.items = items
        self.set_items = set(items)

    def __str__(self):
        return "Set(%s)"%(self.items,)
    __repr__ = __str__

    def __contains__(self, item):
        return item in self.set_items

    def __iter__(self):
        return iter(self.items)

    # __eq__ is '=='

    @property
    def i(self):
        send_items = {i:i 
        return Function(self, self, send_items)



class Function(object):
    def __init__(self, tgt, src, send_items, check=True): # because __mul__ is right-to-left
        assert isinstance(tgt, Set)
        assert isinstance(src, Set)
        send_items = dict(send_items)
        self.data = src, tgt, send_items
        self.src = src
        self.tgt = tgt
        self.send_items = send_items
        if check:
            self.check()

    def check(self):
        src, tgt, send_items = self.data
        for i in src:
            assert i in send_items
            assert send_items[i] in tgt

    def __mul__(self, other):
        assert self.src == other.tgt
        send_items = {i:self.send_items[other.send_items[i]] for i in other.src}
        return Function(self.tgt, other.src, send_items)



def test():
    X = Set([0,1,2,3])
    Y = Set([0,1])
    f = Function(Y, X, {0:0, 1:0, 2:1, 3:0})
    g = Function(Y, Y, {0:1, 1:0})
    g*f



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






