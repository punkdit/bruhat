#!/usr/bin/env python3

"""
Build a pair of adjoint functors (galois connection) out of a relation R <= AxB.

"""


from functools import reduce, cache
from operator import add, mul
from random import shuffle, choice

import numpy # why not

from bruhat.argv import argv


# modified from bruhat.species
class Set:
    def __init__(self, items=[]):
        items = list(items)
        items.sort(key = str) # <--- use str to canonicalize
        self.items = tuple(items) # canonical order
        self.set_items = set(items)
        assert len(self.set_items) == len(self.items), self.items

    @classmethod
    def promote(cls, items):
        if isinstance(items, Set):
            return items
        items = Set(items)
        return items

    def __str__(self):
        return "{%s}" % (', '.join(str(x) for x in self))
    __repr__ = __str__

    #def __repr__(self):
    #    return "Set(%s)"%(str(list(self.items)))

    @classmethod
    def promote(cls, items):
        if isinstance(items, Set):
            return items
        #if type(items) is GeneratorType:
        return Set(items)

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

    def __contains__(self, item):
        return item in self.set_items

    # use promote below ? XXX
    def __eq__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items == b.items

    def __ne__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items != b.items

    def __lt__(a, b):
        assert isinstance(b, Set), repr(b)
        lhs = a.set_items
        rhs = b.set_items
        return lhs != rhs and lhs.issubset(rhs) 

    def __le__(a, b):
        assert isinstance(b, Set), repr(b)
        lhs = a.set_items
        rhs = b.set_items
        return lhs.issubset(rhs)

    def __hash__(self):
        return hash(self.items)

    def __add__(a, b): # XXX what about universal arrows ?
        # categorical coproduct
        items = [(0, ai) for ai in a] + [(1, bi) for bi in b]
        return Set(items)

    def union(a, b):
        items = a.set_items.union(b.set_items)
        return Set(items)

    def __mul__(a, b): # XXX what about universal arrows ?
        # categorical product
        items = [(ai, bi) for ai in a for bi in b]
        return Set(items)
        
    @cache
    def power(self):
        items = self.items
        n = len(items)
        jtems = []
        for bits in numpy.ndindex((2,)*n):
            subset = [items[i] for i in range(n) if bits[i]]
            subset = Set(subset)
            jtems.append(subset)
        return Set(jtems)

    def ident(self):
        pairs = [(a,a) for a in self]
        return Rel(self, self, pairs)

    @cache
    def all_rels(A, B):
        rels = []
        AB = A*B
        for pairs in AB.power():
            rel = Rel(A, B, pairs)
            rels.append(rel)
        return rels
        


class Rel:
    # rename tgt,src as lhs,rhs ??? 
    def __init__(self, tgt, src, pairs):
        assert isinstance(tgt, Set)
        assert isinstance(src, Set)
        pairs = Set.promote(pairs)
        for (a,b) in pairs:
            assert a in tgt
            assert b in src
        self.pairs = pairs
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src)

    def __len__(self):
        return len(self.pairs)

    def __contains__(self, pair):
        return pair in self.pairs

    def __iter__(self):
        return iter(self.pairs)

    def __repr__(self):
        return "Rel(%s, %s, %s)"%(self.tgt, self.src, self.pairs)

    def __str__(self):
        s = ", ".join("%s~%s"%(a,b) for (a,b) in self.pairs)
        return "%s%s[%s]"%(self.tgt, self.src, s)

    def __mul__(self, other):
        assert self.src == other.tgt
        pairs = set()
        for (i,j) in self.pairs:
          for (k,l) in other.pairs:
            if j==k:
                pairs.add((i,l))
        return Rel(self.tgt, other.src, pairs)

    def __eq__(self, other):
        assert self.hom == other.hom, " update __hash__ if you remove this "
        return self.pairs == other.pairs

    def __hash__(self):
        return hash(self.pairs) # hmm, forget the hom?

    def __le__(self, other):
        assert self.hom == other.hom
        return self.pairs <= other.pairs

    def __lt__(self, other):
        assert self.hom == other.hom
        return self.pairs < other.pairs

    @property
    @cache
    def op(self):
        pairs = [(b,a) for (a,b) in self.pairs]
        return Rel(self.src, self.tgt, pairs)

    def send(rel):
        "also known as 'image'"
        F = {}
        A, B = rel.hom
        for A1 in A.power():
            B1 = []
            for b in B:
                for a in A1:
                    if (a,b) in rel:
                        B1.append(b)
                        break
            F[A1] = Set(B1)
        return F

    def pull(rel):
        "also known as 'inverse image'"
        G = {}
        A, B = rel.hom
        for Y in B.power():
            X = []
            for a in A:
                succ = True
                for b in B:
                    if (a,b) in rel and b not in Y:
                        succ = False
                if succ:
                    X.append(a)
            G[Y] = Set(X)
        return G


        





def test():

    A = Set("a0 a1 a2".split())
    #print(A)
    PA = A.power()
    assert len(PA) == 8

    B = Set("b0 b1 b2".split())
    PB = B.power()
    #print(B)

    i = A.ident()
    #print(i)
    assert i*i == i
    assert i.op == i

    # ------------------------------

    #A = Set("a0 a1".split())
    #PA = A.power()
    #B = Set("b0 b1".split())
    #PB = B.power()

    rels = A.all_rels(B)
    #shuffle(rels)
    #rel = rels[5]

    for rel in rels:

        SEND = rel.send()
        PULL = rel.pull()
    
        for X in PA:
            for Y in PB:
                lhs = int(SEND[X] <= Y)
                rhs = int(X <= PULL[Y])
                if lhs == rhs:
                    continue

                print("R =", rel)
                print("SEND =", SEND)
                print("PULL =", PULL)
    
                print("X =", X)
                print("Y =", Y)
                print(SEND[X], Y, SEND[X]<=Y)
                print(X, PULL[Y], X<=PULL[Y])
    
                assert 0



    

if __name__ == "__main__":
    from time import time
    start_time = time()

    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))


