#!/usr/bin/env python3

"""
Category of finite sets and relations
"""

import sys, os

from bruhat.element import Type, GenericElement



class TpSet(Type):

    def eq(self, a, b):
        return a.set_items == b.set_items

    def ne(self, a, b):
        return a.set_items != b.set_items


class TpRel(Type):

    def eq(self, a, b):
        return a.set_items == b.set_items

    def ne(self, a, b):
        return a.set_items != b.set_items

    def mul(self, a, b):
        "tensor"

    def matmul(self, a, b):
        "dot"

    def promote(self, a):
        if isinstance(a, Rel):
            return a
        items = [(x, x) for x in a]
        return Rel(items)


class Set(object):
    def __init__(self, items):
        for x in items:
            assert type(x) is tuple
        items = list(items)
        items.sort() # eeeeck: watch this ! 
        self.items = tuple(items)
        self.set_items = set(items)
        pairs = [(x, x) for x in items]
        self.ident = Rel(pairs, self, self)

    def __str__(self):
        return "Set(%s)"%(str(list(self.items)))
    __repr__ = __str__

    def __contains__(self, item):
        return item in self.set_items

    def __eq__(a, b):
        return a.items == b.items

    def __ne__(a, b):
        return a.items != b.items

    def __hash__(self):
        return hash(self.items)

    def __mul__(a, b):
        "tensor"
        items = [x+y for x in a.items for y in b.items]
        return Set(items)
        
    def __matmul__(a, b):
        if a==b:
            return a
        else:
            assert 0

    @property
    def left_unitor(self):
        src = Set.one * self
        tgt = self
        items = [(Set.star+x, x) for x in self.items]
        return Rel(items, src, tgt)

    @property
    def right_unitor(self):
        src = self * Set.one
        tgt = self
        items = [(x+Set.star, x) for x in self.items]
        return Rel(items, src, tgt)

    @property
    def cap(self):
        "the unit, eta: I --> X*X"
        src = Set.one
        star = Set.star
        tgt = self * self
        items = [(star, (x + x)) for x in self.items]
        return Rel(items, src, tgt)

    @property
    def cup(self):
        "the co-unit, epsilon: X*X --> I"
        tgt = Set.one
        star = Set.star
        src = self * self
        items = [((x + x), star) for x in self.items]
        return Rel(items, src, tgt)



class Rel(object):

    def __init__(self, items, src, tgt):
        assert isinstance(src, Set)
        assert isinstance(tgt, Set)

        items = list(items)
        items.sort() # eeeeck
        self.items = tuple(items)
        self.set_items = set(items)
        for i, j in items:
            assert i in src, "%s not found in %s" % (i, src) # col
            assert j in tgt, "%s not found in %s" % (j, tgt) # row
        self.src = src
        self.tgt = tgt

    def __str__(self):
        return "Rel(%s)"%(str(list(self.items)))
    __repr__ = __str__

    def __eq__(a, b):
        return a.items == b.items

    def __ne__(a, b):
        return a.items != b.items

    def __hash__(self):
        return hash(self.items)

    def __contains__(self, item):
        return item in self.set_items

    def __mul__(a, b):
        "tensor"
        src = a.src * b.src
        tgt = a.tgt * b.tgt
        items = [(x+y, u+v) for (x, u) in a.items for (y, v) in b.items]
        return Rel(items, src, tgt)

    def __matmul__(a, b):
        assert b.tgt == a.src, "%s != %s" % (b.tgt, a.src)

        items = set()
        for j, i in a.items:
          for k, j1 in b.items:
            if j != j1:
                continue
            items.add((k, i))
        return Rel(items, b.src, a.tgt)

    def dual(self):
        items = [(j, i) for (i, j) in self.items]
        return Rel(items, self.tgt, self.src)


Set.star = ("*",)
Set.one = Set([Set.star])

def dot(*rels):
    idx = 0
    A = rels[idx]
    while idx+1 < len(rels):
        B = rels[idx+1]
        A = A@B
        idx += 1 
    A = A%2
    return A


def compose(*rels):
    rels = list(reversed(rels))
    A = dot(*rels)
    return A



def test():

    A = Set([(x,) for x in "abcd"])
    B = Set([(x,) for x in "uv"])
    C = Set([(x,) for x in "1234"])
    D = Set([(x,) for x in "678"])

    assert (A*B)*C == A*(B*C)

    f = Rel([(('a',), ('u',)), (('a',), ('v',))], A, B) # A--f-->B

    assert f@A.ident == f
    assert B.ident@f == f
    
    g = Rel([(('2',), ('6',)), (('3',), ('8',))], C, D) # C--g-->D

    h = f*g*g

    #print(h)

    # right zig-zag equation
    lhs = A.left_unitor @ ( A.cup * A.ident ) @ ( A.ident * A.cap )
    rhs = A.right_unitor
    assert lhs == rhs

    # left zig-zag equation
    lhs = A.right_unitor @ ( A.ident * A.cup ) @ ( A.cap * A.ident )
    rhs = A.left_unitor
    assert lhs == rhs

    # check dual is transpose
    lhs = A.left_unitor @ ( B.cup * A.ident ) @ (B.ident * f * A.ident ) @ ( B.ident * A.cap )
    rhs = f.dual() @ B.right_unitor
    assert lhs == rhs

    #print(lhs)
    #print(rhs)



if __name__ == "__main__":

    test()



