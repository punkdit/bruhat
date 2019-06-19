#!/usr/bin/python3

"""

some code taken from huffman.py

"""

import sys
from functools import reduce
from operator import add
from math import log, log2
from random import shuffle, choice, randint, seed

from bruhat.argv import argv


EPSILON = 1e-8

def is_close(a, b):
    return abs(a-b) < EPSILON



class Multiset(object):
    "un-normalized probability distribution"

    def __init__(self, cs={}, tp=None):
        for k in cs.keys():
            assert tp is None or tp is type(k)
            tp = type(k)
        self.tp = tp
        items = [(k, v) for (k, v) in cs.items() if v>0]
        cs = dict(items)
        self.cs = cs # map item -> count
        self._len = sum(self.cs.values(), 0)
        keys = list(cs.keys())
        keys.sort() # canonicalize
        self._keys = keys 
        self.val = tuple((k, cs[k]) for k in keys) # for hash, __lt__, ...

    def __str__(self):
        cs = self.cs
        keys = self._keys
        items = reduce(add, [(str(key),)*cs[key] for key in keys], ())
        items = '+'.join(items)
        return '(%s)'%items
    __repr__ = __str__

    def get_str(self):
        cs = self.cs
        keys = self._keys
        items = [(key if cs[key]==1 else "%d%s"%(cs[key], key)) for key in keys]
        items = '+'.join(items)
        return items

    def check_tp(self, other):
        assert self.tp is None or other.tp is None or self.tp is other.tp, "%s =!= %s"%(self, other)

    def __eq__(self, other):
        self.check_tp(other)
        return self.cs == other.cs

    def __ne__(self, other):
        self.check_tp(other)
        return self.cs != other.cs

    def __lt__(self, other):
        self.check_tp(other)
        return self.val < other.val

    def __hash__(self):
        return hash(self.val)

    def __mul__(X, Y): 
        "cartesian product of multisets"
        if not isinstance(Y, Multiset):
            return NotImplemented

        if not X.cs or not Y.cs:
            return Multiset() # zero

        if X.tp is str and Y.tp is str:
            xcs, ycs = X.cs, Y.cs
            cs = dict(((x+y), xcs[x]*ycs[y]) for x in xcs for y in ycs)
            return Multiset(cs, str)

        if X.tp is Box and Y.tp is str:
            item = Multiset({}, Box)
            for box in X.terms():
                item = item + box*Y
            return item

        if X.tp is str and Y.tp is Box:
            item = Multiset({}, Box)
            for box in Y.terms():
                item = item + Multiset({X*box : 1}, Box)
            return item

        if X.tp is Box and Y.tp is Box:
            item = Multiset({}, Box)
            for x_box in X.terms():
              for y_box in Y.terms():
                item = item + x_box * y_box
            return item

        assert 0, "%s =X= %s" % (X, Y)


    def __rmul__(self, r):
        "left multiplication by a number"
        assert int(r) == r
        assert r >= 0
        cs = dict((k, r*v) for (k, v) in self.cs.items())
        return Multiset(cs, self.tp)

    def __add__(X, Y):
        # WARNING: not disjoint union (coproduct)
        Y = Multiset.promote(Y)
        X.check_tp(Y)
        xcs, ycs = X.cs, Y.cs
        cs = dict(xcs)
        for k, v in ycs.items():
            cs[k] = cs.get(k, 0) + v
        return Multiset(cs, X.tp)

    def __getitem__(self, k):
        return self.cs[k]

    def keys(self):
        return self._keys

    def items(self):
        return self.cs.items()

    def terms(self):
        cs = self.cs
        for k in self._keys:
            for i in range(cs[k]):
                yield k

    def disjoint(X, Y):
        # We only keep non-zero keys, so this works
        lhs = set(X.cs.keys())
        rhs = set(Y.cs.keys())
        return not bool(lhs.intersection(rhs))

    def contains(self, other):
        "self contains other"
        cs = self.cs
        for k,v in other.cs.items():
            if v > cs.get(k, 0):
                return False
        return True

    def __len__(self):
        return self._len

    def isomorphic(self, other):
        if self._len != other._len:
            return False
        lhs = set(self.cs.values())
        rhs = set(other.cs.values())
        return lhs == rhs

    @classmethod
    def promote(cls, item):
        if isinstance(item, Multiset):
            pass
        else:
            item = Multiset({item:1})
        return item


class Box(object):
    def __init__(self, X):
        assert isinstance(X, Multiset)
        assert X.tp is None or X.tp is str
        self.X = X

    def __str__(self):
        return "Box%s"%(self.X,)
    __repr__ = __str__

    def __eq__(self, other):
        return self.X == other.X

    def __ne__(self, other):
        return self.X != other.X

    def __lt__(self, other):
        return self.X < other.X

    def __hash__(self):
        return hash(self.X)

    def right_multiply(self, Y):
        "right multiply by a multiset"
        assert isinstance(Y, Multiset)
        X = self.X
        cs = {}
        for k, v in Y.items():
            k = X*Multiset({k:1})
            box = Box(k)
            cs[box] = v
        return Multiset(cs)

    def __add__(self, other):
        self = Multiset.promote(self)
        other = Multiset.promote(other)
        return self + other

    def __mul__(self, other):
        if isinstance(other, Multiset):
            item = self.right_multiply(other)
        elif isinstance(other, Box):
            item = self.X * other + self * other.X
        elif isinstance(other, int):
            item = Multiset({self : other})
        else:
            assert 0
        return item

    def __rmul__(self, other):
        if isinstance(other, Multiset):
            item = Box(other*self.X)
        elif isinstance(other, int):
            item = Box(other*self.X)
        else:
            assert 0
        return item



def main():

    X = Multiset({"a":3, "b":1})

    assert (X+X) == Multiset({"a":6, "b":2})
    assert (X+X) == 2*X
    #print(X, X.entropy())

    XX = X*X

    Y = Multiset({"a":2, "b":2})
    #print(Y, Y.entropy())
    assert str(Y) == "(a+a+b+b)"

    a = Multiset({"a" : 1})
    b = Multiset({"b" : 1})
    c = Multiset({"c" : 1})
    d = Multiset({"d" : 1})
    e = Multiset({"e" : 1})
    f = Multiset({"f" : 1})
    g = Multiset({"g" : 1})

    assert a.disjoint(b)
    assert not (a+b).disjoint(b)

    assert list((a+2*b).terms()) == ["a", "b", "b"]

    assert not a.contains(b)
    assert (a+b).contains(b)
    assert not (a+b).contains(2*b)

    def mkrand():
        X = randint(0,2)*a + randint(0,2)*b + randint(0,2)*c + randint(0,2)*d
        return X

    zero = Multiset()

    X = a+2*b
    BX = Box(X)
    Y = c
    #print("%s*%s = %s"%(Y, BX, Y*BX))
    #print("%s*%s = %s"%(BX, Y, BX*Y))
    #print("%s*%s = %s"%(BX, (c+c), BX*(c+c)))

    # test distributivity
    for trial in range(100):
        X = mkrand()
        Y = mkrand()
        Z = mkrand()
        assert (X*Y) * Box(Z) == X*(Y*Box(Z))

        lhs = Box(X) * (Y*Z)
        rhs = (Box(X)*Y)*Z
        assert lhs == rhs

        if X!=zero:
            lhs = (X*Box(Y))*Z
            rhs = X*(Box(Y)*Z)
            assert lhs == rhs, "%s != %s"%(lhs, rhs)

    def mkbox(count=3):
        item = Multiset({}, Box)
        for i in range(randint(0, count)):
            b = mkrand()
            if b == zero:
                continue
            item += Box(b)
        return item

    for trial in range(100):
        A = mkbox()
        B = mkbox()
        C = mkbox()

        lhs = (A+B) * C
        rhs = A*C + B*C
        assert lhs == rhs

        lhs = A*(B+C)
        rhs = A*B + A*C
        assert lhs == rhs

    def strip(s):
        s = str(s)
        s = s.replace("Box", "")
        s = s.replace("(", "")
        s = s.replace(")", "")
        s = s.split("+")
        s.sort()
        return '+'.join(s)

    #zero = Multiset({}, Box)
    n = 10000
    #seed(0)
    for trial in range(1000):
        A = mkbox(1)
        B = mkbox(1)
        C = mkbox(1)
        #if A==zero or B==zero or C==zero:
            #continue
        lhs = (A*B)*C
        rhs = A*(B*C)
        #print(lhs)
        #print(rhs)
        if lhs == rhs:
            continue
        assert strip(lhs) == strip(rhs)
#        s = str(lhs)
#        if len(s) < n:
#            n = len(s)
#            #print(A, B, C)
#            print(strip(lhs))
#            print(strip(rhs))
#            print()

    a = Multiset({"a" : 1})
    b = Multiset({"b" : 1})
    c = Multiset({"c" : 1})
    d = Multiset({"d" : 1})
    A = Multiset({Box(a):1})
    B = Multiset({Box(b):1})
    C = Multiset({Box(c+d):1})
    print((A*B)*C)
    print(A*(B*C))


if __name__ == "__main__":

    main()


