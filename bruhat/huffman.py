
"""
https://golem.ph.utexas.edu/category/2019/03/how_much_work_can_it_be_to_add.html#c055688

see also entropy.py

"""

import sys
from functools import reduce
from operator import add
from math import log
from random import shuffle, choice

#import numpy
#from matplotlib import pyplot

#from bruhat.gelim import row_reduce, shortstr, kernel
#from qupy.dev import linalg

from bruhat.argv import argv

EPSILON = 1e-8

def is_close(a, b):
    return abs(a-b) < EPSILON

def W(items):

    sitems = sum(items)
    r = 0.
    for n in items:
        r += n * log(n)
    return -1*(r - sitems*log(sitems))


class Multiset(object):
    "un-normalized probability distribution"

    def __init__(self, cs={}):
        self.cs = dict(cs) # map item -> count

    def __str__(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()
        items = reduce(add, [ (str(key),)*cs[key] for key in keys ])
        items = ','.join(items)
        return '{%s}'%items

    def __eq__(self, other):
        return self.cs == other.cs

    def __mul__(X, Y): # cartesian product... yes?
        xcs, ycs = X.cs, Y.cs
        cs = dict((x+y, xcs[x]*ycs[y]) for x in xcs for y in ycs)
        return Multiset(cs)

    def __rmul__(self, r):
        assert int(r) == r
        cs = dict((k, r*v) for (k, v) in self.cs.items())
        return Multiset(cs)

    def __add__(X, Y): # WARNING: not disjoint union
        xcs, ycs = X.cs, Y.cs
        cs = dict(xcs)
        for k, v in ycs.items():
            cs[k] = cs.get(k, 0) + v
        return Multiset(cs)

    def __len__(self):
        return sum(self.cs.values(), 0)

    def entropy(self):
        "un-normalized entropy"
        cs = self.cs
        items = [n for n in cs.values() if n>0]
        return W(items)

    def huffman(self):
        cs = self.cs
        keys = list(cs.keys())

        # build a tree, start with the leaves:
        nodes = [Node(Multiset({key : cs[key]})) for key in keys]

        while len(nodes) > 1:
            n = len(nodes)
            best = (0, 1)
            value = nodes[0].cost() + nodes[1].cost()
            for i in range(n):
              for j in range(i+1, n):
                w = nodes[i].cost() + nodes[j].cost()
                if w < value:
                    best = (i, j)
                    value = w
            i, j = best
            assert i < j, (i, j)
            right = nodes.pop(j)
            left = nodes.pop(i)
            node = Node(left.X + right.X, left, right)
            nodes.append(node)

        return nodes[0]

    def total_length(self):
        n = sum([len(k)*v for (k, v) in self.cs.items()], 0)
        return n
    


class Node(object):
    def __init__(self, X, left=None, right=None):
        self.X = X
        self.left = left
        self.right = right

    def cost(self):
        return len(self.X)

    def encode(self):
        " the (un-normalized) distribution of encoded words "
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            return Multiset({'' : len(X)})
        left = left.encode()
        right = right.encode()
        left = Multiset(dict(('0'+k, v) for (k, v) in left.cs.items()))
        right = Multiset(dict(('1'+k, v) for (k, v) in right.cs.items()))
        return left + right

    def __str__(self):
        X = self.X
        left = self.left
        right = self.right

        if left is None and right is None:
            return str(X)
        assert left and right
        return "(%s : %s)" % (left, right)



def main():

    X = Multiset({"A":3, "B":1})

    assert (X+X) == Multiset({"A":6, "B":2})
    assert (X+X) == 2*X
    #print(X, X.entropy())

    XX = X*X
    #print(XX, XX.entropy())

    Y = Multiset({"A":2, "B":2})
    #print(Y, Y.entropy())

    #print( ((X*Y).entropy(), len(X)*Y.entropy() + len(Y)*X.entropy()))
    assert is_close((X*Y).entropy(), len(X)*Y.entropy() + len(Y)*X.entropy())

    tree = X.huffman()
    assert tree.X == X
    print(tree)
    print(tree.encode())

    tree = (X*X).huffman()
    print(tree.encode())
    assert tree.encode().total_length() == 27




if __name__ == "__main__":

    main()


