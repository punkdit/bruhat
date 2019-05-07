
"""
https://golem.ph.utexas.edu/category/2019/03/how_much_work_can_it_be_to_add.html#c055688

see also entropy.py

"""

import sys
from functools import reduce
from operator import add
from math import log, log2
from random import shuffle, choice, randint, seed

#import numpy
#from matplotlib import pyplot

#from bruhat.gelim import row_reduce, shortstr, kernel
#from qupy.dev import linalg

from bruhat.argv import argv


EPSILON = 1e-8

def is_close(a, b):
    return abs(a-b) < EPSILON


def entropy(items, base=2):
    "un-normalized entropy"
    k = log(base)
    sitems = sum(items)
    r = 0.
    for n in items:
        r += n * log(n) / k
    return -1*(r - sitems*log(sitems) / k)


def entropy(items):
    "un-normalized entropy"
    sitems = sum(items)
    r = 0.
    for n in items:
        r += n * log2(n)
    return -1*(r - sitems*log2(sitems))


class Multiset(object):
    "un-normalized probability distribution"

    def __init__(self, cs={}):
        self.cs = dict(cs) # map item -> count
        self._len = sum(self.cs.values(), 0)

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
        return self._len

    def entropy(self):
        "un-normalized entropy"
        cs = self.cs
        items = [n for n in cs.values() if n>0]
        return entropy(items)

    def huffman(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()

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

    def product_tree(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()
        for k in keys:
            assert len(k) == 2
            # fail..

    def total_length(self):
        n = sum([len(k)*v for (k, v) in self.cs.items()], 0)
        return n

    def W(self): # brain fart the name
        return self.huffman().encode().total_length()
    

def W(X):
    return X.W()


class Node(object):
    def __init__(self, X, left=None, right=None):
        self.X = X
        self._cost = len(X)
        self.left = left
        self.right = right

    def cost(self):
        return self._cost

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

#    def __mul__(S, T):
#        # argh...



def main():

    X = Multiset({"A":3, "B":1})

    assert (X+X) == Multiset({"A":6, "B":2})
    assert (X+X) == 2*X
    #print(X, X.entropy())

    XX = X*X

    Y = Multiset({"A":2, "B":2})
    #print(Y, Y.entropy())

    #print( ((X*Y).entropy(), len(X)*Y.entropy() + len(Y)*X.entropy()))
    assert is_close((X*Y).entropy(), len(X)*Y.entropy() + len(Y)*X.entropy())

    tree = X.huffman()
    assert tree.X == X

    #assert str(tree) == "({B} : {A,A,A})", repr(str(tree))
    #assert str(tree.encode()) == "{0,1,1,1}"

    tree = XX.huffman()
    assert str(tree.encode()) == "{0,0,0,0,0,0,0,0,0,10,10,10,110,110,110,111}"

    assert XX.W() == 27

    A = Multiset({"A" : 1})
    B = Multiset({"B" : 1})
    C = Multiset({"C" : 1})
    D = Multiset({"D" : 1})
    E = Multiset({"E" : 1})
    F = Multiset({"F" : 1})

    def mkrand(a=1, b=3):
        Z = randint(a, b)*A + randint(a, b)*B + randint(a, b)*C + randint(a, b)*D + randint(a, b)*E
        return Z

    for trial in range(100):
        X = mkrand()
        n = randint(2, 5)
        assert n*X.W() == (n*X).W()
    #    print(Z.entropy(), Z.W())

    X = 3*A + B
    #print(X.huffman())
    lhs, rhs = (X*X).W(), len(X)*X.W() + len(X)*X.W()
    #print(lhs, rhs)
    assert lhs < rhs
    assert lhs == 27
    assert rhs == 32

    for trial in range(100):
        X = mkrand(1, 3)
        Y = mkrand(1, 3)
        lhs, rhs = (X*Y).W(), len(X)*Y.W() + len(Y)*X.W()
        assert lhs<=rhs

        lhs, rhs = (X*Y).entropy(), len(X)*Y.entropy() + len(Y)*X.entropy()
        assert is_close(lhs, rhs)

#    Z = 25*A + 25*B + 20*C + 15*D + 15*E

    def mkrand(items, a=1, b=3):
        Z = Multiset()
        for A in items:
            Z = Z + randint(a, b)*A
        return Z

    for trial in range(100):

        X = mkrand([A, B, C])
        Y = mkrand([D, E, F])

        #print(X, Y)
        #print(X+Y)

        lhs = W(X+Y)
        rhs = W(X + len(Y)*D) + W(Y)
        #print(lhs, rhs)
        assert lhs <= rhs

        lhs = (X+Y).entropy()
        rhs = (X + len(Y)*D).entropy() + (Y).entropy()
        assert is_close(lhs, rhs)
        #print(lhs, rhs)
        #print()


        #break
        
    
    #X = 5*A + 1*B + C + 1*D
    X = 1*A + 1*B + 1*C
    X1 = 1
    h = X.entropy()
    print(h)
    n = 1
    while 1:
        X1 = X1*X
        lhs, rhs = (X1.W(), X1.entropy())
        r = n * (len(X)**(n-1))
        assert is_close(rhs/r, h)
        print('\t', lhs / r)
#        if len(X1) > 10000000:
#            break
        n += 1

    print(len(X1))




if __name__ == "__main__":

    main()


