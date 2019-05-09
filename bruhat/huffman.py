#!/usr/bin/python3

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
        items = [(k, v) for (k, v) in cs.items() if v>0]
        cs = dict(items)
        self.cs = dict(cs) # map item -> count
        self._len = sum(self.cs.values(), 0)
        keys = list(cs.keys())
        keys.sort() # canonicalize
        self.keys = keys 

    def __str__(self):
        cs = self.cs
        keys = self.keys
        items = reduce(add, [(str(key),)*cs[key] for key in keys], ())
        items = ','.join(items)
        return '{%s}'%items

    __repr__ = __str__

    def __eq__(self, other):
        return self.cs == other.cs

    def __ne__(self, other):
        return self.cs != other.cs

    def __mul__(X, Y): 
        "cartesian product of multisets"
        if isinstance(Y, Multiset):
            xcs, ycs = X.cs, Y.cs
            cs = dict((x+y, xcs[x]*ycs[y]) for x in xcs for y in ycs)
            return Multiset(cs)
        return NotImplemented

    def __rmul__(self, r):
        "left multiplication by a number"
        assert int(r) == r
        assert r >= 0
        cs = dict((k, r*v) for (k, v) in self.cs.items())
        return Multiset(cs)

    def __add__(X, Y):
        # WARNING: not disjoint union (coproduct)
        xcs, ycs = X.cs, Y.cs
        cs = dict(xcs)
        for k, v in ycs.items():
            cs[k] = cs.get(k, 0) + v
        return Multiset(cs)

    def terms(self):
        cs = self.cs
        return [Multiset({k:cs[k]}) for k in self.keys]

    def disjoint(X, Y):
        lhs = set(X.cs.keys())
        rhs = set(Y.cs.keys())
        #print("disjoint:", lhs, rhs)
        #print("disjoint:", lhs.intersection(rhs))
        return not bool(lhs.intersection(rhs))

    def __len__(self):
        return self._len

    def isomorphic(self, other):
        if self._len != other._len:
            return False
        lhs = set(self.cs.values())
        rhs = set(other.cs.values())
        return lhs == rhs

    def entropy(self):
        "un-normalized entropy"
        cs = self.cs
        items = [n for n in cs.values() if n>0]
        return entropy(items)

    def huffman(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()

        if not keys:
            return Node(self) # the empty tree

        # build a tree, start with the leaves:
        nodes = [Node(Multiset({key : cs[key]})) for key in keys]

        while len(nodes) > 1:
            shuffle(nodes)
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
    

def W(item):
    return item.W()


class Node(object):
    "A tree over a multiset"
    "mutable !!"

    def __init__(self, X, left=None, right=None):
        self.X = X
        self._cost = len(X)
        self.left = left
        self.right = right
        assert self.check(), str(self)

    def cost(self):
        return self._cost

    def W(self):
        return self.encode().total_length()

    def check(self):
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            return True
        if not left.X.disjoint(right.X):
            return False
        return X == left.X + right.X and left.check() and right.check()

    def __eq__(self, other):
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            return self.X == other.X
        return self.X == other.X and (
            self.left == other.left and self.right == other.right or
            self.right == other.left and self.left == other.right)

    def __ne__(self, other):
        return not (self==other)
    
    def clone(self):
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            return Node(self.X) # X is immutable.... for now..
        return Node(self.X, left.clone(), right.clone())

    def __getitem__(self, idx):
        if type(idx) is int:
            assert idx==0 or idx==1
            node = [self.left, self.right][idx]
        elif type(idx) is tuple:
            node = self
            for i in idx:
                node = node[i] # recurse
        else:
            raise TypeError
        return node

    def __setitem__(self, idx, node):
        assert isinstance(node, Node), node
        if type(idx) is tuple and len(idx)==1:
            idx = idx[0]
        if type(idx) is int:
            assert idx==0 or idx==1
            assert self.has_children
            child = [self.left, self.right][idx]
            assert node.X == child.X
            if idx==0:
                self.left = node
            else:
                self.right = node
        elif type(idx) is tuple:
            assert len(idx)>1
            child = self
            for i in idx[:-1]:
                child = child[i]
            child[idx[-1]] = node # recurse
        else:
            raise TypeError
        assert self.check()

    @property
    def has_children(self):
        return self.left is not None and self.right is not None

    def all_isos(self, other):
        #if len(self.X) != len(other.X):
        if not self.X.isomorphic(other.X):
            return
        if not self.has_children and not other.has_children:
            yield 1
            return
        elif not self.has_children:
            return
        elif not other.has_children:
            return
        for l_isos in self.left.all_isos(other.left):
          for r_isos in self.right.all_isos(other.right):
            yield 1
        for l_isos in self.right.all_isos(other.left):
          for r_isos in self.left.all_isos(other.right):
            yield 1

    def isomorphic(self, other):
        "up to multiset isomorphism.."
        for iso in self.all_isos(other):
            return True
        return False

    def _subtrees(self):
        yield self
        yield Node(self.X)
        if not self.has_children:
            return
        X = self.X
        left = self.left
        right = self.right
        lsubs = list(left.subtrees())
        rsubs = list(right.subtrees())
        for sub in lsubs:
            yield sub
        for sub in rsubs:
            yield sub
        for l in lsubs:
          for r in rsubs:
            if l.X + r.X == X:
                yield Node(X, l, r)

    def subtrees(self):
        found = set()
        for sub in self._subtrees():
            key = str(sub)
            if key in found:
                continue
            found.add(key)
            yield sub

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
            s = str(X)
            assert s[0]=="{"
            assert s[-1]=="}"
            return "(%s)"%s[1:-1]
        assert left and right
        return "(%s : %s)" % (left, right)
    __repr__ = __str__

    def idxs(self): # dict .keys()
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            yield ()
        else:
            for idx, X in left.idxs():
                yield (0,)+idx
            for idx, X in right.idxs():
                yield (1,)+idx

    def leaves(self): # dict .values()
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            yield X
        else:
            for X in left.leaves():
                yield X
            for X in right.leaves():
                yield X

    def items(self): # dict .items()
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            yield ((), X)
        else:
            for idx, X in left.items():
                yield ((0,)+idx, X)
            for idx, X in right.items():
                yield ((1,)+idx, X)

    def __rmul__(self, r):
        " left multiplication by a Multiset "
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            return Node(r*X)
        return Node(r*X, r*left, r*right)

    def __lmul__(self, r):
        " right multiplication by a Multiset "
        X = self.X
        left = self.left
        right = self.right
        if left is None and right is None:
            return Node(X*r)
        return Node(X*r, left*r, right*r)

    def __mul__(TX, TY):
        if not isinstance(TY, Node):
            return TX.__lmul__(TY)

        X = TX.X
        Y = TY.X
        #print("__mul__", TX, TY)

        XTY = X * TY
        #print("__mul__", XTY)

        #TXY = TX * Y
        top = XTY
        for (idx, r) in TY.items():
            # glue
            top[idx] = TX*r # recurse

        return top


def main():

    X = Multiset({"A":3, "B":1})

    assert (X+X) == Multiset({"A":6, "B":2})
    assert (X+X) == 2*X
    #print(X, X.entropy())

    XX = X*X

    Y = Multiset({"A":2, "B":2})
    #print(Y, Y.entropy())
    assert str(Y) == "{A,A,B,B}"

    A = Multiset({"A" : 1})
    B = Multiset({"B" : 1})
    C = Multiset({"C" : 1})
    D = Multiset({"D" : 1})
    E = Multiset({"E" : 1})
    F = Multiset({"F" : 1})
    G = Multiset({"G" : 1})

    assert A.disjoint(B)
    assert not (A+B).disjoint(B)

    assert (A+2*B).terms() == [A, 2*B]

    # ---------------------------------------------------------------------

    assert Node(A+B, Node(A), Node(B)) == Node(B+A, Node(B), Node(A))

    lhs, rhs = (Node(A+B+C, Node(A+B, Node(A), Node(B)), Node(C)),
        Node(A+B+C, Node(A), Node(B+C, Node(B), Node(C))))

    assert lhs.isomorphic(rhs)

    T = Node(A+B+C, Node(A+B, Node(A), Node(B)), Node(C))
    subs = list(T.subtrees())
    assert len(subs) == 8

    # test left multiplication
    for r in [2, A, A+B]:
        assert r*T == Node(r*A+r*B+r*C, Node(r*A+r*B, Node(r*A), Node(r*B)), Node(r*C))

    T = Node(A+B+C+D, Node(A+B, Node(A), Node(B)), Node(C+D, Node(C), Node(D)))
    subs = list(T.subtrees())
    assert len(subs) == 13, len(subs)

    S, T = Node(A+B, Node(A), Node(B)), Node(B)
    assert S[0] != T
    assert S[1] == T
    U = S.clone()
    assert U==S
    U[1] = T
    assert U[0] != T
    assert U[1] == T
    assert U==S

    T = Node(A+B+C+D, Node(A+B, Node(A), Node(B)), Node(C+D, Node(C), Node(D)))
    assert T[0] == S
    T[0] = Node(A+B)

    T = Node(2*A+B+C+D+E+F, Node(2*A+B+E), Node(C+D+F))
    U = T.clone()
    U[0] = Node(2*A+B+E, Node(2*A), Node(B+E))
    U[0, 1] = Node(B+E, Node(B), Node(E))
    assert U.clone() == U

    T = Node(A+B, Node(A), Node(B))
    S = Node(A+B+2*C, Node(A+B, Node(A), Node(B)), Node(2*C))
    assert str(T*S) == "((((AA) : (BA)) : ((AB) : (BB))) : ((AC,AC) : (BC,BC)))"

    # ---------------------------------------------------------------------

    #print( ((X*Y).entropy(), len(X)*Y.entropy() + len(Y)*X.entropy()))
    assert is_close((X*Y).entropy(), len(X)*Y.entropy() + len(Y)*X.entropy())

    tree = X.huffman()
    assert tree.X == X

    #assert str(tree) == "({B} : {A,A,A})", repr(str(tree))
    #assert str(tree.encode()) == "{0,1,1,1}"

    tree = XX.huffman()
    #assert str(tree.encode()) == "{0,0,0,0,0,0,0,0,0,10,10,10,110,110,110,111}"

    assert XX.W() == 27

    def mkrand(a=1, b=3):
        Z = randint(a, b)*A + randint(a, b)*B + randint(a, b)*C  #+ randint(a, b)*D + randint(a, b)*E
        return Z

    #seed(0)

    for trial in range(1000):
        X = mkrand(1, 5)
        lhs = X.huffman()
        rhs = X.huffman()
        assert lhs.isomorphic(rhs) # huffman is unique up to isomorphism ? this can't be right..

    for trial in range(100):
        X = mkrand(1, 3)
        T = X.huffman()
        for S in T.subtrees():
            assert S.check()

    assert W(Multiset()) == 0

    for trial in range(100):
        a = randint(1, 3)
        b = randint(1, 3)
        c = randint(1, 3)
        X = a*A + b*B + c*C
        lhs = W(X*X)
        rhs = 2*len(X)*W(X)
        assert lhs <= rhs
        #if lhs==rhs: # no nice characterization of this
        #    print(X)
        #else:
        #    print("*")

    for trial in range(100):
        X = mkrand()
        Y = mkrand()
        S = X.huffman()
        T = Y.huffman()
        ST = (X*Y).huffman()
        lhs = W(ST) 
        rhs = len(X)*W(T) + W(S)*len(Y)
        #print(lhs, rhs)
        assert lhs<=rhs
    
    def mkdyadic(a=0, b=4, terms=[A, B, C, D, E]):
        while 1:
            cs = [2**randint(a, b) for t in terms]
            c = sum(cs)
            if bin(c).count('1')==1: # is power of 2
                break
        Z = reduce(add, [c*term for (c, term) in zip(cs, terms)])
        return Z

    for trial in range(100):
        X = mkdyadic()
        Y = mkdyadic()
        #print(X, Y)
        S = X.huffman()
        T = Y.huffman()
        ST = (X*Y).huffman()
        lhs = W(ST) 
        rhs = len(X)*W(T) + W(S)*len(Y)
        #print(lhs, rhs)
        assert lhs==rhs
        assert X.entropy() == W(X)
        assert Y.entropy() == W(Y)
        assert (X*Y).entropy() == lhs

    return
    
    for trial in range(1000):
        a = randint(1, 3)
        b = randint(1, 3)
        c = randint(1, 3)
        X = a*A + b*B + c*C
        lhs = W(X)
        for aa in range(a+1):
         for bb in range(b+1):
          for cc in range(c+1):
            Y = aa*A + bb*B + cc*C
            XY = (a-aa)*A + (b-bb)*B + (c-cc)*C
            assert XY + Y == X
            rhs = W(XY + len(Y)*D) + W(Y)
            assert lhs <= rhs
            if len(Y)==0:
                assert XY == X
                assert XY + len(Y)*D == X
                assert lhs == rhs
    

    return

    for trial in range(100):
        X = mkrand()
        n = randint(2, 5)
        assert n*X.W() == (n*X).W()
        print(X)
        lhs, rhs = n*X.huffman(), (n*X).huffman()
        print(lhs, rhs)
        print()
        #assert n*X.huffman() == (n*X).huffman()
        assert lhs.isomorphic(rhs)
        assert X.huffman().check()
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

        #assert (X*Y) == (Y*X) # nope ( not on the nose.. )
        assert (X*Y).W() == (Y*X).W()
        #assert (X*Y).huffman() == (Y*X).huffman() # nope ( not on the nose.. )

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

    for trial in range(100):

        X0 = mkrand([A, B, C], 1, 3)
        Y = mkrand([D, E, F], 1, 3)

        print(X, Y)
        for a in range(1, 10):
            X = a*X0
            lhs = W(X+Y)
            rhs = W(X + len(Y)*G) + W(Y)
            print(lhs, rhs)
            assert lhs <= rhs
            if lhs==rhs:
                break
        else:
            fail
        print()


    return

    seed(0)

    while 1:
        X = mkrand([A, B, C], 1, 3)
        Y = mkrand([D, E, F], 1, 3)

        print(X, Y)

        a = 1
        while 1:
            print("[%s]"%a, end="", flush=True)
            #aX = a*X
            #lhs, rhs = W(aX*Y), len(aX)*W(Y) + W(aX)*len(Y)
            aY = a*Y
            lhs, rhs = W(X*aY), len(X)*W(aY) + W(X)*len(aY)
            if lhs==rhs:
                break
            print(lhs, rhs)
            assert lhs == a*W(X*Y)
            assert rhs == a*(len(X)*W(Y) + W(X)*len(Y))
            a += 1
            #assert a<10
            if a>10:break
        print(".", end="", flush=True)

    return
        
    found = set()
    #for trial in range(100):
    while 1:
        i = 2**randint(0, 3)
        j = 2**randint(0, 3)
        k = 2**randint(0, 3)
        X = i*A + j*B + k*C
        #X = mkrand([A, B, C], 1, 8)
        lhs, rhs = X.entropy(), X.W()
        if is_close(lhs, rhs):
            #vals = list(X.cs.values())
            vals = [i, j, k]
            vals.sort()
            vals = tuple(vals)
            print(vals)
            if vals not in found:
                print(vals)
                found.add(vals)

    return
        
    
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


