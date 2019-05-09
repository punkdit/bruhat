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

    def get_str(self):
        cs = self.cs
        keys = self.keys
        items = [(key if cs[key]==1 else "%d%s"%(cs[key], key)) for key in keys]
        items = '+'.join(items)
        return items

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

    def entropy(self):
        "un-normalized entropy"
        cs = self.cs
        items = [n for n in cs.values() if n>0]
        return entropy(items)

    def huffman(self, sort=False):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()

        if not keys:
            return Node(self) # the empty tree

        # build a tree, start with the leaves:
        nodes = [Node(Multiset({key : cs[key]})) for key in keys]

        while len(nodes) > 1:
            if not sort:
                shuffle(nodes)
            else:
                nodes.sort(key = str)
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

    def W(self):
        return self.encode().total_length()

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

    def __add__(self, other):
        assert self.X.contains(other.X)
        X = self.X
        left = self.left
        right = self.right
        if not self.has_children:
            assert self.X == other.X
            return other
        elif left.X == other.X:
            assert not left.has_children
            left = other
        elif right.X == other.X:
            assert not right.has_children
            right = other
        elif left.X.contains(other.X):
            left = left+other # recurse
        elif right.X.contains(other.X):
            right = right+other # recurse
        else:
            assert 0, (self, other)

        return Node(X, left, right)

    # ------------- rendering ----------------------

    def get_bbox(self, R=1.0):
        "find (width,height) of bounding box for render"
        X = self.X
        left = self.left
        right = self.right
        if not self.has_children:
            s = X.get_str()
            W = (1 + 0.1*len(s)) * R # fudge this
            return (W, R) # square
        lbb = left.get_bbox(R)
        rbb = right.get_bbox(R)
        W = lbb[0] + rbb[0] 
        H = max(lbb[1], rbb[1]) + R
        return (W, H)

    def render(self, x=0, y=0, R=1.0, r=0.2, can=None, name=None):
        "(x, y) is top center of this tree"

        if can is None:
            can = pyx.canvas.canvas()

        X = self.X
        left = self.left
        right = self.right

        can.fill(path.circle(x, y, r), [white])
        can.stroke(path.circle(x, y, r))

        if not self.has_children:
            can.text(x, y-2.4*r, X.get_str(), south)

        else:

            w0, h0 = left.get_bbox(R)
            w1, h1 = right.get_bbox(R)

            w, h = self.get_bbox(R)
            x0 = x-0.5*w+0.5*w0
            x1 = x+0.5*w-0.5*w1
            y0 = y1 = y-R

            # render self...
            can.fill(path.circle(x, y, 0.4*r), [black])
            can.stroke(path.line(x0, y0, x, y), st_Thick)
            can.stroke(path.line(x1, y1, x, y), st_Thick)

            left.render(x0, y0, can=can)
            right.render(x1, y1, can=can)

        if name is not None:
            can.writePDFfile(name)

class Box(object):
    pass

class HBox(Box):
    def __init__(self, items, align="center"):
        self.items = list(items)
        self.align = align

    def render(self, x=0., y=0., can=None, name=None, **kw):
        "(x, y) is top center of this box"
        if can is None:
            can = pyx.canvas.canvas()
        items = self.items
        boxs = [item.get_bbox(**kw) for item in items]
        w = sum(b[0] for b in boxs) # sum widths
        h = max([b[1] for b in boxs]) # max of heights
        x0 = x - 0.5*w
        y0 = y - 0.5*h

        align = self.align
        for i, item in enumerate(self.items):
            b = boxs[i]
            if align == "center":
                item.render(x0 + 0.5*b[0], y0+0.5*b[1], can=can, **kw)
            x0 += b[0]

        if name is not None:
            can.writePDFfile(name)


class TextBox(Box):
    def __init__(self, s, w=1, h=1):
        self.s = s
        self.w = w
        self.h = h

    def get_bbox(self):
        return self.w, self.h

    def render(self, x=0., y=0., can=None, name=None, **kw):
        can.text(x, y-self.h, self.s, south)


def render():
    head = "transversal2018/"
    seed(0)

    a = Multiset({"a" : 1})
    b = Multiset({"b" : 1})
    c = Multiset({"c" : 1})
    d = Multiset({"d" : 1})
    d = Multiset({"d" : 1})
    e = Multiset({"e" : 1})
    f = Multiset({"f" : 1})
    g = Multiset({"g" : 1})
    
    def mkrand(items, a=1, b=3):
        Z = Multiset()
        for A in items:
            Z = Z + randint(a, b)*A
        return Z

    T = Node(a+b+2*c, Node(a+b, Node(a), Node(b)), Node(2*c))
    #T.render(name="pic_a_b_2c.pdf")

    S = ((a+b)*T)
    S.render(name=head+"pic_left_a_b_2c.pdf")
    #T = T*T

    U = Node(a*a+b*a, Node(a*a), Node(b*a))
    SU = S+U
    box = HBox([S, TextBox("$+$"), U, TextBox("$=$"), SU])
    box.render(name=head+"pic_add.pdf")
    
    S = Node(a+b, Node(a), Node(b))
    #(S*T).render(name=head+"pic_prod.pdf")
    
    box = HBox([S, TextBox(r"$\times$"), T, TextBox("$=$"), S*T])
    box.render(name=head+"pic_prod.pdf")

    X = mkrand([a,b,c,d,e,f,g])
    TX = X.huffman(sort=True)
    print(W(TX))
    box = HBox([TX])
    box.render(name=head+"pic_huffman.pdf")

    X = a + b + 2*c + 4*d
    TX = X.huffman(sort=True)
    box = HBox([TX])
    box.render(name=head+"pic_dyadic.pdf")

    print("OK")


def randtree(X):
    "random complete tree on X"
    nodes = [Node(Y) for Y in X.terms()]
    while len(nodes)>1:
        left = nodes.pop(randint(0, len(nodes)-1))
        right = nodes.pop(randint(0, len(nodes)-1))
        node = Node(left.X+right.X, left, right)
        nodes.append(node)
    return nodes[0]


def main():

    X = Multiset({"a":3, "b":1})

    assert (X+X) == Multiset({"a":6, "b":2})
    assert (X+X) == 2*X
    #print(X, X.entropy())

    XX = X*X

    Y = Multiset({"a":2, "b":2})
    #print(Y, Y.entropy())
    assert str(Y) == "{a,a,b,b}"

    A = Multiset({"a" : 1})
    B = Multiset({"b" : 1})
    C = Multiset({"c" : 1})
    D = Multiset({"d" : 1})
    E = Multiset({"e" : 1})
    F = Multiset({"f" : 1})
    G = Multiset({"g" : 1})

    assert A.disjoint(B)
    assert not (A+B).disjoint(B)

    assert (A+2*B).terms() == [A, 2*B]

    assert not A.contains(B)
    assert (A+B).contains(B)
    assert not (A+B).contains(2*B)

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
    assert str(T*S) == "((((aa) : (ba)) : ((ab) : (bb))) : ((ac,ac) : (bc,bc)))"

    def randmultiset(a=0, b=4):
        Z = randint(a, b)*A + randint(a, b)*B + randint(a, b)*C + randint(a, b)*D + randint(a, b)*E
        return Z

    for trial in range(100):
        X = randmultiset()
        TX = randtree(X)
        assert X.entropy() <= W(X) <= W(TX)
        
        Y = randmultiset()
        TY = randtree(Y)

        # W is a derivation on complete trees
        assert W(TX*TY) == len(X)*W(TY) + W(TX)*len(Y)

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

    if argv.render:
        # yay globals...
        import pyx
        from pyx import path, deco, trafo, style, text, color, deformer
        from pyx.color import rgb, cmyk
        from pyx.color import rgbfromhexstring as rgbhex
    
        black = rgb(0., 0., 0.)
        blue = rgb(0., 0., 0.8)
        lred = rgb(1., 0.4, 0.4)
        red = rgb(1., 0.0, 0.0)
        white = rgb(1., 1., 1.)

        grey = rgb(0.75, 0.75, 0.75)
        shade = grey
        shade0 = rgb(0.25, 0.25, 0.25)
        shade1 = rgb(0.80, 0.80, 0.80)
        shade2 = rgb(0.85, 0.85, 0.85)
        
        light_shade = rgb(0.85, 0.65, 0.1)
        light_shade = rgb(0.9, 0.75, 0.4)
        
        
        north = [text.halign.boxcenter, text.valign.top]
        northeast = [text.halign.boxright, text.valign.top]
        northwest = [text.halign.boxleft, text.valign.top]
        south = [text.halign.boxcenter, text.valign.bottom]
        southeast = [text.halign.boxright, text.valign.bottom]
        southwest = [text.halign.boxleft, text.valign.bottom]
        east = [text.halign.boxright, text.valign.middle]
        west = [text.halign.boxleft, text.valign.middle]
        center = [text.halign.boxcenter, text.valign.middle]
        
        
        st_dashed = [style.linestyle.dashed]
        st_dotted = [style.linestyle.dotted]
        st_round = [style.linecap.round]
        #st_mitre = [style.linecap.square]
        
        st_thick = [style.linewidth.thick]
        st_Thick = [style.linewidth.Thick]
        st_THick = [style.linewidth.THick]
        st_THIck = [style.linewidth.THIck]
        st_THICk = [style.linewidth.THICk]
        st_THICK = [style.linewidth.THICK]

        render()
    else:
        main()


