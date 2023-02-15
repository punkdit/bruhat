#!/usr/bin/env python3

"""
Implement Thompson's group F

https://conf.math.illinois.edu/~mobeidin/thompson.pdf
"""

from random import choice, seed
from functools import total_ordering

from bruhat.action import mulclose
from bruhat import unify

try:
    from huygens import config
    config(text="pdflatex", latex_header=r"""
    \usepackage{amsmath}
    \usepackage{amssymb}
    """)
    from huygens.namespace import *
except ImportError:
    print("could not find huygens package")


@total_ordering
class Dyadic(object):
    """
        Dyadic rational 
        https://en.wikipedia.org/wiki/Dyadic_rational
    """

    def __init__(self, top, exp=0):
        while top%2==0 and exp>0:
            #print(top, exp)
            top //= 2
            exp -= 1
        #print(top, exp)
        self.top = int(top)
        self.exp = int(exp)
        # top / 2^exp

    def __str__(self):
        return "%s/%s"%(self.top, 2**self.exp)

    @classmethod
    def promote(cls, item):
        if isinstance(item, Dyadic):
            return item
        return Dyadic(item)

    def __eq__(self, other):
        other = self.promote(other)
        return (self.top, self.exp) == (other.top, other.exp)

    def __lt__(self, other):
        other = self.promote(other)
        # self < other
        a = other - self
        return a.top > 0

    def __add__(self, other):
        other = self.promote(other)
        a, b = self.top, self.exp
        c, d = other.top, other.exp
        top = 2**(d-min(b,d)) * a + 2**(b-min(b,d)) * c
        exp = max(b, d)
        return Dyadic(top, exp)

    def __sub__(self, other):
        other = self.promote(other)
        a, b = self.top, self.exp
        c, d = other.top, other.exp
        top = 2**(d-min(b,d)) * a - 2**(b-min(b,d)) * c
        exp = max(b, d)
        return Dyadic(top, exp)

    def __neg__(self):
        return Dyadic(-self.top, self.exp)

    def __mul__(self, other):
        other = self.promote(other)
        a, b = self.top, self.exp
        c, d = other.top, other.exp
        return Dyadic(a*c, b+d)

    def __truediv__(self, other):
        other = self.promote(other)
        top = other.top
        s = bin(top)[2:]
        if s.count('1') != 1:
            return NotImplemented
        assert top > 0, "todo..."
        exp = len(s) - 1
        top = self.top * (2**other.exp)
        exp = self.exp + exp
        return Dyadic(top, exp)


class Interval(object):
    def __init__(self, a, b):
        assert a<b, "not an interval: %s>=%s"%(a, b)
        a = Dyadic.promote(a)
        b = Dyadic.promote(b)
        self.a = a
        self.b = b

    def __str__(self):
        return "[%s, %s]"%(self.a, self.b)

    def __eq__(self, other):
        return self.a==other.a and self.b==other.b


#class Key(object):
#    def __init__(self, idx=()):
#        self.idx = idx
#    def __eq__(self, other):
#        if type(other) is tuple:
#            other = Key(other)
#        return self.idx == other.idx
#    def __str__(self):
#        return str(self.idx)
#    __repr__ = __str__
#    def __getitem__(self, i):
#        return self.idx[i]
#    def __len__(self):
#        return len(self.idx)
#    def __add__(self, other):
#        if type(other) is tuple:
#            other = Key(other)
#        return Key(self.idx+other.idx)
#    def __lt__(lhs, rhs):
#        # depth first search
#        if lhs.idx[:len(rhs.idx)] == rhs.idx: # idx startswith jdx
#            return True
#        return lhs.idx < rhs.idx

#
#def join(lhs, rhs):
#    lkeys = lhs.keys().sort()
#    rkeys = rhs.keys().sort()
#    li = ri = 0
#    while 1:
#        if li < len(lkeys) and ri < len(rkeys):
#            lk = lkeys[li]
#            rk = rkeys[ri]
#            lexpr = lhs[lk]
#            rexpr = rhs[rk]
#            # fark !!
#        elif li < len(lkeys):
#            lk = lkeys[li]
#        elif ri < len(rkeys):
#            rk = rkeys[ri]
#        else:
#            break
#            

class Tree(object):
    tp = unify.Theory("T")
    tp.add_op("*", 2)

    def __mul__(self, other):
        assert isinstance(other, Tree), other
        return Pair(self, other)
    def to_expr(self, prefix="v", tp=None):
        tp = self.tp
        ks = self.keys()
        #ks.sort() # depth first search
        send = {}
        count = 0
        for k in ks:
            tree = self[k]
            if isinstance(tree, Leaf):
                send[k] = tp.get_var("%s_%s"%(prefix,count))
                count += 1
            elif isinstance(tree, Pair):
                a, b = tree.a, tree.b
                expr = send[k+(0,)] * send[k+(1,)]
                send[k] = expr
            else:
                assert 0
        return send[()]
    @classmethod
    def from_expr(cls, expr):
        vs = expr.all_vars()
        ns = dict((str(v), Leaf()) for v in vs)
        tree = eval(str(expr), ns)
        return tree
    def find_carets(self):
        items = {}
        count = 0
        for k,tree in self.items():
            tree = self[k]
            if type(tree) == Leaf:
                count += 1
            elif tree.is_caret():
                assert count>1
                items[count] = k
        return items
    def log_render(self, dx=1, dy=-1, depth=None):
        if depth is None:
            depth = self.get_depth()
        cvs = Canvas()
        mx = 4
        cvs.stroke(path.circle(0, 0, 0.1), [red]) # origin
        for key, value in self.items():
            if isinstance(value, Leaf):
                continue
            W = 2**(depth-1)*dx
            H = 2**(depth-1)*dy
            x0, y0 = 0, 0 # find the root position of this Pair
            for k in key:
                y0 += H
                W //= 2
                x0 += W*(k-1/2)
                H //= 2
            m = 0
            for i in range(depth-len(key)):
                m += H
                H //= 2
            cvs.stroke(path.line( mx*x0, y0, mx*(x0-m*(1/4)), y0+m ), st_THIck)
            cvs.stroke(path.line( mx*x0, y0, mx*(x0+m*(1/4)), y0+m ), st_THIck)
        #print(cvs.get_bound_box())
        #print(depth)
        # height is 2**depth - 1
        return cvs

    def render(self, w=1, h=-1):
        n = self.nleaves()
        layout = {} # coordinates
        height = (n-1)/2
        items = list(self.items())
        i = 0
        count = 0
        origin = None # origin of canvas

        # First we layout the leaves.
        while i < len(items):
            key, tree = items[i]
            if isinstance(tree, Pair):
                i += 1
                continue
            items.pop(i)
            layout[key] = (count*w, height*h)
            if origin is None:
                origin = layout[key]
            count += 1
        assert len(layout) == n

        # Now we layout all the internal nodes of the tree,
        # working backwards from the leaves.
        m = h/w # slope of our lines
        while len(items):
            i = 0
            while i < len(items):
                key, tree = items[i]
                assert isinstance(tree, Pair)
                ka, kb = key+(0,), key+(1,)
                if ka not in layout or kb not in layout:
                    i += 1
                    continue # not ready
                items.pop(i)
                xa, ya = layout[ka]
                xb, yb = layout[kb]
                assert xb>xa
                # Yes i actually did some algebra:
                dx = xb-xa
                dy = yb-ya
                b = (1/2)*(m*dx - dy)
                c = dy/m
                d = (1/2)*(dx-c)
                x = xa+c+d
                y = ya+dy+b
                assert tree not in layout, tree
                #print(tree, "at", x, y)
                layout[key] = x, y

        cvs = Canvas()
        x0, y0 = origin
        #cvs.stroke(path.circle(0,0,0.2), [red]) # show origin
        cvs.append(Translate(-x0, -y0))
        for key in layout:
            if isinstance(self[key], Leaf):
                continue
            x0, y0 = layout[key]
            ka, kb = key+(0,), key+(1,)
            for (x1, y1) in [layout[ka], layout[kb]]:
                cvs.stroke(path.line(x0, y0, x1, y1), st_THIck+st_round)
        return cvs


def all_trees(n_or_items):
    if type(n_or_items) is int:
        items = [Leaf() for i in range(n_or_items)]
    else:
        items = list(n_or_items)
    n = len(items)
    if n==1:
        yield items[0]
        return
    elif n==2:
        yield Pair(*items)
        return
    elif n==0:
        return
    found = set()
    for i in range(n-1):
        a, b = items[i:i+2]
        trees = items[:i] + [Pair(a, b)] + items[i+2:]
        for tree in all_trees(trees):
            if tree not in found:
                yield tree
                found.add(tree)



class Leaf(Tree):
    def __init__(self, name="_"):
        self.name = name
    def __str__(self):
        return self.name
    def __eq__(self, other):
        return isinstance(other, Leaf)
    def __hash__(self):
        return 0
    def copy(self):
        return Leaf(self.name)
    def __getitem__(self, key):
        #assert type(key) is Key, key
        if key==():
            return self
        raise IndexError
    def keys(self):
        #return [Key()]
        return [()]
    def items(self):
        yield (), self
    def sup(lhs, rhs):
        if isinstance(rhs, Leaf):
            return lhs.copy()
        return rhs.copy()
#    def join(lhs, rhs, send={}, idx=0):
#        send[0] = rhs
#        idx += 1
#        return send, idx
    def nleaves(self):
        return 1
    def is_caret(self):
        return False
    def rewrite(self, key, tree, depth=0):
        #print("  "*depth+"send", self, key, tree)
        assert key==(), "%s --> %s"%(key, tree)
        return tree
    def internal_nodes(self, root=True):
        return
        yield
    def reassoc(self, key):
        assert 0
    def get_depth(self):
        return 0

class Pair(Tree):
    def __init__(self, a, b):
        assert isinstance(a, Tree), a
        assert isinstance(b, Tree), b
        self.a = a.copy()
        self.b = b.copy()
    def __str__(self):
        return "(%s*%s)"%(self.a, self.b)
    def __eq__(self, other):
        return isinstance(other,Pair) and self.a==other.a and self.b==other.b
    def __hash__(self):
        a = self.a.__hash__()
        b = self.b.__hash__()
        return hash( (a, b) )
    def copy(self):
        return Pair(self.a, self.b)
    def __getitem__(self, key):
        #assert type(key) is Key, key
        if key==():
            return self
        head, tail = key[0], key[1:]
        return [self.a, self.b][head][tail]
    def keys(self):
        #ks = [Key()]
        ks = []
        for key in self.a.keys():
            #ks.append(Key((0,))+key)
            ks.append((0,)+key)
        for key in self.b.keys():
            #ks.append(Key((1,))+key)
            ks.append((1,)+key)
        ks.append(())
        return ks
    def items(self):
        for (k,v) in self.a.items():
            #assert self.a[k] == v
            yield (0,)+k, v
        for (k,v) in self.b.items():
            #assert self.b[k] == v
            yield (1,)+k, v
        yield (), self
    def sup(lhs, rhs):
        if isinstance(rhs, Leaf):
            return lhs.copy()
        a = lhs.a.sup(rhs.a)
        b = lhs.b.sup(rhs.b)
        return Pair(a, b)
#    def join(lhs, rhs, send={}, idx=0):
#        #la, lb = lhs.a, lhs.b
#        #ra, rb = rhs.a, rhs.b
#        #???
#        return send, idx
    def nleaves(self):
        return self.a.nleaves() + self.b.nleaves()
    def is_caret(self):
        return self.a.__class__ is Leaf and self.b.__class__ is Leaf
    def rewrite(self, key, tree, depth=0):
        #print("  "*depth+"send", self, key, tree)
        assert self[key]
        if key==():
            return tree
        head, tail = key[0], key[1:]
        children = [self.a, self.b]
        children[head] = children[head].rewrite(tail, tree, depth+1)
        return Pair(*children)
    def internal_nodes(self, root=True):
        if type(self.a) is Pair:
            for k in self.a.internal_nodes(False):
                yield (0,)+k
        if type(self.b) is Pair:
            for k in self.b.internal_nodes(False):
                yield (1,)+k
        if not root:
            yield ()
    def reassoc(self, key):
        assert len(key)
        child = self[key]
        tail, idx = key[:-1], key[-1]
        parent = self[tail]
        if idx == 0:
            a = child[0,]
            b = child[1,]
            c = parent[(1-idx,)]
            tree = Pair(a, Pair(b, c))
        else:
            a = parent[(1-idx,)]
            b = child[0,]
            c = child[1,]
            tree = Pair(Pair(a, b), c)
        tree = self.rewrite(tail, tree)
        return tree
    def get_depth(self):
        return 1 + max(self.a.get_depth(), self.b.get_depth())


def cancel_carets(lhs, rhs):
    while 1:
        #print("cancel_carets", lhs, rhs)
        l = lhs.find_carets()
        r = rhs.find_carets()
        #print("\t", l, r)
        items = set(l.keys()) & set(r.keys())
        if not items:
            break
        i = min(items) # deterministic
        lhs = lhs.rewrite( l[i], Leaf() ) # remove caret
        rhs = rhs.rewrite( r[i], Leaf() ) # remove caret
    return lhs, rhs
    

class Assoc(object):
    "Element of Thompson group F"
    def __init__(self, tgt, src):
        assert isinstance(tgt, Tree)
        assert isinstance(src, Tree)
        assert tgt.nleaves() == src.nleaves()
        tgt, src = cancel_carets(tgt, src)
        self.tgt = tgt
        self.src = src
    def __str__(self):
        return "%s<--%s"%(self.tgt, self.src)
    def __eq__(self, other):
        return self.tgt == other.tgt and self.src == other.src
    def __hash__(self):
        return hash( (self.tgt.__hash__(), self.src.__hash__() ) )
    def __mul__(l, r):
        # l.tgt <--- l.src == r.tgt <--- r.src
        ltgt = l.tgt.to_expr("l")
        lsrc = l.src.to_expr("l")
        rtgt = r.tgt.to_expr("r")
        rsrc = r.src.to_expr("r")
        #print("Assoc.__mul__", ltgt, "<---", lsrc, "===", rtgt, "<---", rsrc)
    
        send = unify.Expr.unify(lsrc, rtgt)
        #print(send)
        tgt, src = ltgt.subs(send), rsrc.subs(send)
        #print("Assoc.__mul__", tgt, "<---", src)
        tgt = Tree.from_expr(tgt)
        src = Tree.from_expr(src)
        #print("Assoc.__mul__", tgt, "<---", src)
        return Assoc(tgt, src)
    def __invert__(self):
        return Assoc(self.src, self.tgt)
    def render(self):
        top = self.tgt.render(h=-1)
        bot = self.src.render(h=+1)
        cvs = Canvas([top, bot])
        return cvs
    


def bracket(g, h):
    return (~g)*(~h)*g*h


def test():

    a = Dyadic(2)
    b = a/4
    assert b == Dyadic(2, 2)
    assert b == Dyadic(1, 1)
    assert str(b) == "1/2"
    assert Dyadic(3*256, 4) == 3*256//(2**4)

    a = Dyadic(7)
    b = Dyadic(1)/2
    assert a+b == Dyadic(15)/2
    assert a-b == Dyadic(13)/2
    assert a*b == Dyadic(7)/2

    assert b < a
    assert b <= a
    assert -b > -a

    assert str(Interval(b, a)) == "[1/2, 7/1]"


def test_tree():
    x = Leaf()
    lhs = x*(x*x)
    rhs = (x*x)*x
    #print(lhs, rhs)
    #for key in lhs.keys():
    #    print(key, lhs[key])
    #print(lhs.sup(rhs))
    lhs, rhs = x*((x*x)*x), (x*x)*x
    #print(lhs.sup(rhs))

    lhs, rhs = x*((x*x)*x), (x*x)*x
    ks = lhs.keys()
    assert ks == [(0,), (1, 0, 0), (1, 0, 1), (1, 0), (1, 1), (1,), ()]
    assert lhs.find_carets() == {3: (1,0)}

    s = ((x*x)*x)*x
    t = (x*x)*(x*x)

    assert t.find_carets() == {2: (0,), 4: (1,)}
    assert s.find_carets() == {2: (0, 0)}
    t, s = cancel_carets(t, s)
    assert t == x*(x*x), t
    assert s == (x*x)*x, s

    assert s.reassoc((0,)) == t
    assert t.reassoc((1,)) == s

    tree = s*s
    assert list(tree.internal_nodes()) == [(0, 0), (0,), (1, 0), (1,)]

    catalan = [1, 1, 2, 5, 14, 42, 132]
    for i in range(7):
        n = i+1
        trees = list(all_trees(n))
        assert len(trees) == catalan[i]
        for tree in trees:
            assert tree.nleaves() == n
            keys = list(tree.internal_nodes())
            assert len(keys) == max(0, n-2), (n, len(keys))
            for k in keys:
                other = tree.reassoc(k)
                assert other != tree
                assert other in trees

    I = Assoc( x, x )
    assert I == I
    assert I*I == I
    F = Assoc( (x*x)*x, (x*x)*x )
    assert F == I

    A = Assoc( (x*x)*x, x*(x*x) )
    B = Assoc( x*((x*x)*x), x*(x*(x*x)) )
    #print(A, B)

    assert A*A == Assoc( (((x*x)*x)*x), (x*(x*(x*x))) )

    assert I*A == A
    assert A*I == A
    assert A*~A == I
    assert ~A*A == I

    lhs = bracket( A*~B, ~A*B*A )
    rhs = bracket( A*~B, (~A)*(~A)*B*A*A )
    assert lhs == rhs

    G = mulclose([A, B, ~A, ~B], maxsize=100)
    assert len(G) >= 100

    #for g in G:
    #    print(g.diamond_str())

    for n in range(1, 6):
        trees = list(all_trees(n))
        found = []
        for src in trees:
            for k in src.internal_nodes():
                tgt = src.reassoc(k)
                f = Assoc(tgt, src)
                found.append(f)
        uniq = set(found)
        #print(len(trees), len(found), len(set(found)))
        assert len(uniq) == max(0, 2**(n-1) - 2) # woah !
        uniq = list(uniq)
        uniq.sort(key = str)
        for f in uniq:
            assert ~f in uniq
            #print(f)
        #print()

#def make_eq(l_cvs, r_cvs, space=2.0):
#    cvs = Canvas()
#    l_bb = l_cvs.get_bound_box()
#    r_bb = r_cvs.get_bound_box()
#    lw = max(l_bb.width, space)
#    rw = max(r_bb.width, space)
#    cvs.insert(-0.7*lw, 0., l_cvs)
#    cvs.insert(+0.7*rw, 0., r_cvs)
#    bb = cvs.get_bound_box()
#    x, y = bb.center
#    cvs.text(x, y, r"$=$", [Scale(2.)]+st_center)
#    return cvs

def make_eq(items):
    cvs = Canvas()
    x = 0
    for item in items:
        if type(item) is str:
            item = Canvas().text(0, 0, item, [Scale(3.)]+st_center)
        cvs.insert(x, 0, item)
        x += 1.1*item.get_bound_box().width + 0.4
    return cvs
    

def test_render():
    x = Leaf()
    t = ((x*x)*x)*x
    s = ((x*x)*(x*x))

    #cvs = s.render(h=-1)
    #cvs.writePDFfile("tree.pdf")
    #return
    A = Assoc( (x*x)*x, x*(x*x) )
    B = Assoc( x*((x*x)*x), x*(x*(x*x)) )
    G = mulclose([A, B, ~A, ~B], maxsize=60)
    G = list(G)

    seed(2)

    question = Canvas().text(0, 0, "?", [Scale(7.)]+st_center)
    cvs = Canvas()
    x, y = 0, 0
    for i in range(4):
        f = choice(G)
        g = choice(G)
    
        row = [f.render(), r"$\times$", g.render(), r"$=$", (f*g).render()]
        #if i==0:
        #    row[-1] = question
        row = make_eq(row)
        cvs.insert(x, y, row)
        y += 1.1*row.get_bound_box().height

    bb = cvs.get_bound_box()
    bg = Canvas()
    m = 0.1
    bg.fill(path.rect(bb.llx-m, bb.lly-m, bb.width+2*m, bb.height+2*m), [white])
    cvs = bg.append(cvs)
    if 0:
        cvs.writePDFfile("tree-puzzle.pdf")
        cvs.writePNGfile("tree-puzzle.png")
    else:
        cvs.writePDFfile("tree-soln.pdf")
        cvs.writePNGfile("tree-soln.png")


if __name__ == "__main__":
    test()
    test_tree()
    test_render()
    print("OK\n")


