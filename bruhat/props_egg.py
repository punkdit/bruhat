#!/usr/bin/env python
"""

pip3 install snake-egg

"""


from __future__ import annotations

from typing import NamedTuple, Union
from snake_egg import EGraph, Rewrite, Var, vars
#import snake_egg as egg

from bruhat.argv import argv


class Ident(NamedTuple):
    def __str__(self):
        return "I"

class Mul(NamedTuple):
    l: Expr
    r: Expr
    def __str__(self):
        return "(%s*%s)"%self

Expr = Union[Mul, Ident, Var]


def test_egg():
    a, b, c = vars("a b c")
    rules = [
        Rewrite(Mul(Mul(a, b), c), Mul(a, Mul(b, c)), "_assoc"),
        Rewrite(Mul(Ident(), a), a, "_left_ident"),
        Rewrite(Mul(a, Ident()), a, "_right_ident"),
    ]

    egraph = EGraph()

    a, b, c, d = "abcd"
    rules.append( Rewrite( Mul(a, b), Mul(b, a), "com" ) )

    def equ(lhs, rhs):
        l = egraph.add(lhs)
        r = egraph.add(rhs)
        egraph.run(rules)
        return egraph.equiv(l, r)

    lhs = Mul(c, Mul(b, Mul(a,a)))
    rhs = Mul(Mul(c, b), Mul(a,a))
    i = Ident()

    assert equ(lhs, rhs) 
    assert not equ(Mul(b, c), Mul(c, b)) 
    assert equ(Mul(b, Mul(a, c)), Mul(Mul(a, b), c)) 
    assert equ(Mul(i, a), a)
    assert not equ(Mul(b, i), a)
    assert equ(Mul(Mul(b,c), i), Mul(b,c))



class Morphism(object):
    def __init__(self, cat, expr, tgt=None, src=None):
        self.cat = cat
        self.expr = expr
        self.tgt = tgt
        self.src = src
        self.ref = cat.egraph.add(expr)
        #print("Morphism", expr, tgt, src)

    def __str__(self):
        return str(self.expr)

    def __mul__(self, other):
        assert isinstance(other, Morphism)
        assert other.cat is self.cat
        return self.cat.mul(self, other)
        #assert other.tgt == self.src
        #expr = Mul(self.expr, other.expr)
        #return Morphism(self.cat, expr, self.tgt, other.src)

    def __eq__(self, other):
        return self.cat.is_equal(self, other)

    def equate(self, other):
        self.cat.equate(self, other)


class Category(object):
    def __init__(self):
        self.egraph = EGraph()
        a, b, c = vars("a b c")
        self.rules = [
            Rewrite(Mul(Mul(a, b), c), Mul(a, Mul(b, c)), "_assoc"),
            Rewrite(Mul(Ident(), a), a, "_left_ident"),
            Rewrite(Mul(a, Ident()), a, "_right_ident"),
        ]
        self.cache = {}
        self.update = True

    def morphism(self, expr, tgt=None, src=None):
        key = (expr, tgt, src)
        g = self.cache.get(key)
        if g is None:
            g = Morphism(self, expr, tgt, src)
            self.cache[key] = g
            self.update = True
        assert g.tgt == tgt
        assert g.src == src
        return g

    def identity(self, tgt=None):
        return self.morphism(Ident(), tgt, tgt)

#    def get(self, expr):
#        return self.cache[expr]

    def equate(self, lhs, rhs):
        assert isinstance(lhs, Morphism)
        assert isinstance(rhs, Morphism)
        rule = Rewrite(lhs.expr, rhs.expr, "%s=%s"%(lhs.expr, rhs.expr))
        self.rules.append(rule)
        self.update = True

    def mul(self, lhs, rhs):
        assert isinstance(lhs, Morphism)
        assert isinstance(rhs, Morphism)
        assert lhs.src == rhs.tgt, (
            "can't compose %s*%s: %s!=%s"%(lhs, rhs, lhs.src, rhs.tgt))
        expr = Mul(lhs.expr, rhs.expr)
        return self.morphism(expr, lhs.tgt, rhs.src)

    def is_equal(self, lhs, rhs):
        assert isinstance(lhs, Morphism)
        assert isinstance(rhs, Morphism)
        egraph = self.egraph
        if self.update:
            egraph.run(self.rules)
            self.update = False
        return egraph.equiv(lhs.ref, rhs.ref)



def test_assoc():
    terms = "((a*b)*c)*d (a*(b*c))*d a*((b*c)*d) a*(b*(c*d)) (a*b)*(c*d)".split()

    cat = Category()
    morphism = cat.morphism
    identity = cat.identity
    a, b, c, d = [morphism(_) for _ in 'abcd']

    assert a*a is a*a
    assert a*(b*c) == (a*b)*c
    assert a*(b*c) != (b*a)*c
    for lhs in terms:
      for rhs in terms:
        l = eval(lhs, locals())
        r = eval(rhs, locals())
        assert l == r

    assert a*identity(a.src) == a
    assert (b*a)*identity(b.tgt) == b*a

    f = morphism("f", 1, 0)
    assert f*identity(f.src) == f
    assert identity(f.tgt)*f == f


def test_group():
    cat = Category()
    n = 4
    gens = [cat.morphism(c) for c in "abcdefgh"[:n]]
    I = cat.identity()
    for g in gens:
        (g*g).equate(I)
    normal = []
    for i in range(n-1):
        a, b = gens[i:i+2]
        normal.append( a*b*a*b*a*b )
        (a*b*a).equate(b*a*b)
        #(a*b*a*b).equate(b*a)
        #(a*b*a*b*a).equate(b)
        #(a*b*a*b*a*b).equate(I)
        #(b*a*b*a).equate(a*b)
        #(b*a*b*a*b).equate(a)
        #(b*a*b*a*b*a).equate(I)
    for i in range(n):
      for j in range(i+2, n):
        a, b = gens[i], gens[j]
        normal.append( a*b*a*b )
        (a*b).equate(b*a)

    #G = list(gens)
    G = [I]
    bdy = [I]
    while bdy:
        N = len(G)
        print(len(G), len(bdy))
        _bdy = []
        #for g in list(G):
        for g in bdy:
          for h in gens:
            gh = g*h
            if gh not in G:
                print("[%s]"%len(G), end="", flush=True)
                G.append(gh)
                _bdy.append(gh)
                #assert len(G) < 50
                for k in normal:
                    (gh*k).equate(k*gh)
            else:
                print("/", end="", flush=True)
            #for i,a in enumerate(G):
            #  for b in G[i+1:]:
            #    assert a!=b, "%s = %s"%(a,b)
            #if len(G) == 49:
            #    print()
            #    for g in G:
            #        print(g)
        print()
        #if len(G)==N:
        #    break
        bdy = _bdy
    print(len(G))


def test_inj():
    cat = Category()

    N = 5
    rows = []
    lookup = {}
    for i in range(N):
      row = []
      for j in range(i):
        name = "pop(%d,%d)"%(i, j)
        op = cat.morphism(name, i, i-1)
        lookup[i,j] = op
        row.append(op)
      rows.append(row)

    for i in range(N):
      for j in range(i):
        op = rows[i][j]
        print(op, end=" ")
      print()

    pop = lambda i,j: lookup[i,j]

    for tgt in range(1, N):
      for i in range(tgt):
       for j in range(tgt-1):
        e = lookup[tgt, i] * lookup[tgt-1, j]

    for src in range(2, N):
     for i in range(src):
      for j in range(src-1):
        lhs = pop(src, i) * pop(src-1, j)
        model = list(range(src))
        model.pop(i)
        model.pop(j)
        nodel = list(range(src))
        if i<j:
            rhs = pop(src, j+1) * pop(src-1, i)
            nodel.pop(j+1)
            nodel.pop(i)
            assert model == nodel, (src, i, j)
        elif i>j:
            rhs = pop(src, j) * pop(src-1, i-1)
            nodel.pop(j)
            nodel.pop(i-1)
            assert model == nodel, (src, i, j)
        else:
            assert i==j
            rhs = pop(src, j+1) * pop(src-1, i)
            nodel.pop(j+1)
            nodel.pop(i)
            assert model == nodel, (src, i, j)
        #print(lhs, "==", rhs)
        lhs.equate(rhs)

    #return

    if 0:
      for src in range(2, N):
        interp = []
        for i in range(src):
         for j in range(src-1):
           model = list(range(src))
           model.pop(i)
           model.pop(j)
           g = pop(src, i) * pop(src-1, j)
           interp.append((g, model))
        for lhs in interp:
          for rhs in interp:
            if lhs[1] == rhs[1]:
                lhs[0].equate(rhs[0])

    print(pop(4,1)*pop(3,0)*pop(2,0) )
    assert pop(3,0)*pop(2,0) == pop(3,1)*pop(2,0)
    assert pop(4,3)*pop(3,0) == pop(4,0)*pop(3,2)
    assert pop(4,1)*pop(3,0)*pop(2,0) == pop(4,1)*pop(3,1)*pop(2,0)

    assert pop(4,1)*pop(3,0) != pop(4,3)*pop(3,1)
    assert pop(4,1)*pop(3,0) != pop(4,3)*pop(3,1)
    assert pop(4,1)*pop(3,0)*pop(2,0) != pop(4,3)*pop(3,1)*pop(2,0)




def test():
    test_egg()
    test_assoc()
    test_group()
    test_inj()


if __name__ == "__main__":

    from time import sleep, time
    start_time = time()
    profile = argv.profile
    name = argv.next() or "test"
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    else:
        fn = eval(name)
        fn()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))



