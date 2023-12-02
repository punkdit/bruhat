#!/usr/bin/env python
"""


"""


from __future__ import annotations

from typing import NamedTuple, Union
from snake_egg import EGraph, Rewrite, Var, vars
#import snake_egg as egg

from bruhat.argv import argv



class Mul(NamedTuple):
    l: Expr
    r: Expr
Expr = Union[Mul, Var]

def test_egg():
    a, b, c = vars("a b c")
    rules = [
        Rewrite(Mul(Mul(a, b), c), Mul(a, Mul(b, c)))
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

    assert equ(lhs, rhs) 
    assert not equ(Mul(b, c), Mul(c, b)) 
    assert equ(Mul(b, Mul(a, c)), Mul(Mul(a, b), c)) 



class Morphism(object):
    def __init__(self, solver, expr, tgt=None, src=None):
        self.solver = solver
        self.expr = expr
        self.tgt = tgt
        self.src = src
        self.ref = solver.egraph.add(expr)

    def __str__(self):
        return str(self.expr)

    def __mul__(self, other):
        assert isinstance(other, Morphism)
        assert other.solver is self.solver
        assert other.tgt == self.src
        expr = Mul(self.expr, other.expr)
        return Morphism(self.solver, expr, self.tgt, other.src)

    def __eq__(self, other):
        return self.solver.is_equal(self, other)

    def equate(self, other):
        self.solver.equate(self, other)


class Category(object):
    def __init__(self):
        self.egraph = EGraph()
        a, b, c = vars("a b c")
        self.rules = [ Rewrite(Mul(Mul(a, b), c), Mul(a, Mul(b, c)), "_assoc") ]
        self.cache = {}

    def gen(self, name, tgt=None, src=None):
        g = self.cache.get(name)
        if g is None:
            g = Morphism(self, name, tgt, src)
            self.cache[name] = g
        assert g.tgt == tgt
        assert g.src == src
        return g

    def get(self, expr):
        return self.cache[expr]

    def equate(self, lhs, rhs):
        assert isinstance(lhs, Morphism)
        assert isinstance(rhs, Morphism)
        rule = Rewrite(lhs.expr, rhs.expr, "%s=%s"%(lhs.expr, rhs.expr))
        self.rules.append(rule)

    def mul(self, lhs, rhs):
        assert isinstance(lhs, Morphism)
        assert isinstance(rhs, Morphism)
        expr = Mul(lhs.expr, rhs.expr)
        return self.gen(expr)

    def is_equal(self, lhs, rhs):
        assert isinstance(lhs, Morphism)
        assert isinstance(rhs, Morphism)
        egraph = self.egraph
        egraph.run(self.rules)
        return egraph.equiv(lhs.ref, rhs.ref)



def test_assoc():
    terms = "((a*b)*c)*d (a*(b*c))*d a*((b*c)*d) a*(b*(c*d)) (a*b)*(c*d)".split()

    solver = Category()
    gen = solver.gen
    a, b, c, d = [gen(_) for _ in 'abcd']

    assert a*(b*c) == (a*b)*c
    assert a*(b*c) != (b*a)*c
    for lhs in terms:
      for rhs in terms:
        l = eval(lhs, locals())
        r = eval(rhs, locals())
        assert l == r


def test_inj():
    cat = Category()

    N = 5
    rows = []
    lookup = {}
    for i in range(N):
      row = []
      for j in range(i):
        name = "pop(%d,%d)"%(i, j)
        op = cat.gen(name, i, i-1)
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

    assert pop(3,0)*pop(2,0) == pop(3,1)*pop(2,0)
    assert pop(4,3)*pop(3,0) == pop(4,0)*pop(3,2)
    assert pop(4,1)*pop(3,0)*pop(2,0) == pop(4,1)*pop(3,1)*pop(2,0)

    assert pop(4,1)*pop(3,0) != pop(4,3)*pop(3,1)
    assert pop(4,1)*pop(3,0) != pop(4,3)*pop(3,1)
    assert pop(4,1)*pop(3,0)*pop(2,0) != pop(4,3)*pop(3,1)*pop(2,0)




def test():
    test_egg()
    test_assoc()
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



