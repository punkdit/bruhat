#!/usr/bin/env python

import sys
from time import time
start_time = time()

from bruhat.theory_z3 import Sort, Expr, Const, Operator, Theory, Z3Solver

# https://ncatlab.org/nlab/show/globular+set
BARS = ["-", "=", "≡", "≣"] # ...?
OPS = ["*", "<<", "@", "??"]

class NotComposable(Exception):
    pass

INDENT = "  "

def build_cell(self, tgt, src, name=""):
    # this is a bit clunky... a covariant inheritance problem
    assert (tgt is None) == (src is None), (tgt, src)
    if tgt is None:
        dim = 0
        desc = name
    else:
        assert isinstance(tgt, Expr), tgt
        assert isinstance(src, Expr), src
        dim = src.dim+1
        assert src.dim == tgt.dim
        #desc = "(%s <%s%s%s %s)"%(tgt.name, BARS[dim-1], name, BARS[dim-1], src.name)
        desc = "(%s <%s%s%s %s)"%(tgt, BARS[dim-1], name, BARS[dim-1], src)
        # check globular conditions
        assert src.src is tgt.src
        assert src.tgt is tgt.tgt
    self.dim = dim
    self.src = src
    self.tgt = tgt
    self.desc = desc
    theory = self.theory


class Cell(Const):
    def __init__(self, theory, name, sort, tgt=None, src=None):
        Const.__init__(self, theory, name, sort)
        build_cell(self, tgt, src, name)

    def __str__(self):
        return self.desc

    def __repr__(self):
        if self.tgt is not None:
            return "Cell(%s, %s)"%(self.tgt, self.src)
        else:
            return "Cell(%r)"%(self.name,)

    def _deepstr(self, depth=0):
        if self.tgt is None:
            lines = [INDENT*depth + self.desc]
        else:
            lines = [INDENT*depth + "["]
            lines += self.tgt._deepstr(depth+1)
            lines += self.src._deepstr(depth+1)
            lines += [INDENT*depth + "]"]
        return lines

    def deepstr(self, depth=0):
        lines = self._deepstr(depth)
        return '\n'.join(lines)

    @property
    def shape(self):
        return (self.dim,)

    @property
    def is_endo(self):
        assert self.dim > 0, "wah?"
        return self.src == self.tgt

    @property
    def is_identity(self):
        assert self.dim > 0, "wah?"
        return self.src == self.tgt and self == self.src.identity

#    def __mul__(lhs, rhs):
#        assert lhs.dim == rhs.dim 
#        n = lhs.dim
#        offset = lhs.cat.dim+lhs.cat.codim
#        assert n>=offset
#        try:
#            pair = lhs.topcat.Pair(n-offset, lhs, rhs)
#        except NotComposable:
#            print("__mul__: NotComposable")
#            print("lhs =")
#            print(lhs.deepstr())
#            print("rhs =")
#            print(rhs.deepstr())
#            raise
#        return pair
#
#    def __lshift__(lhs, rhs):
#        assert lhs.dim == rhs.dim 
#        n = lhs.dim
#        offset = lhs.cat.dim-1+lhs.cat.codim
#        assert n>=offset
#        return lhs.topcat.Pair(n-offset, lhs, rhs)
#
#    def __matmul__(lhs, rhs):
#        assert lhs.dim == rhs.dim 
#        n = lhs.dim
#        offset = lhs.cat.dim-2+lhs.cat.codim
#        assert n>=offset
#        return lhs.topcat.Pair(n-offset, lhs, rhs)

class Identity(Operator):
    def __init__(self, theory):
        Operator.__init__(self, theory, "identity", theory.cell1, (theory.cell0,), 
            postfix=True)

    def __call__(self, cell):
        expr = Operator.__call__(self, cell)
        build_cell(expr, cell, cell)
        return expr


class Pair(Operator):
    def __init__(self, theory, name, sort, codim):
        Operator.__init__(self, theory, name, sort, (sort, sort), inline=True)
        assert codim >= 0
        self.codim = codim

    def __call__(self, lhs, rhs):
        assert lhs.dim == rhs.dim, "use whisker first?"
        dim = lhs.dim
        #print("Pair:", lhs, rhs)
        codim = self.codim
        if codim == 0:
            if rhs.tgt != lhs.src:
                raise NotComposable("%s != %s"%(lhs.src, rhs.tgt))
            tgt = lhs.tgt
            src = rhs.src
        elif codim == 1:
            if lhs.src.src != rhs.src.tgt:
                raise NotComposable("%s != %s"%(lhs.src.src, rhs.src.tgt))
            # lhs.tgt.src == rhs.tgt.tgt
            tgt = pair(0, lhs.tgt, rhs.tgt)
            src = pair(0, lhs.src, rhs.src)
        elif codim == 2:
            if lhs.src.src.src != rhs.src.src.tgt:
                raise NotComposable("%s != %s"%(lhs.src.src.src, rhs.src.src.tgt))
            # -> lhs.tgt.src.src == rhs.tgt.src.tgt
            # -> lhs.tgt.tgt.src == rhs.tgt.tgt.tgt
            tgt = pair(1, lhs.tgt, rhs.tgt)
            src = pair(1, lhs.src, rhs.src)
        else:
            src, tgt = lhs, rhs
            for _ in range(codim):
                src = src.src
                tgt = tgt.src
            if rhs.tgt != lhs.src:
                raise NotComposable("%s != %s"%(lhs.src, rhs.tgt))
            tgt = pair(codim-1, lhs.tgt, rhs.tgt)
            src = pair(codim-1, lhs.src, rhs.src)
        expr = Operator.__call__(self, lhs, rhs)
        build_cell(expr, tgt, src)
        expr.lhs = lhs
        expr.rhs = rhs
        return expr


class Category(Theory):
    def __init__(self):
        Theory.__init__(self, Z3Solver())
        cell0 = Sort("cell0")
        cell1 = Sort("cell1")
        self.cell0 = cell0
        self.cell1 = cell1

        mul = Pair(self, "*", cell1, 0)
        self.add_operator(mul)
        #self.Operator("identity", cell1, [cell0], postfix=True)
        identity = Identity(self)
        self.add_operator(identity)

        # build theory
        Cell = self.Cell
        l = Cell(name="l")
        m = Cell(name="m")
        n = Cell(name="n")
        o = Cell(name="o")
    
        f = Cell(m, l)
        g = Cell(n, m)
        h = Cell(o, n)

        self.Rewrite( l.identity*l.identity, l.identity )
        self.Rewrite( f*l.identity, f )
        self.Rewrite( m.identity*f, f )
        self.Equation( (h*g)*f, h*(g*f) )

    def Cell(self, tgt=None, src=None, name=""):
        assert (tgt is None) == (src is None), (tgt, src)
        sort = self.cell0 if tgt is None else self.cell1
        cell = Cell(self, name, sort, tgt, src)
        self.Const(name, sort, cell) # brrr
        return cell

    def Nullary(self, tgt=None, src=None, name="c"):
        sort = self.cell0 if tgt is None else self.cell1
        expr = Theory.Nullary(self, name, sort)
        build_cell(expr, tgt, src, name)
        return expr

    def Iso(self, tgt, src, name="c"):
        cell = self.Const(tgt, src, name)
        inv = self.Const(src, tgt, name+"_i")
        cell.inv = inv
        inv.inv = cell
        self.Rewrite( inv*cell, src.identity )
        self.Rewrite( cell*inv, tgt.identity )
        return cell


def test_category():

    print("test_category")
    theory = Category()
    Cell = theory.Cell

    l = Cell(name="l")
    m = Cell(name="m")
    n = Cell(name="n")
    o = Cell(name="o")

    e = Cell(l, l)
    f = Cell(m, l)
    g = Cell(n, m)
    h = Cell(o, n)

    assert (h*g)*f == h*(g*f)
    assert h*n.identity == h
    assert n.identity*n.identity == n.identity
    assert h*n.identity*g*m.identity == h*g

    assert e*e != e
    assert f != f*e
    assert e*l.identity == e


def test_category_more():

    print("test_category_more")

    # -----------------------------------------
    # Test operations in a category

    theory = Category()
    Cell = theory.Cell

    # 0-cells
    l, m, n, o, p = [Cell(name=ch) for ch in 'lmnop']

    Il = l.identity
    assert Il*Il == Il

    # 1-cells
    A = Cell(m, l, "A")
    B = Cell(n, m, "B")
    C = Cell(o, n, "C")
    D = Cell(p, o, "D")
    In = n.identity
    Im = m.identity

    BA = B*A
    assert BA == B*A
    #assert BA != A
    assert BA.src == l
    assert BA.tgt == n

    assert A*Il == A
    assert Im*A == A
    assert Im*Im == Im

    assert (C*B)*A == C*(B*A)

    B*A
    C*B
    D*C
    C*(B*A)
    (C*B)*A
    D*(C*B)
    D*C*B
    
    cell = D*C*B*A
    
    assert cell == ((D*C)*B)*A
    assert cell == (D*C)*(B*A)
    assert cell == D*(C*(B*A))
    assert cell == D*((C*B)*A)
    assert cell == (D*(C*B))*A

#    theory.DEBUG = True
#    theory.solver.DEBUG = True
#    E = theory.Iso(m, l)
#    assert E*E.inv == Im
#    assert E.inv*E == Il
#
#    lhs = (B*E)*E.inv 
#
##    assert (B*E)*E.inv == B*(E*E.inv)
##    assert (B*E)*E.inv == B*m.identity
#    assert B*E*E.inv == B
#
#    return

    # 1-cell
    A = Cell(m, m)
    assert A != A*A


if __name__ == "__main__":
    print("\n\n")
    test_category()
    test_category_more()
    print("OK\n")

