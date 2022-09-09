#!/usr/bin/env python

import sys
sys.setrecursionlimit(40)
from time import time
start_time = time()

from bruhat.theory import Sort, Expr, Const, Operator, Theory, Debug, distinct

# https://ncatlab.org/nlab/show/globular+set
BARS = ["-", "=", "≡", "≣"] # ...?
OPS = ["*", "<<", "@", "??"]

class NotComposable(Exception):
    pass

INDENT = "  "


class Cell(Debug):
    def __init__(self, cat, tgt=None, src=None, name="?", expr=None):
        assert (tgt is None) == (src is None), (tgt, src)
        if tgt is None:
            dim = 0
            desc = name
        else:
            assert isinstance(tgt, Cell), tgt
            assert isinstance(src, Cell), src
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
        self.cat = cat
        assert expr is not None
        self.expr = expr

    def __str__(self):
        return self.desc

    def __repr__(self):
        if self.tgt is not None:
            return "Cell(%s, %s)"%(self.tgt, self.src)
        else:
            return "Cell(%r)"%(self.name,)

    def __eq__(lhs, rhs):
        return lhs.expr == rhs.expr

    def __hash__(self):
        assert 0

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

    def __mul__(lhs, rhs):
        lhs.info("__mul__", lhs, rhs)
        expr = lhs.expr * rhs.expr
        n = lhs.dim
        offset = lhs.cat.dim+lhs.cat.codim
        assert n>=offset
        try:
            cell = Pair(n-offset, lhs, rhs, expr)
        except NotComposable:
            print("__mul__: NotComposable")
            print("lhs =")
            print(lhs.deepstr())
            print("rhs =")
            print(rhs.deepstr())
            raise
        return cell

    def __lshift__(lhs, rhs):
        expr = lhs.expr << rhs.expr
        n = lhs.dim
        offset = lhs.cat.dim-1+lhs.cat.codim
        assert n>=offset
        cell = Pair(n-offset, lhs, rhs, expr)
        return cell

    def __matmul__(lhs, rhs):
        expr = lhs.expr @ rhs.expr
        n = lhs.dim
        offset = lhs.cat.dim-2+lhs.cat.codim
        assert n>=offset
        cell = Pair(n-offset, lhs, rhs, expr)
        return cell

    @property
    def identity(self):
        expr = self.expr.identity
        return Cell(self.cat, self, self, str(expr), expr)

    def __getattr__(self, attr):
        method = getattr(self.cat, attr)
        value = method(self)
        setattr(self, attr, value)
        return value

    def rewrite(self, other):
        self.expr.rewrite(other.expr)


class Pair(Cell):
    """
        A _composable pair of cell's.
    """

    def __init__(self, codim, lhs, rhs, expr):
        self.info("Pair", codim, lhs, rhs, expr)
        cat = lhs.cat
        assert codim >= 0
        assert lhs.dim == rhs.dim, "use whisker first"
        if codim == 0:
            if rhs.tgt != lhs.src:
                raise NotComposable("%s != %s"%(lhs.src, rhs.tgt))
            tgt = lhs.tgt
            src = rhs.src
        elif codim == 1:
            if lhs.src.src != rhs.src.tgt:
                raise NotComposable("%s != %s"%(lhs.src.src, rhs.src.tgt))
            # lhs.tgt.src == rhs.tgt.tgt
            #tgt = Pair(0, lhs.tgt, rhs.tgt)
            #src = Pair(0, lhs.src, rhs.src)
            tgt = lhs.tgt << rhs.tgt
            src = lhs.src << rhs.src
        elif codim == 2:
            assert 0, "todo"
            if lhs.src.src.src != rhs.src.src.tgt:
                raise NotComposable("%s != %s"%(lhs.src.src.src, rhs.src.src.tgt))
            # -> lhs.tgt.src.src == rhs.tgt.src.tgt
            # -> lhs.tgt.tgt.src == rhs.tgt.tgt.tgt
            tgt = Pair(1, lhs.tgt, rhs.tgt)
            src = Pair(1, lhs.src, rhs.src)
        else:
            assert 0, "todo"
            src, tgt = lhs, rhs
            for _ in range(codim):
                src = src.src
                tgt = tgt.src
            if rhs.tgt != lhs.src:
                raise NotComposable("%s != %s"%(lhs.src, rhs.tgt))
            tgt = Pair(codim-1, lhs.tgt, rhs.tgt)
            src = Pair(codim-1, lhs.src, rhs.src)
        self.info("\t", tgt, src)
        self.codim = codim
        self.lhs = lhs
        self.rhs = rhs
        #print("Pair:", self.expr)
#        self.key = (codim, lhs.key, rhs.key)
        Cell.__init__(self, cat, tgt, src, expr=expr)

    @property
    def shape(self):
        return (self.codim, self.dim)

    def __getitem__(self, i):
        return [self.lhs, self.rhs][i]

    def __str__(self):
        return "(%s *%d %s)" % (self.lhs, self.codim, self.rhs)

    def _deepstr(self, depth=0):
        lines  = [INDENT*depth + "(" + str(self.shape)]
        lines += self.lhs._deepstr(depth+1)
        lines += [INDENT*depth + "*%d"%self.codim]
        lines += self.rhs._deepstr(depth+1)
        lines += [INDENT*depth + ")"]
        return lines


class Globular(Debug):
    """
    dim = 0 # set
    dim = 1 # category
    dim = 2 # bicategory
    dim = 3 # tricategory
    """

    def __init__(self, dim, sorts=None):
        assert dim>=0
        self.dim = dim
        if sorts is None:
            sorts = [Sort("cell%d"%i) for i in range(dim+1)]
        self.sorts = sorts
        self.codim = len(sorts)-dim-1
        self.theory = Theory()

    def Cell(self, tgt=None, src=None, name="?"):
        assert (tgt is None) == (src is None), (tgt, src)
        sorts = self.sorts
        if tgt is None:
            sort = sorts[0] 
        else:
            assert tgt.dim < self.dim
            sort = sorts[tgt.dim+1]
        expr = self.theory.Const(name, sort)
        cell = Cell(self, tgt, src, name, expr)
        return cell


class Category(Globular):
    def __init__(self, dim=1):
        Globular.__init__(self, dim)
        cell0, cell1 = self.sorts[-2:]
        theory = self.theory
        Variable = theory.Variable
        Operator = theory.Operator
        Rewrite = theory.Rewrite
        Equation = theory.Equation

        Operator("identity", cell1, [cell0], postfix=True)
        Operator("*", cell1, [cell1, cell1], inline=True)

        # build theory
        l = Variable("l", cell0)
        m = Variable("m", cell0)
        n = Variable("n", cell0)
        o = Variable("o", cell0)
        f = Variable("f", cell1) # m <-- l
        g = Variable("g", cell1) # n <-- m
        h = Variable("h", cell1) # o <-- n

        Rewrite( l.identity*l.identity, l.identity )
        Rewrite( f*l.identity, f )
        Rewrite( m.identity*f, f )
        Equation( (h*g)*f, h*(g*f) )

    def Iso(self, tgt, src, name):
        cell = self.Cell(tgt, src, name)
        inv = self.Cell(src, tgt, "~"+name)
        cell.inv = inv
        inv.inv = cell
        ( inv*cell ).rewrite( src.identity )
        ( cell*inv ).rewrite( tgt.identity )
        return cell


class Bicategory(Category):
    def __init__(self, dim=2):
        Category.__init__(self, dim)

        cell0, cell1, cell2 = self.sorts[-3:]
        theory = self.theory

        Variable = theory.Variable
        Operator = theory.Operator
        Rewrite = theory.Rewrite
        Equation = theory.Equation

        Operator("identity", cell1, [cell0],               postfix=True)
        Operator("<<",       cell1, [cell1, cell1],        inline =True)
        Operator("<<",       cell2, [cell2, cell2],        inline =True)

        # These need to be Iso's, not sure if we can handle that using Rewrite's..?
        #Operator("lunitor",  cell2, [cell1],               postfix=True)
        #Operator("runitor",  cell2, [cell1],               postfix=True)
        #Operator("reassoc",  cell2, [cell1, cell1, cell1], )

        # build theory
        l = Variable("l", cell0)
        m = Variable("m", cell0)
        n = Variable("n", cell0)
        o = Variable("o", cell0)

        A = A0 = Variable("A0", cell1) # m <--A-- l
        A1     = Variable("A1", cell1) # m <--A-- l
        A2     = Variable("A2", cell1) # m <--A-- l
        B = B0 = Variable("B0", cell1) # m <--B-- l
        B1     = Variable("B1", cell1) # m <--B-- l
        B2     = Variable("B2", cell1) # m <--B-- l
        C      = Variable("C", cell1) # o <--C-- n

        f0 = Variable("f0", cell2) #  A1 <----- A0
        f1 = Variable("f1", cell2) #  A2 <----- A1
        g0 = Variable("g0", cell2) #  B1 <----- B0
        g1 = Variable("g1", cell2) #  B2 <----- B1

        Equation( (B<<A).identity, B.identity << A.identity )
        #Rewrite( (B<<A).identity, B.identity << A.identity )

        # Equation here causes infinite recursion...
        Rewrite( (g1*g0) << (f1*f0) , (g1 << f1) * (g0 << f0) )
        Rewrite( (g1 << f1) * (g0 << f0), (g1*g0)<<(f1*f0) )

    def lunitor(self, cell):
        assert cell.dim == 1
        lunitor = self.Iso(cell, cell.tgt.identity << cell, "lunitor")
        return lunitor

    def runitor(self, cell):
        assert cell.dim == 1
        if cell == cell.src.identity:
            runitor = cell.lunitor
        else:
            runitor = self.Iso(cell, cell << cell.src.identity, "runitor")
        return runitor

    def reassoc(self, cell):
        assert cell.dim == 1
        assert cell.shape == (0, 1), cell
        CB, A = cell
        assert CB.shape == (0, 1), cell
        C, B = CB
        reassoc = self.Iso(C<<(B<<A), cell, "reassoc")
        return reassoc


def test_category():

    print("test_category")
    cat = Category()
    Cell = cat.Cell

    l = Cell(name="l")
    m = Cell(name="m")
    n = Cell(name="n")
    o = Cell(name="o")

    e = Cell(l, l, "e")
    f = Cell(m, l, "f")
    g = Cell(n, m, "g")
    h = Cell(o, n, "h")


    assert (h*g)*f == h*(g*f)
    assert h*n.identity == h
    cells = [l, m, n, o]
    for cell in cells:
        assert cell.identity*cell.identity == cell.identity
    assert distinct([cell.identity for cell in cells])

    for morphism in [e, f, g, h]:
        assert morphism*morphism.src.identity == morphism
        assert morphism.tgt.identity*morphism == morphism
        assert morphism.tgt.identity*morphism*morphism.src.identity == morphism

    assert h*n.identity*g*m.identity == h*g

    assert e*e != e
    assert f != f*e
    assert e*l.identity == e


def test_category_more():

    print("test_category_more")

    # -----------------------------------------
    # Test operations in a category

    cat = Category()
    Cell = cat.Cell

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

    E = cat.Iso(m, l, "E")
    assert E*E.inv == Im
    assert E.inv*E == Il

    assert (B*E)*E.inv == B*(E*E.inv)
    assert (B*E)*E.inv == B*m.identity
    assert B*E*E.inv == B

    # 1-cell
    A = Cell(m, m)
    assert A != A*A


def test_bicategory():

    print("test_bicategory")
    cat = Bicategory()
    theory = cat.theory

    Cell = cat.Cell

    # 0-cells
    l = Cell(name="l")
    m = Cell(name="m")
    n = Cell(name="n")
    o = Cell(name="o")
    p = Cell(name="p")
    q = Cell(name="q")

    # identity 1-cells
    Il = l.identity
    assert Il<<Il != Il

    I = l.identity
    unitor = I.lunitor
    assert unitor == I.runitor

    i = I.identity
    assert i*i == i
    f = Cell(I, I)
    assert f*i == f
    assert i*f == f

    assert i*unitor == unitor
    II = I<<I
    assert II.identity == i<<i
    assert unitor*(i<<i) == unitor # FAIL

    assert II != I

    assert II.identity


    # 1-cells
    A0 = A = Cell(m, l, "A")
    A1 = Cell(m, l, "A1")
    A2 = Cell(m, l, "A2")
    A3 = Cell(m, l, "A3")
    B0 = B = Cell(n, m, "B")
    B1 = Cell(n, m, "B1")
    B2 = Cell(n, m, "B2")
    B3 = Cell(n, m, "B3")
    C0 = C = Cell(o, n, "C")
    C1 = Cell(o, n, "C1")
    D = Cell(p, o, "D")
    E = Cell(q, p, "E")

    assert A0.dim == 1

    BA = B0<<A0

    A0.identity # better than hasattr..
    A0.lunitor
    A0.runitor

    u = A0.identity
    uu = u*u
    assert uu == u

    BAI = BA << A.src.identity
    BA.identity
    BAI.reassoc
    BA.lunitor

    assert BA.identity == B0.identity<<A0.identity

    i = BA.identity
    assert i * i == i

    lhs = (B0.identity << A0.identity) * (B0.identity << A0.identity)
    rhs = (B0.identity * B0.identity) << (A0.identity * A0.identity)
    assert lhs == rhs
    assert lhs == BA.identity

    assert BA.lunitor.tgt == BA
    assert BA.lunitor.src == BA.tgt.identity << BA
    assert BA.runitor.tgt == BA
    assert BA.runitor.src == BA << BA.src.identity

    CBA = C<<BA
    assert CBA.lunitor.tgt == CBA
    assert CBA.lunitor.src == CBA.tgt.identity << CBA
    assert CBA.runitor.tgt == CBA
    assert CBA.runitor.src == CBA << CBA.src.identity

    # 2-cells
    f = f0 = Cell(A1, A0, "f0")
    f1 = Cell(A2, A1, "f1")
    f2 = Cell(A3, A2, "f2")
    f22 = Cell(A3, A2, "f22")
    g = g0 = Cell(B1, B0, "g0")
    g1 = Cell(B2, B1, "g1")
    g2 = Cell(B3, B2, "g2")
    h = h0 = Cell(C1, C0, "h0")

    cell = Cell(B1<<A1, B0<<A0, "fg")
    assert cell * (B0.identity<<A0.identity) == cell
    assert cell * (B0<<A0).identity == cell
    assert cell * cell.src.identity == cell
    assert cell.src.identity == (B0<<A0).identity

    f21 = f2*f1
    f210 = f21*f0
    f10 = f1*f0
    rhs = f2*f10
    assert f210 == (f2*f1)*f0
    assert f210 == f2*(f1*f0)

    assert f0.dim == 2

    assert f2 != f22
    assert f2*f1 != f22*f1

    ff = f1*f0
    assert ff == f1*f0

    assert A1.identity*f0 == f0*A0.identity
    assert (f2*f1)*f0 == f2*(f1*f0)

    assert (g0<<f0).src == g0.src << f0.src

    lhs = g0 << f0
    i = g0.src.identity << f0.src.identity
    assert lhs == lhs * i

    i = (g0.src << f0.src).identity
    assert lhs == lhs * i

#    gf = g<<f
#    lhs = gf * gf.src.lunitor 
#    rhs = gf.tgt.lunitor * (gf.src.tgt.identity.identity << gf)
#    assert lhs == rhs # FAIL

    gf = g0<<f0
    assert gf == g0<<f0

    assert gf.src == B0<<A0
    assert gf * BA.identity == gf

    lhs = (g1*g0) << (f1*f0)
    rhs = (g1<<f1) * (g0<<f0)
    assert lhs == rhs

    g210 = g2*g1*g0
    f210 = f2*f1*f0
    lhs = g210 << f210

    rhs = (g2<<f2) * (g1<<f1) * (g0<<f0)
    assert lhs == rhs

    rrhs = (g2<<f2) * ((g1<<f1) * (g0<<f0))
    assert rhs == rrhs

    def test_identity(f):
        assert f * f.src.identity == f
        assert f.tgt.identity * f == f
    test_identity(g<<f)
    test_identity((g*g.src.lunitor)<<f)
    test_identity((g1*g0)<<f)

    def test_assoc(f2, f1, f0):
        assert (f2*f1)*f0 == f2*(f1*f0)
    test_assoc(f2, f1, f0)
    test_assoc(g2<<f2, g1<<f1, g0<<f0)

    def test_iso(f):
        i = cat.Iso(f.src, f.src, "iso")
        assert i * i.inv == f.src.identity
        assert f * f.src.identity == f
        assert f * i * i.inv == f
#    test_iso(f)
#    test_iso(g<<f)
#    test_iso((g*g.src.lunitor)<<f)

    # naturality of unitors
    def test_unitors(f):
        lhs = f * f.src.lunitor 
        rhs = f.tgt.lunitor * (f.src.tgt.identity.identity << f)
        assert lhs == rhs
        lhs = f * f.src.runitor 
        rhs = f.tgt.runitor * (f << f.src.src.identity.identity)
        assert lhs == rhs

    Theory.DEBUG = True
    theory.dump()
    test_unitors(f)
    test_unitors(g)

    test_unitors(g<<f) # FAIL

    lhs = g * g.src.lunitor 
    rhs = g.tgt.lunitor * (g.src.tgt.identity.identity << g)
    assert lhs == rhs
    lhs = lhs << f
    rhs = rhs << f
    assert lhs == rhs

    gf = g<<f
    lhs = gf * gf.src.lunitor 
    I = g.src.tgt.identity
    assert gf.src.lunitor.src == I << (g.src << f.src)
    reassoc = ((I<<g.src)<<f.src).reassoc
    #test_assoc((g*g.src.lunitor) << f, reassoc.inv, reassoc) # Works here !!!
    assert lhs == ((g*g.src.lunitor) << f) * reassoc.inv
    assert reassoc.inv * reassoc == reassoc.src.identity
    assert lhs*reassoc == ((g*g.src.lunitor) << f) * reassoc.inv * reassoc
    test_assoc((g*g.src.lunitor) << f, reassoc.inv, reassoc) # FAIL !!!
    assert lhs*reassoc == ((g*g.src.lunitor) << f) * (reassoc.inv * reassoc)
    assert lhs*reassoc == ((g*g.src.lunitor) << f) * reassoc.src.identity
    test_identity( ((g*g.src.lunitor) << f) )
    assert lhs*reassoc == ((g*g.src.lunitor) << f)
    assert lhs == ((g.tgt.lunitor*(I.identity << g)) << f) * reassoc.inv
    rhs = gf.tgt.lunitor * (gf.src.tgt.identity.identity << gf)
    #assert rhs * reassoc 

    #assert lhs == rhs # FAIL 

    # _associators
    lhs = (C << B) << A
    rhs = C << (B << A)
    assert lhs.reassoc.src == lhs
    assert lhs.reassoc.tgt == rhs
    assert not hasattr(rhs, "reassoc")

    cell = (D<<C) << (B<<A)
    assert cell.reassoc

    # naturality of reassoc
    hgf = (h0 << g0) << f0
    assert hgf.tgt.reassoc * hgf == (h0<<(g0<<f0)) * hgf.src.reassoc

    # triangle equation
    def test_triangle(B, A):
        lhs = B.runitor << A.identity
        rhs = B.identity << A.lunitor
        rhs = rhs * lhs.src.reassoc
        assert lhs == rhs
    test_triangle(B, A)
    test_triangle(C<<B, A)
    test_triangle(C, B<<A)

    # fooling around with identity's
    I = l.identity
    i = I.identity
    assert I.lunitor == I.runitor
    mul = I.lunitor
    comul = mul.inv
    assert comul == I.runitor.inv
    assert mul * comul == i
    ii = i<<i
    assert comul * mul == ii
    lhs = mul*(mul << i)
    rhs = mul*(i<<mul) 
    reassoc = ((I<<I)<<I).reassoc
    rhs = rhs * reassoc
    assert lhs == rhs
    #assert (mul << i) * reassoc.inv * (i << comul) == ii # FAIL
    #assert (i<<mul) * reassoc * (comul << i) == ii # FAIL

    # pentagon equation
    def test_pentagon(D, C, B, A):
        src = ((D<<C)<<B)<<A
        tgt = D<<(C<<(B<<A))
        lhs = src.reassoc 
        lhs = lhs.tgt.reassoc * lhs
        rhs = ((D<<C)<<B).reassoc << A.identity
        rhs = rhs.tgt.reassoc * rhs
        rhs = (D.identity << ((C<<B)<<A).reassoc) * rhs
        assert lhs == rhs
    test_pentagon(D, C, B, A)
    test_pentagon(E<<D, C, B, A)
    test_pentagon(E, D<<C, B, A)

    # Eckman - Hilton
    U = m.identity
    u = U.identity
    f = Cell(U, U)
    g = Cell(U, U)
    assert f!=g

    assert f*u == f == u*f
    assert g*u == g == u*g
    assert f<<g == (f*u) << (u*g)
    assert         (f*u) << (u*g) == (f<<u)*(u<<g)
    #assert                           (f<<u)*(u<<g) == ???
    #assert f*g == g*f
    # not yet...

    


if __name__ == "__main__":
    print("\n\n")
    test_category()
    test_category_more()
    test_bicategory()
    print("OK\n")


