#!/usr/bin/env python

"""
previous version: globular.py

"""

from time import time
start_time = time()

import z3


class Solver(object):
    def __init__(self):
        self.solver = z3.Solver()
        self.sort = z3.DeclareSort("THE_SORT")
        self.v_count = 0
        self.ops = set()

    def get_var(self, stem="v"):
        name = stem+str(self.v_count)
        self.v_count += 1
        v = z3.Const(name, self.sort)
        return v

    def get_op(self, name, arity=2):
        assert name not in self.ops
        f = z3.Function(name, [self.sort]*(arity+1))
        self.ops.add(name)
        return f

    def equate(self, lhs, rhs):
        self.solver.add( lhs == rhs )

    def is_equal(self, lhs, rhs):
        solver = self.solver
        solver.push()
        solver.add( lhs != rhs )
        result = solver.check()
        solver.pop()
        return result == z3.unsat


def test_solver():

    solver = Solver()
    a = solver.get_var("a")
    b = solver.get_var("b")
    c = solver.get_var("c")
    d = solver.get_var("d")

    solver.equate(a, b)
    assert solver.is_equal(a, b)
    solver.equate(b, c)
    assert solver.is_equal(a, c)
    assert not solver.is_equal(a, d)

    mul = solver.get_op("*", 2)

    lhs = mul(mul(a, b), c)
    rhs = mul(a, mul(b, c))
    solver.equate( lhs, rhs )

    assert solver.is_equal( mul(lhs, d), mul(rhs, d) )
    assert not solver.is_equal( mul(lhs, d), mul(d, rhs) )



# https://ncatlab.org/nlab/show/globular+set

BARS = ["-", "=", "≡", "≣"] # ...?
OPS = ["*", "<<", "@", "??"]

class NotComposable(Exception):
    pass

INDENT = "  "


class Cell(object):
    """
        These are the generating elements, 
        out of which further expressions (Pair's) are built.
    """

    inv = None # for isomorphisms

    def __init__(self, cat, tgt=None, src=None, name="", parent=None):
        assert isinstance(cat, Globular), cat
        assert (tgt is None) == (src is None), (tgt, src)
        if tgt is None:
            dim = 0
            desc = name
        else:
            assert isinstance(src, Cell), src
            assert isinstance(tgt, Cell), tgt
            dim = src.dim+1
            assert src.dim == tgt.dim
            #desc = "(%s <%s%s%s %s)"%(tgt.name, BARS[dim-1], name, BARS[dim-1], src.name)
            desc = "(%s <%s%s%s %s)"%(tgt, BARS[dim-1], name, BARS[dim-1], src)
            # check globular conditions
            assert src.src == tgt.src
            assert src.tgt == tgt.tgt
        self.cat = cat
        self.topcat = cat.topcat
        self.dim = dim
        self.src = src
        self.tgt = tgt
        self.desc = desc
        self.parent = parent
        self._hash = None
        self.is_decorated = 0
        self.expr = cat.solver.get_var("cell")
        self.solver = cat.solver
        self.key = id(self)

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

    def __hash__(self):
        assert 0

    def __eq__(lhs, rhs):
        return lhs.solver.is_equal( lhs.expr, rhs.expr )

    def equate(lhs, rhs):
        solver = lhs.solver
        if lhs != rhs:
            solver.equate( lhs.expr , rhs.expr )

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
        assert lhs.dim == rhs.dim 
        n = lhs.dim
        offset = lhs.cat.dim+lhs.cat.codim
        assert n>=offset
        try:
            pair = lhs.topcat.Pair(n-offset, lhs, rhs)
        except NotComposable:
            print("__mul__: NotComposable")
            print("lhs =")
            print(lhs.deepstr())
            print("rhs =")
            print(rhs.deepstr())
            raise
        return pair

    def __lshift__(lhs, rhs):
        assert lhs.dim == rhs.dim 
        n = lhs.dim
        offset = lhs.cat.dim-1+lhs.cat.codim
        assert n>=offset
        return lhs.topcat.Pair(n-offset, lhs, rhs)

    def __matmul__(lhs, rhs):
        assert lhs.dim == rhs.dim 
        n = lhs.dim
        offset = lhs.cat.dim-2+lhs.cat.codim
        assert n>=offset
        return lhs.topcat.Pair(n-offset, lhs, rhs)


class Pair(Cell):
    """
        A _composable pair of cell's.
    """

    def __init__(self, op, codim, lhs, rhs):
        cat = lhs.cat
        pair = cat.Pair
        assert codim >= 0
        assert lhs.dim == rhs.dim, "use whisker first"
        #print("Pair:", op, lhs, rhs)
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
        Cell.__init__(self, cat, tgt, src)
        self.expr = op(lhs.expr, rhs.expr)
        self.codim = codim
        self.lhs = lhs
        self.rhs = rhs
        #print("Pair:", self.expr)
        self.key = (codim, lhs.key, rhs.key)

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

    _identity = None
    @property
    def identity(self):
        if self._identity is None:
            cat = self.cat.supercat
            self._identity = cat.Identity(self)
        return self._identity

    _lunitor = None
    @property
    def lunitor(self):
        if self._lunitor is None:
            cat = self.cat.supercat
            self._lunitor = cat.LeftUnitor(self)
        return self._lunitor

    _runitor = None
    @property
    def runitor(self):
        if self._runitor is None:
            cat = self.cat.supercat
            self._runitor = cat.RightUnitor(self)
        return self._runitor

    _inv = None
    @property
    def inv(self):
        if self._inv is None:
            cat = self.cat.supercat
            self._inv = cat.Inv(self)
        return self._inv




class Globular(object):
    """
    dim = 0 # set
    dim = 1 # category
    dim = 2 # bicategory
    dim = 3 # tricategory
    """

    DEBUG = False

    def __init__(self, dim, supercat=None, topcat=None):
        # I am a (the) homcat in my supercat
        self.dim = dim
        self.codim = 0 if supercat is None else supercat.codim+1
        self.supercat = supercat
        if topcat is None:
            # i am the topcat, meow!
            self.topcat = self
            self.solver = Solver()
            ops = []
            for d in range(dim):
                ops.append(self.solver.get_op(OPS[d]))
            self.ops = ops
            self.cache = {}
        else:
            self.topcat = topcat
            self.solver = topcat.solver
            self.ops = topcat.ops
            self.cache = topcat.cache

    def info(self, *msg):
        if self.DEBUG:
            print(*msg)

    def Cell(self, tgt=None, src=None, name=""):
        assert (tgt is None) == (src is None), (tgt, src)
        cell = Cell(self, tgt, src, name)
        return cell

    @classmethod
    def whisker(cls, lhs, rhs):
        while lhs.dim < rhs.dim:
            lhs = self.Cell(lhs, lhs)
        while rhs.dim < lhs.dim:
            rhs = self.Cell(rhs, rhs)
        assert lhs.dim == rhs.dim
        return (lhs, rhs)

    def Pair(self, codim, lhs, rhs):
        assert isinstance(lhs, Cell)
        assert isinstance(rhs, Cell)
        pair = Pair(self.ops[codim], codim, lhs, rhs)
        if pair.key in self.cache:
            pair = self.cache[pair.key]
        else:
            self.cache[pair.key] = pair
        return pair

    def equate(self, lhs, rhs):
        lhs.equate(rhs)


class Set(Globular):
    def __init__(self, supercat=None, topcat=None):
        Globular.__init__(self, 0, supercat, topcat)

    # put the equate method here? do we care?


class Category(Globular):
    def __init__(self, supercat=None, topcat=None):
        Globular.__init__(self, 1, supercat, topcat)
        self.homcat = Set(self, self.topcat)
        #self.homs = {}

    def Object(self, name):
        cell = Globular.Object(self, name)
        return cell

    def decorate_pair(self, pair):
        compose = lambda lhs, rhs : Globular.Pair(self.topcat, pair.codim, lhs, rhs)
        lhs, rhs = pair
        self.info("Category.decorate_pair:", pair)
        if pair.codim + pair.dim != self.codim + self.dim:
            self.info("Category.decorate_pair: bailing")
            return
        if isinstance(rhs, Pair) and rhs.shape == pair.shape:
            # reassoc to the left
            self.info("Category.decorate_pair: left reassoc")
            a = lhs
            b, c = rhs
            other = compose(a*b, c)
            self.equate(pair, other)
        elif isinstance(lhs, Pair) and lhs.shape == pair.shape:
            # reassoc to the right
            self.info("Category.decorate_pair: right reassoc")
            a, b = lhs
            c = rhs
            bc = b*c
            other = compose(a, bc)
            self.equate(pair, other)
            assert pair == other
            if self.DEBUG:
                pair.debug_eq()
                other.debug_eq()

    def decorate(self, cell):
        assert isinstance(cell, Cell)
        self.info("Category.decorate: cell.is_decorated=%d" % cell.is_decorated)
        if cell.is_decorated >= self.dim:
            return
        cell.is_decorated = self.dim
        codim = self.codim
        if cell.dim < codim:
            assert 0, cell.dim
        elif isinstance(cell, Pair):
            self.decorate_pair(cell)
        elif cell.dim == codim:
            assert not hasattr(cell, "identity"), "wup"
            cell.identity = Cell(self, cell, cell)
            e = cell.identity * cell.identity
            self.equate(cell.identity*cell.identity, cell.identity)
        elif cell.dim == codim+1 and cell.is_identity:
            assert 0, "pass!"
        elif cell.dim == codim+1:
            if isinstance(cell, Pair) and (cell.lhs.is_identity or cell.rhs.is_identity):
                pass
            else:
                src = cell.src.identity
                self.equate(cell, cell*src)
                self.equate(src, src*src)
                tgt = cell.tgt.identity
                self.equate(cell, tgt*cell)
                self.equate(tgt, tgt*tgt)
        else:
            print("cell =", cell)
            assert 0, cell.dim
        #if isinstance(cell, Pair):
        #    self.decorate_pair(cell)

    def Cell(self, tgt=None, src=None, name=""):
        assert (tgt is None) == (src is None), (tgt, src)
        # todo: I own the Object's, my homcat owns the higher Cell's (?)
        cell = Globular.Cell(self, tgt, src, name)
        self.decorate(cell)
        return cell

    def Pair(self, codim, lhs, rhs):
        assert codim>=0
        assert lhs.dim == rhs.dim
        if self.codim + codim + 1 > lhs.dim:
            return self.supercat.Pair(codim, lhs, rhs)
        self.info("  Category.Pair", self.codim, codim, lhs.dim, lhs, rhs)
        compose = lambda lhs, rhs : Globular.Pair(self, codim, lhs, rhs)
        pair = compose(lhs, rhs)
        self.decorate(pair)
        self.info("  Category.Pair: shape=", pair.shape)
        return pair

    def Iso(self, tgt, src, name=""):
        "make an invertible tgt<--src"
        cell = self.Cell(tgt, src, name)
        inv = self.Cell(src, tgt, "~"+name)
        self.equate(src.identity, inv*cell)
        self.equate(tgt.identity, cell*inv)
        cell.inv = inv
        inv.inv = cell
        return cell


class Bicategory(Globular):
    def __init__(self, supercat=None, topcat=None):
        Globular.__init__(self, 2, supercat, topcat)

        # we just bundle all the hom categories up into one python object
        self.homcat = Category(self, self.topcat)

    def decorate_cell2(self, cell):
        assert cell.dim == 2
        self.info("Bicategory.decorate_cell2:", cell)

        lhs = cell * cell.src.lunitor 
        rhs = cell.tgt.lunitor * (cell.src.tgt.identity.identity << cell)
        self.equate(lhs, rhs)
        
        lhs = cell * cell.src.runitor 
        rhs = cell.tgt.runitor * (cell << cell.src.src.identity.identity)
        self.equate(lhs, rhs)

    def decorate_pair(self, pair):
        lhs, rhs = pair
        self.info("Bicategory.decorate_pair:")
        self.info("\tpair.shape =", pair.shape)
        self.info("\tlhs.shape =", lhs.shape)
        self.info("\trhs.shape =", rhs.shape)
        homcat = self.homcat
        found = False
        if pair.shape == (0, 2) and lhs.shape == (1, 2) and rhs.shape == (1, 2):
            # interchange law: 
            # pair = (a<<b) * (c<<d)
            a, b = lhs
            c, d = rhs
            other = (a*c) << (b*d)
            self.equate(pair, other)
            found = True
        if pair.shape == (1, 2) and lhs.shape == (0, 2) and rhs.shape == (0, 2):
            # interchange law: 
            # pair = (a*b) << (c*d)
            a, b = lhs
            c, d = rhs
            other = (a<<c) * (b<<d)
            self.equate(pair, other)
            found = True
        if pair.shape == (1, 2) and lhs.shape == (1, 2) \
                and (not pair.is_endo or pair != pair.src.identity): # arfff! 
            # naturality of reassoc
            h, g = lhs
            f = rhs
            C, B = lhs.src
            A = rhs.src
            src = pair.src
            assert src == ((C<<B)<<A)
            assert src.shape == (0,1) and src.lhs.shape == (0,1)
            assert "reassoc" in src.__dict__
            reassoc = src.reassoc
            self.equate(pair.tgt.reassoc * pair, (h<<(g<<f)) * src.reassoc)
        if pair.shape == (0, 1) and lhs.shape == (0, 1):
            # reassoc to the right
            self.info("decorate_pair: pair has reassoc", pair)
            C, B = lhs
            A = rhs
            other = C << (B<<A)
            pair.reassoc = homcat.Iso(other, pair)
            if B.is_endo and B == C.src.identity:
                # triangle equation
                lhs = C.runitor << A.identity
                rhs = (C.identity << A.lunitor) * pair.reassoc
                self.equate(lhs, rhs)
            found = True
        if pair.shape == (0, 1) and lhs.shape == (0, 1) and lhs.lhs.shape == (0, 1):
            # pentagon equation
            A = pair.rhs
            B = pair.lhs.rhs
            C = pair.lhs.lhs.rhs
            D = pair.lhs.lhs.lhs
            src = ((D<<C)<<B)<<A
            tgt = D<<(C<<(B<<A))
            lhs = src.reassoc 
            lhs = lhs.tgt.reassoc * lhs
            rhs = ((D<<C)<<B).reassoc << A.identity
            rhs = rhs.tgt.reassoc * rhs
            rhs = (D.identity << ((C<<B)<<A).reassoc) * rhs
            self.equate(lhs, rhs)
            found = True
        #if pair.dim == 2:
        #    self.decorate_cell2(pair)
        if not found:
            self.info("Bicategory.decorate_pair: decorate not found")

    def decorate(self, cell):
        assert isinstance(cell, Cell)
        assert self.codim == 0 # todo
        if cell.is_decorated >= self.dim:
            return # <------------ return
        cell.is_decorated = self.dim
        if isinstance(cell, Pair):
            self.decorate_pair(cell)
            return # <------------ return
        codim = self.codim
        homcat = self.homcat
        if cell.dim == 0:
            # this cell is an Object: make its identity 1-cell
            assert not hasattr(cell, "identity"), "wup"
            identity = homcat.Cell(cell, cell)
            cell.identity = identity
            lunitor = homcat.Iso(identity, identity<<identity)
            identity.lunitor = identity.runitor = lunitor
        elif cell.dim == 1:
            # this cell is a 1-cell, make its unit'ors
            if cell.__class__ is Cell:
                assert hasattr(cell, "identity"), "wup"
                cell.lunitor = homcat.Iso(cell, cell.tgt.identity << cell, "l")
                cell.runitor = homcat.Iso(cell, cell << cell.src.identity, "r")
            # argh, we have to do this lazily:
            #elif cell.shape == (0, 1):
            #    lhs, rhs = cell
            #    cell.lunitor = lhs.lunitor << rhs.identity
            #    cell.runitor = lhs.identity << rhs.runitor
            #else:
            #    assert 0, cell
        elif cell.dim == 2:
            self.decorate_cell2(cell)
        else:
            assert 0, cell

    def Identity(self, cell):
        assert isinstance(cell, Pair)
        assert cell.shape == (0, 1), cell
        lhs, rhs = cell
        identity = lhs.identity << rhs.identity
        return identity

    def LeftUnitor(self, cell):
        assert isinstance(cell, Pair)
        assert cell.shape == (0, 1), cell
        lhs, rhs = cell
        lunitor = lhs.lunitor << rhs.identity
        lunitor = lunitor * lunitor.src.reassoc.inv
        return lunitor

    def RightUnitor(self, cell):
        assert isinstance(cell, Pair)
        assert cell.shape == (0, 1), cell
        reassoc = (cell<<cell.src.identity).reassoc
        lhs, rhs = cell
        runitor = lhs.identity << rhs.runitor
        runitor = runitor * reassoc
        return runitor

    def Inv(self, cell):
        assert isinstance(cell, Pair)
        assert cell.shape == (1, 2), cell
        lhs, rhs = cell
        inv = None
        if lhs.inv and rhs.inv:
            inv = lhs.inv << rhs.inv
        return inv

    def Cell(self, tgt=None, src=None, name=""):
        assert (tgt is None) == (src is None), (tgt, src)
        # I own the Object's, my homcat owns the higher Cell's
        codim = self.codim
        homcat = self.homcat
        if tgt is None:
            assert codim == 0
            cell = Globular.Cell(self, tgt, src, name)
        else:
            assert codim == 0 # todo: lift this restriction
            cell = homcat.Cell(tgt, src, name)
        self.decorate(cell)
        return cell

    def Pair(self, codim, lhs, rhs):
        self.info("Bicategory.Pair:")
        self.info("\t", lhs)
        self.info("\t", rhs)
        #if codim > 0:
        #    return self.supercat.Pair(codim, lhs, rhs)
        compose = lambda lhs, rhs : Globular.Pair(self, codim, lhs, rhs)
        pair = compose(lhs, rhs)
        self.info("Bicategory.Pair: pair.is_decorated=%s"%pair.is_decorated)
        self.homcat.decorate(pair)
        self.decorate(pair)
        return pair


class Tricategory(Globular):
    def __init__(self, supercat=None, topcat=None):
        Globular.__init__(self, 3, supercat, topcat)

        # we just bundle all the hom bicategories up into one python object
        self.homcat = Bicategory(self, self.topcat)


def test_category():

    print("test_category")

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

    E = cat.Iso(m, l)
    assert E*E.inv == Im
    assert E.inv*E == Il

    assert B*E*E.inv == B

    # 1-cell
    A = Cell(m, m)
    assert A != A*A


def test_bicategory():

    print("test_bicategory")

    # -----------------------------------------
    # Test operations in a Bicategory

    cat = Bicategory()
    homcat = cat.homcat
    assert homcat.topcat is cat
    Cell = cat.Cell

    # 0-cells
    l = Cell(name="l")
    m = Cell(name="m")
    n = Cell(name="n")
    o = Cell(name="o")
    p = Cell(name="p")
    q = Cell(name="q")

    # identity 1-cells
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
    assert unitor*(i<<i) == unitor

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

    assert A0.cat is cat.homcat

    BA = B0<<A0
    assert BA.cat is cat.homcat

    assert A0.is_decorated == 2
    A0.identity # better than hasattr..
    A0.lunitor
    A0.runitor

    u = A0.identity
    uu = u*u
    assert uu == u

    BAI = BA << A.src.identity
    BAI.reassoc
    BA.identity
    BA.lunitor

    assert BA.identity == B0.identity<<A0.identity

    i = BA.identity
    assert i.topcat is cat
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

    f21 = f2*f1
    f210 = f21*f0
    f10 = f1*f0
    rhs = f2*f10
    assert f210 is (f2*f1)*f0
    assert f210 == f2*(f1*f0)

    assert f0.dim == 2

    assert f2 != f22
    assert f2*f1 != f22*f1

    ff = f1*f0
    assert ff is f1*f0
    assert ff.key in cat.cache # this Pair lives in the cat

    assert A1.identity*f0 == f0*A0.identity
    assert (f2*f1)*f0 == f2*(f1*f0)

    gf = g0<<f0
    assert gf is g0<<f0
    assert gf.key in cat.cache

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
        i = homcat.Iso(f.src, f.src)
        assert i * i.inv == f.src.identity
        assert f * f.src.identity == f
        assert f * i * i.inv == f
    test_iso(f)
    test_iso(g<<f)
    test_iso((g*g.src.lunitor)<<f)

    # naturality of unitors
    def test_unitors(f):
        lhs = f * f.src.lunitor 
        rhs = f.tgt.lunitor * (f.src.tgt.identity.identity << f)
        assert lhs == rhs
        lhs = f * f.src.runitor 
        rhs = f.tgt.runitor * (f << f.src.src.identity.identity)
        assert lhs == rhs

    test_unitors(f)
    test_unitors(g)
    #test_unitors(g<<f) # FAIL

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


def test_globular():
    Globular.DEBUG = False

    print("test_globular")

    cat = Globular(3)
    Cell = cat.Cell

    a, b, c = [Cell(name=ch) for ch in 'abc']
    assert str(a) == 'a'

    assert a==a
    assert a != Cell(name="a")
    assert a != b

    f = Cell(a, b)
    g = Cell(a, b)
    u = Cell(f, g)
    assert str(u) == "((a <-- b) <== (a <-- b))"

    assert u==u
    assert u != Cell(f, g)

    u = Cell(u, u)
    u = Cell(u, u)

    g = Cell(c, b)
    f = Cell(b, a)
    uu = cat.Pair(0, g, f) 

    # -----------------------------------------
    # Test operations in a 0-category aka set

    cat = Globular(0)
    Cell = cat.Cell

    # 0-cells
    l, m = [Cell(name=ch) for ch in 'lm']

    # there are no operations apart from ==
    assert l==l
    assert l!=Cell(name="l")
    assert l!=m


    # -----------------------------------------
    # Test operations in a bicategory

    cat = Globular(2)
    Cell = cat.Cell

    # 0-cells
    l, m, n, o = [Cell(name=ch) for ch in 'lmno']

    u = Cell(l, l)
    uu = u<<u

    # 1-cells
    A = Cell(m, l)
    A1 = Cell(m, l)
    A2 = Cell(m, l)
    B = Cell(n, m)
    B1 = Cell(n, m)
    C = Cell(o, n)

    reassoc = Cell(
        (C<<B)<<A,
        C<<(B<<A))

    # 2-cells
    f = Cell(A1, A)
    g = Cell(B1, B)
    f1 = Cell(A2, A1)

    BA = B<<A
    assert BA.tgt == n
    assert BA.src == l

    gf = g<<f
    assert gf.tgt == (B1<<A1)
    assert gf.src == (B<<A)
    assert gf.codim == 1

    ff = f1*f
    assert ff.tgt == A2
    assert ff.src == A

    assert ff.codim == 0
    assert gf.codim == 1
    assert BA.codim == 0

    # -----------------------------------------
    # Test operations in a one object tricategory == monoidal bicategory

    cat = Globular(3)
    Cell = cat.Cell

    # 0-cell
    star = Cell(name="*")

    # 1-cells
    l, m, n = [Cell(star, star, ch) for ch in 'lmn']

    # 2-cells
    A = Cell(m, l)
    A1 = Cell(m, l)
    A2 = Cell(m, l)
    B = Cell(n, m)
    B1 = Cell(n, m)

    assert B.dim == 2

    # 3-cells
    f = Cell(A1, A)
    g = Cell(B1, B)
    f1 = Cell(A2, A1)

    # operations on 1-cell's
    mn = m@n

    # operations on 2-cell's
    AA = A@A1
    assert AA.src == l@l
    assert AA.tgt == m@m

    BA = B<<A
    assert BA.tgt == n, BA.tgt
    assert BA.src == l, BA.src

    # operations on 3-cell's
    ff = f@f
    assert ff.src == A@A
    assert ff.tgt == A1@A1

    gf = g<<f
    assert gf.tgt == (B1<<A1)
    assert gf.src == (B<<A)

    ff = f1*f
    assert ff.tgt == A2
    assert ff.src == A





if __name__ == "__main__":

    print("\n\n")
    test_solver()
    test_category()
    test_bicategory()
    test_globular()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







