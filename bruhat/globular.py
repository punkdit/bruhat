#!/usr/bin/env python

# https://ncatlab.org/nlab/show/globular+set

bars = ["-", "=", "≡", "≣"] # ...?

class Cell(object):
    def __init__(self, cat, *args):
        """
        call via either:
            Cell(tgt, src)
            Cell(name)
        """
        assert isinstance(cat, NCategory), cat
        tgt = None
        src = None
        if len(args)==1:
            name = args[0]
            n = 0
        else:
            assert len(args) in [2,3]
            tgt, src = args[:2]
            assert isinstance(src, Cell), src
            assert isinstance(tgt, Cell), tgt
            n = src.n+1
            assert src.n == tgt.n
            name = args[2] if args[2:] else ""
            #name = "(%s <%s%s%s %s)"%(tgt.name, bars[n-1], name, bars[n-1], src.name)
            name = "(%s <%s%s%s %s)"%(tgt, bars[n-1], name, bars[n-1], src)
        if n>1:
            # check globular conditions
            assert src.src == tgt.src
            assert src.tgt == tgt.tgt
        self.cat = cat
        self.n = n
        self.src = src
        self.tgt = tgt
        self.name = name

    def __str__(self):
        return self.name

    def __repr__(self):
        if self.tgt is not None:
            return "Cell(%s, %s)"%(self.tgt, self.src)
        else:
            return "Cell(%r)"%(self.name,)

    def __mul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = lhs.cat.level
        assert n>=offset
        return Pair(n-offset, lhs, rhs)

    def __lshift__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = lhs.cat.level-1
        assert n>=offset
        return Pair(n-offset, lhs, rhs)

    def __matmul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = lhs.cat.level-2
        assert n>=offset
        return Pair(n-offset, lhs, rhs)



class Pair(Cell):

    @classmethod
    def whisker(cls, lhs, rhs):
        while lhs.n < rhs.n:
            lhs = Cell(lhs, lhs)
        while rhs.n < lhs.n:
            rhs = Cell(rhs, rhs)
        assert lhs.n == rhs.n
        return (lhs, rhs)

    def __new__(cls, codim, lhs, rhs):
        assert isinstance(lhs, Cell)
        assert isinstance(rhs, Cell)
        cat = lhs.cat
        assert cat == rhs.cat
        paircache = cat.paircache
        key = (codim, lhs, rhs)
        if key in paircache:
            return paircache[key]
        ob = object.__new__(cls)
        paircache[key] = ob
        assert codim >= 0
        assert lhs.n == rhs.n, "use whisker first"
        if codim == 0:
            assert rhs.tgt == lhs.src, "not composable"
            tgt = lhs.tgt
            src = rhs.src
        elif codim == 1:
            assert lhs.src.src == rhs.src.tgt
            # lhs.tgt.src == rhs.tgt.tgt
            tgt = Pair(0, lhs.tgt, rhs.tgt)
            src = Pair(0, lhs.src, rhs.src)
        elif codim == 2:
            assert lhs.src.src.src == rhs.src.src.tgt
            # -> lhs.tgt.src.src == rhs.tgt.src.tgt
            # -> lhs.tgt.tgt.src == rhs.tgt.tgt.tgt
            tgt = Pair(1, lhs.tgt, rhs.tgt)
            src = Pair(1, lhs.src, rhs.src)
        else:
            tgt = Pair(codim-1, lhs.tgt, rhs.tgt)
            src = Pair(codim-1, lhs.src, rhs.src)
        Cell.__init__(ob, cat, tgt, src)
        ob.codim = codim
        ob.lhs = lhs
        ob.rhs = rhs
        ob.key = key
        return ob

    def __init__(self, codim, lhs, rhs):
        pass

    def __str__(self):
        return "(%s *%d %s)" % (self.lhs, self.codim, self.rhs)


class NCategory(object):
    """
    level = 0 # set
    level = 1 # category
    level = 2 # bicategory
    level = 3 # tricategory
    """

    def __init__(self, level=1):
        self.level = level
        self.paircache = {} # Pair object cache
        self.unitcache = {}

    def cell(self, *args):
        c = Cell(self, *args)
        return c

    def unit(self, cell):
        assert cell.cat is self
        cache = self.unitcache
        if cell in cache:
            return cache[cell]
        unit = Cell(self, cell, cell)
        cache[cell] = unit
        return unit

    def equate(self, tgt, src):
        if tgt == src:
            return
        assert isinstance(src, Pair), src
        paircache = self.paircache
        assert paircache[src.key] is src # do we care?
        paircache[src.key] = tgt
        return tgt


class Category(NCategory):
    def __init__(self):
        NCategory.__init__(self, 1)

    def cell(self, *args):
        cell = NCategory.cell(self, *args)
        assert 0<=cell.n<=1
        if cell.n > 0:
            src = self.unit(cell.src)
            self.equate(cell, cell*src)
            tgt = self.unit(cell.tgt)
            self.equate(cell, tgt*cell)
        return cell

    def unit(self, cell):
        I = NCategory.unit(self, cell)
        self.equate(I, I*I)
        return I


class Bicategory(NCategory):
    def __init__(self):
        NCategory.__init__(self, 2)

    #def reunit(self, A):

    def reassoc(self, *cells):
        n = len(cells)
        for cell in cells:
            assert cell.n == 0, "expected 0-cell: %s"%(cell,)
        if n<2:
            return None # ?
        if n==2:
            x, y = cells
            


def main():

    cat = NCategory(3)
    cell = cat.cell

    a, b, c = [cell(ch) for ch in 'abc']
    assert str(a) == 'a'

    f = cell(a, b)
    g = cell(a, b)
    u = cell(f, g)
    assert str(u) == "((a <-- b) <== (a <-- b))"

    assert u==u
    assert u != cell(f, g)

    u = cell(u, u)
    u = cell(u, u)

    g = cell(c, b)
    f = cell(b, a)
    uu = Pair(0, g, f)

    assert uu is Pair(0, g, f)
    assert uu == Pair(0, g, f)

    # -----------------------------------------
    # Test operations in a 0-category aka set

    cat = NCategory(0)
    cell = cat.cell

    # 0-cells
    l, m, n = [cell(ch) for ch in 'lmn']

    # there are no operations apart from ==
    assert l==l
    assert l!=m

    # -----------------------------------------
    # Test operations in a category

    cat = Category()
    cell = cat.cell

    # 0-cells
    l, m, n = [cell(ch) for ch in 'lmn']

    # 1-cells
    A = cell(m, l)
    B = cell(n, m)
    In = cat.unit(n)
    Im = cat.unit(m)
    Il = cat.unit(l)

    assert In == cat.unit(n)

    BA = B*A
    assert BA.src == l
    assert BA.tgt == n

    assert A*Il == A
    assert Im*A == A
    assert Im*Im == Im

    # -----------------------------------------
    # Test operations in a bicategory

    cat = NCategory(2)
    cell = cat.cell

    # 0-cells
    l, m, n, o = [cell(ch) for ch in 'lmno']

    # 1-cells
    A = cell(m, l)
    A1 = cell(m, l)
    A2 = cell(m, l)
    B = cell(n, m)
    B1 = cell(n, m)
    C = cell(o, n)

    reassoc = cell(
        (C<<B)<<A,
        C<<(B<<A))
    print(reassoc)

    # 2-cells
    f = cell(A1, A)
    g = cell(B1, B)
    f1 = cell(A2, A1)

    BA = B<<A
    assert BA.tgt == n
    assert BA.src == l

    gf = g<<f
    assert gf.tgt == (B1<<A1)
    assert gf.src == (B<<A)

    ff = f1*f
    assert ff.tgt == A2
    assert ff.src == A

    assert ff.codim == 0
    assert gf.codim == 1
    assert BA.codim == 0

    # -----------------------------------------
    # Test operations in a one object tricategory == monoidal bicategory

    cat = NCategory(3)
    cell = cat.cell

    # 0-cell
    star = cell("*")

    # 1-cells
    l, m, n = [cell(star, star, ch) for ch in 'lmn']

    # 2-cells
    A = cell(m, l)
    A1 = cell(m, l)
    A2 = cell(m, l)
    B = cell(n, m)
    B1 = cell(n, m)

    assert B.n == 2

    # 3-cells
    f = cell(A1, A)
    g = cell(B1, B)
    f1 = cell(A2, A1)

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

    main()

    print("OK\n")



