#!/usr/bin/env python

# https://ncatlab.org/nlab/show/globular+set

bars = ["-", "=", "≡", "≣"] # ...?

class Cell(object):
    def __init__(self, *args):
        """
        call via either:
            Cell(tgt, src)
            Cell(name)
        """
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
            name = "(%s <%s%s%s %s)"%(tgt.name, bars[n-1], name, bars[n-1], src.name)
        if n>1:
            # check globular conditions
            assert src.src == tgt.src
            assert src.tgt == tgt.tgt
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

    #level = 1 # categories
    #level = 2 # bicategories
    level = 3 # tricategories

    def __mul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = Cell.level
        assert n>=offset
        return Pair(n-offset, lhs, rhs)

    def __lshift__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = Cell.level-1
        assert n>=offset
        return Pair(n-offset, lhs, rhs)

    def __matmul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = Cell.level-2
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

    cache = {}
    def __new__(cls, codim, lhs, rhs):
        assert isinstance(lhs, Cell)
        assert isinstance(rhs, Cell)
        key = (codim, lhs, rhs)
        if key in cls.cache:
            return cls.cache[key]
        ob = object.__new__(cls)
        cls.cache[key] = ob
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
        Cell.__init__(ob, tgt, src)
        ob.codim = codim
        ob.lhs = lhs
        ob.rhs = rhs
        return ob

    def __init__(self, codim, lhs, rhs):
        pass

    def __str__(self):
        return "(%s *%d %s)" % (self.lhs, self.codim, self.rhs)



def main():

    a, b, c = [Cell(ch) for ch in 'abc']
    assert str(a) == 'a'

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
    uu = Pair(0, g, f)

    assert uu is Pair(0, g, f)
    assert uu == Pair(0, g, f)

    # -----------------------------------------
    # Test operations in a category

    Cell.level = 1

    # 0-cells
    l, m, n = [Cell(ch) for ch in 'lmn']

    # 1-cells
    A = Cell(m, l)
    B = Cell(n, m)

    BA = B*A
    assert BA.src == l
    assert BA.tgt == n

    # -----------------------------------------
    # Test operations in a bicategory

    Cell.level = 2

    # 0-cells
    l, m, n = [Cell(ch) for ch in 'lmn']

    # 1-cells
    A = Cell(m, l)
    A1 = Cell(m, l)
    A2 = Cell(m, l)
    B = Cell(n, m)
    B1 = Cell(n, m)

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

    ff = f1*f
    assert ff.tgt == A2
    assert ff.src == A

    assert ff.codim == 0
    assert gf.codim == 1
    assert BA.codim == 0

    # -----------------------------------------
    # Test operations in a one object tricategory == monoidal bicategory

    Cell.level = 3

    # 0-cell
    star = Cell("*")

    # 1-cells
    l, m, n = [Cell(star, star, ch) for ch in 'lmn']

    # 2-cells
    A = Cell(m, l)
    A1 = Cell(m, l)
    A2 = Cell(m, l)
    B = Cell(n, m)
    B1 = Cell(n, m)

    assert B.n == 2

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

    main()

    print("OK\n")



