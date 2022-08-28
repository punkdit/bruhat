#!/usr/bin/env python

# https://ncatlab.org/nlab/show/globular+set

bars = ["-", "=", "≡", "≣"] # ...?


class Cell(object):
    """
        These are the generating elements, 
        out of which further expressions (Pair's) are built.
    """

    def __init__(self, cat, tgt=None, src=None, name="", parent=None):
        assert isinstance(cat, NCategory), cat
        assert (tgt is None) == (src is None), (tgt, src)
        if tgt is None:
            n = 0
            desc = name
        else:
            assert isinstance(src, Cell), src
            assert isinstance(tgt, Cell), tgt
            n = src.n+1
            assert src.n == tgt.n
            #desc = "(%s <%s%s%s %s)"%(tgt.name, bars[n-1], name, bars[n-1], src.name)
            desc = "(%s <%s%s%s %s)"%(tgt, bars[n-1], name, bars[n-1], src)
            # check globular conditions
            assert src.src == tgt.src
            assert src.tgt == tgt.tgt
        self.cat = cat
        self.n = n
        self.src = src
        self.tgt = tgt
        self.desc = desc
        self.parent = parent
        self._hash = None

    def __str__(self):
        return self.desc

    def __repr__(self):
        if self.tgt is not None:
            return "Cell(%s, %s)"%(self.tgt, self.src)
        else:
            return "Cell(%r)"%(self.name,)

    def get_root(self):
        while self.parent is not None:
            self = self.parent
        return self

    # WARNING:
    # The following __hash__, __eq__ and equate is a delicate dance 
    # !!!!!!!!

    def __hash__(self):
        # Here we exhibit paranoia about hash value changing!
        # Do not hash me before assembling into equivelances!
        dest = self.get_root()
        assert self._hash is None or self._hash == dest._hash
        value = id(dest)
        assert dest._hash is None or dest._hash == value
        dest._hash = value
        self._hash = value
        return value

    def __eq__(lhs, rhs):
        assert lhs.cat is rhs.cat
        lhs = lhs.get_root()
        rhs = rhs.get_root()
        return lhs is rhs

    def equate(lhs, rhs):
        assert lhs.cat is rhs.cat
        lhs = lhs.get_root()
        rhs = rhs.get_root()
        if lhs is rhs:
            return
        assert lhs._hash is None or rhs._hash is None
        if rhs._hash is not None:
            lhs.parent = rhs
        else:
            rhs.parent = lhs

    def __mul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = lhs.cat.dim+lhs.cat.codim
        assert n>=offset
        return lhs.cat.Pair(n-offset, lhs, rhs)

    def __lshift__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = lhs.cat.dim-1+lhs.cat.codim
        assert n>=offset
        return lhs.cat.Pair(n-offset, lhs, rhs)

    def __matmul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        offset = lhs.cat.dim-2+lhs.cat.codim
        assert n>=offset
        return lhs.cat.Pair(n-offset, lhs, rhs)


class Object(Cell):
    def __init__(self, cat, name, parent=None):
        self.cat = cat
        self.n = 0
        self.src = None
        self.tgt = None
        self.name = name
        self.parent = parent
        self._hash = None



class Pair(Cell):
    """
        A _composable pair of cell's.
    """

    def __init__(self, codim, lhs, rhs):
        cat = lhs.cat
        pair = cat.Pair
        assert codim >= 0
        assert lhs.n == rhs.n, "use whisker first"
        if codim == 0:
            assert rhs.tgt == lhs.src, "not _composable"
            tgt = lhs.tgt
            src = rhs.src
        elif codim == 1:
            assert lhs.src.src == rhs.src.tgt
            # lhs.tgt.src == rhs.tgt.tgt
            tgt = pair(0, lhs.tgt, rhs.tgt)
            src = pair(0, lhs.src, rhs.src)
        elif codim == 2:
            assert lhs.src.src.src == rhs.src.src.tgt
            # -> lhs.tgt.src.src == rhs.tgt.src.tgt
            # -> lhs.tgt.tgt.src == rhs.tgt.tgt.tgt
            tgt = pair(1, lhs.tgt, rhs.tgt)
            src = pair(1, lhs.src, rhs.src)
        else:
            tgt = pair(codim-1, lhs.tgt, rhs.tgt)
            src = pair(codim-1, lhs.src, rhs.src)
        Cell.__init__(self, cat, tgt, src)
        self.codim = codim
        self.lhs = lhs
        self.rhs = rhs

    # used for testing only..
    @property
    def key(self):
        return (self.codim, self.lhs, self.rhs)

    def __getitem__(self, i):
        return [self.lhs, self.rhs][i]

    def __str__(self):
        return "(%s *%d %s)" % (self.lhs, self.codim, self.rhs)



class NCategory(object):
    """
    dim = 0 # set
    dim = 1 # category
    dim = 2 # bicategory
    dim = 3 # tricategory
    """

    def __init__(self, dim, supercat=None):
        # I am a (the) homcat in my supercat
        self.dim = dim
        self.codim = 0 if supercat is None else supercat.codim+1
        self.supercat = supercat
        self.cache = {}

    def Object(self, name):
        codim = self.codim
        assert codim == 0
        cell = Object(self, name)
        return cell

    def Cell(self, tgt=None, src=None, name=""):
        assert (tgt is None) == (src is None), (tgt, src)
        cell = Cell(self, tgt, src, name)
        return cell

    @classmethod
    def whisker(cls, lhs, rhs):
        while lhs.n < rhs.n:
            lhs = self.Cell(lhs, lhs)
        while rhs.n < lhs.n:
            rhs = self.Cell(rhs, rhs)
        assert lhs.n == rhs.n
        return (lhs, rhs)

    def Pair(self, codim, lhs, rhs):
        assert isinstance(lhs, Cell)
        assert isinstance(rhs, Cell)
        cache = self.cache
        key = (codim, lhs, rhs)
        if key in cache:
            pair = cache[key]
        else:
            pair = Pair(codim, lhs, rhs)
            cache[key] = pair
        return pair

    def equate(self, lhs, rhs):
        lhs.equate(rhs)


class Set(NCategory):
    def __init__(self, supercat=None):
        NCategory.__init__(self, 0, supercat)

    # put the equate method here? do we care?


class Category(NCategory):
    def __init__(self, supercat=None):
        NCategory.__init__(self, 1, supercat)
        #self.homs = {}
        self.homcat = Set(self)

#    def Hom(self, tgt, src):
#        key = (tgt, src)
#        hom = self.homs.get(key)
#        if hom is None:
#            hom = Set()
#            self.homs[key] = hom
#        return hom

    def Object(self, name):
        cell = NCategory.Object(self, name)
        return cell

    def Cell(self, tgt=None, src=None, name=""):
        assert (tgt is None) == (src is None), (tgt, src)
        codim = self.codim
        # TODO: I own the Object's, my homcat owns the higher Cell's
        cell = NCategory.Cell(self, tgt, src, name)
        if cell.n < codim:
            assert 0, cell.n
        elif cell.n == codim:
            cell.unit = Cell(self, cell, cell)
            self.equate(cell.unit*cell.unit, cell.unit)
        elif cell.n == codim+1:
            src = cell.src.unit
            self.equate(cell, cell*src)
            self.equate(src, src*src)
            tgt = cell.tgt.unit
            self.equate(cell, tgt*cell)
            self.equate(tgt, tgt*tgt)
        else:
            assert 0, cell.n
        return cell

    def Pair(self, codim, lhs, rhs):
        assert codim>=0
        if codim > 0:
            return self.supercat.Pair(codim, lhs, rhs)
        compose = lambda lhs, rhs : NCategory.Pair(self, codim, lhs, rhs)
        pair = compose(lhs, rhs)
        if isinstance(rhs, Pair):
            # reassoc to the left
            a = lhs
            b, c = rhs
            other = compose(a*b, c)
            self.equate(pair, other)
        elif isinstance(lhs, Pair):
            # reassoc to the right
            a, b = lhs
            c = rhs
            other = compose(a, b*c)
            self.equate(pair, other)
        return pair

    def iso(self, tgt, src, name=""):
        "make an iso tgt<--src"
        cell = self.Cell(tgt, src, name)
        inv = self.Cell(src, tgt, "~"+name)
        self.equate(src.unit, inv*cell)
        self.equate(tgt.unit, cell*inv)
        cell.inv = inv
        inv.inv = cell
        return cell


class Bicategory(NCategory):
    def __init__(self, supercat=None):
        NCategory.__init__(self, 2, supercat)

        # we just bundle all the hom categories up into one python object
        self.homcat = Category(self)

    def Cell(self, tgt=None, src=None, name=""):
        assert (tgt is None) == (src is None), (tgt, src)
        # I own the Object's, my homcat owns the higher Cell's
        codim = self.codim
        homcat = self.homcat
        if tgt is None:
            assert codim == 0
            cell = NCategory.Cell(self, tgt, src, name)
            # this cell is an Object: make its unit 1-cell
            cell.unit = homcat.Cell(cell, cell)
        else:
            cell = homcat.Cell(tgt, src, name)
        return cell



def main():

    cat = NCategory(3)
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

    cat = NCategory(0)
    Cell = cat.Cell

    # 0-cells
    l, m = [Cell(name=ch) for ch in 'lm']

    # there are no operations apart from ==
    assert l==l
    assert l!=Cell(name="l")
    assert l!=m

    # -----------------------------------------
    # Test operations in a category

    cat = Category()
    Cell = cat.Cell

    # 0-cells
    l, m, n, o, p = [Cell(name=ch) for ch in 'lmnop']

    Il = l.unit
    assert Il*Il == Il

    # 1-cells
    A = Cell(m, l)
    B = Cell(n, m)
    C = Cell(o, n)
    D = Cell(p, o)
    In = n.unit
    Im = m.unit

    BA = B*A
    assert BA == B*A
    assert BA != A
    assert BA.src == l
    assert BA.tgt == n

    assert A*Il == A
    assert Im*A == A
    assert Im*Im == Im

    assert (C*B)*A == C*(B*A)

    cell = D*C*B*A
    
    assert cell == ((D*C)*B)*A
    assert cell == (D*C)*(B*A)
    assert cell == D*(C*(B*A))
    assert cell == D*((C*B)*A)
    assert cell == (D*(C*B))*A

    E = cat.iso(m, l)
    assert E*E.inv == Im
    assert E.inv*E == Il

    assert B*E*E.inv == B

    # 1-cell
    A = Cell(m, m)
    assert A != A*A

    # -----------------------------------------
    # Test operations in a Bicategory

    cat = Bicategory()
    Cell = cat.Cell

    # 0-cells
    l, m, n, o = [Cell(name=ch) for ch in 'lmno']

    # unit 1-cells
    u = l.unit
    uu = u<<u
    assert uu != u
    w = u.unit
    assert w*w == w

    # 1-cells
    A0 = Cell(m, l)
    A1 = Cell(m, l)
    A2 = Cell(m, l)
    A3 = Cell(m, l)
    B0 = Cell(n, m)
    B1 = Cell(n, m)

    assert A0.cat is cat.homcat

    BA = B0<<A0
    assert BA.cat is cat.homcat

    u = A0.unit
    uu = u*u
    assert uu == u

    # 2-cells
    f0 = Cell(A1, A0)
    f1 = Cell(A2, A1)
    f2 = Cell(A3, A2)
    f22 = Cell(A3, A2)
    g0 = Cell(B1, B0)

    assert f2 != f22
    assert f2*f1 != f22*f1

    gf = g0<<f0
    assert gf is g0<<f0
    assert gf.key in cat.cache

    ff = f1*f0
    assert ff is f1*f0
    assert ff.key not in cat.cache
    assert ff.key in cat.homcat.cache # this Pair lives in the homcat

    assert A1.unit*f0 == f0*A0.unit
    assert (f2*f1)*f0 == f2*(f1*f0)

    #return

    # -----------------------------------------
    # Test operations in a bicategory

    cat = NCategory(2)
    Cell = cat.Cell

    # 0-cells
    l, m, n, o = [Cell(name=ch) for ch in 'lmno']

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

    cat = NCategory(3)
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



