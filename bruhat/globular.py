#!/usr/bin/env python

# https://ncatlab.org/nlab/show/globular+set

bars = ["-", "=", "≡", "≣"] # ...?


class NotComposable(Exception):
    pass

INDENT = "  "


class Cell(object):
    """
        These are the generating elements, 
        out of which further expressions (Pair's) are built.
    """

    inv = None # for isomorphisms
    shape = None # for Pair's

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
            #desc = "(%s <%s%s%s %s)"%(tgt.name, bars[dim-1], name, bars[dim-1], src.name)
            desc = "(%s <%s%s%s %s)"%(tgt, bars[dim-1], name, bars[dim-1], src)
            # check globular conditions
            assert src.src == tgt.src
            assert src.tgt == tgt.tgt
        self.cat = cat
        self.dim = dim
        self.src = src
        self.tgt = tgt
        self.desc = desc
        self.parent = parent
        self._hash = None
        self.is_decorated = 0

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

    def deepstr(self):
        lines = self._deepstr()
        return '\n'.join(lines)

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
        assert lhs.dim == rhs.dim 
        n = lhs.dim
        offset = lhs.cat.dim+lhs.cat.codim
        assert n>=offset
        try:
            pair = lhs.cat.Pair(n-offset, lhs, rhs)
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
        return lhs.cat.Pair(n-offset, lhs, rhs)

    def __matmul__(lhs, rhs):
        assert lhs.dim == rhs.dim 
        n = lhs.dim
        offset = lhs.cat.dim-2+lhs.cat.codim
        assert n>=offset
        return lhs.cat.Pair(n-offset, lhs, rhs)


class Object(Cell):
    def __init__(self, cat, name, parent=None):
        self.cat = cat
        self.dim = 0
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
        self.codim = codim
        self.lhs = lhs
        self.rhs = rhs

    # used for testing only..
    @property
    def key(self):
        return (self.codim, self.lhs, self.rhs)

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

    @property
    def identity(self):
        cat = self.cat.supercat
        return cat.Identity(self)

    @property
    def lunitor(self):
        cat = self.cat.supercat
        return cat.LeftUnitor(self)

    @property
    def runitor(self):
        cat = self.cat.supercat
        return cat.RightUnitor(self)

    @property
    def inv(self):
        cat = self.cat.supercat
        return cat.Inv(self)



class Globular(object):
    """
    dim = 0 # set
    dim = 1 # category
    dim = 2 # bicategory
    dim = 3 # tricategory
    """

    DEBUG = False

    def __init__(self, dim, supercat=None):
        # I am a (the) homcat in my supercat
        self.dim = dim
        self.codim = 0 if supercat is None else supercat.codim+1
        self.supercat = supercat # meow
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
        while lhs.dim < rhs.dim:
            lhs = self.Cell(lhs, lhs)
        while rhs.dim < lhs.dim:
            rhs = self.Cell(rhs, rhs)
        assert lhs.dim == rhs.dim
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


class Set(Globular):
    def __init__(self, supercat=None):
        Globular.__init__(self, 0, supercat)

    # put the equate method here? do we care?


class Category(Globular):
    def __init__(self, supercat=None):
        Globular.__init__(self, 1, supercat)
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
        cell = Globular.Object(self, name)
        return cell

    def decorate(self, cell):
        assert isinstance(cell, Cell)
        if cell.is_decorated >= self.dim:
            return
        cell.is_decorated = self.dim
        codim = self.codim
        if cell.dim < codim:
            assert 0, cell.dim
        elif isinstance(cell, Pair):
            pass
        elif cell.dim == codim:
            assert not hasattr(cell, "identity"), "wup"
            cell.identity = Cell(self, cell, cell)
            self.equate(cell.identity*cell.identity, cell.identity)
        elif cell.dim == codim+1:
            src = cell.src.identity
            self.equate(cell, cell*src)
            self.equate(src, src*src)
            tgt = cell.tgt.identity
            self.equate(cell, tgt*cell)
            self.equate(tgt, tgt*tgt)
        else:
            assert 0, cell.dim

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
        if self.DEBUG:
            print("  Category.Pair", self.codim, codim, lhs.dim, lhs, rhs)
        compose = lambda lhs, rhs : Globular.Pair(self, codim, lhs, rhs)
        pair = compose(lhs, rhs)
        if isinstance(rhs, Pair) and rhs.shape == pair.shape:
            # reassoc to the left
            a = lhs
            b, c = rhs
            other = compose(a*b, c)
            self.equate(pair, other)
        elif isinstance(lhs, Pair) and lhs.shape == pair.shape:
            # reassoc to the right
            a, b = lhs
            c = rhs
            bc = b*c
            other = compose(a, bc)
            self.equate(pair, other)
        if self.DEBUG:
            print("  Category.Pair: shape=", pair.shape)
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
    def __init__(self, supercat=None):
        Globular.__init__(self, 2, supercat)

        # we just bundle all the hom categories up into one python object
        self.homcat = Category(self)

    def decorate(self, cell):
        assert isinstance(cell, Cell)
        assert self.codim == 0 # todo
        if cell.is_decorated >= self.dim:
            return
        cell.is_decorated = self.dim
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
            pass
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
        return lunitor

    def RightUnitor(self, cell):
        assert isinstance(cell, Pair)
        assert cell.shape == (0, 1), cell
        lhs, rhs = cell
        runitor = lhs.identity << rhs.runitor
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
        if self.DEBUG:
            print("Bicategory.Pair", self.codim, codim, lhs.dim, lhs, rhs)
        #if codim > 0:
        #    return self.supercat.Pair(codim, lhs, rhs)
        compose = lambda lhs, rhs : Globular.Pair(self, codim, lhs, rhs)
        pair = compose(lhs, rhs)
        self.homcat.decorate(pair)
        self.decorate(pair)
#        if isinstance(rhs, Pair):
#            # reassoc to the left
#            a = lhs
#            b, c = rhs
#            other = compose(a<<b, c)
#            #self.equate(pair, other)
#        elif isinstance(lhs, Pair):
#            # reassoc to the right
#            a, b = lhs
#            c = rhs
#            other = compose(a, b<<c)
#            #self.equate(pair, other)
        # pair.shape == (0,1) or (1,2)
        if self.DEBUG:
            print("Bicategory.Pair: shape=", pair.shape)
        #send = compose(lhs.get_root(), rhs.get_root())
        #self.homcat.equate(pair, send)
        return pair


class Tricategory(Globular):
    def __init__(self, supercat=None):
        Globular.__init__(self, 3, supercat)

        # we just bundle all the hom bicategories up into one python object
        self.homcat = Bicategory(self)


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
    A = Cell(m, l)
    B = Cell(n, m)
    C = Cell(o, n)
    D = Cell(p, o)
    In = n.identity
    Im = m.identity

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
    Cell = cat.Cell

    # 0-cells
    l = Cell(name="l")
    m = Cell(name="m")
    n = Cell(name="n")
    o = Cell(name="o")

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
    A0 = Cell(m, l)
    A1 = Cell(m, l)
    A2 = Cell(m, l)
    A3 = Cell(m, l)
    B0 = Cell(n, m)
    B1 = Cell(n, m)

    assert A0.dim == 1

    assert A0.cat is cat.homcat

    BA = B0<<A0
    assert BA.cat is cat.homcat

    assert A0.is_decorated == 2
    assert hasattr(A0, "identity")
    assert hasattr(A0, "lunitor")
    assert hasattr(A0, "runitor")

    u = A0.identity
    uu = u*u
    assert uu == u

    assert hasattr(BA, "identity")
    assert hasattr(BA, "lunitor")

    assert BA.identity == B0.identity<<A0.identity

    i = BA.identity
    Globular.DEBUG = True
    i*i
    #assert i * i == i # FAIL
    print(i.cat)
    Globular.DEBUG = False

    # 2-cells
    f0 = Cell(A1, A0)
    f1 = Cell(A2, A1)
    f2 = Cell(A3, A2)
    f22 = Cell(A3, A2)
    g0 = Cell(B1, B0)

    assert f0.dim == 2

    assert f2 != f22
    assert f2*f1 != f22*f1

    ff = f1*f0
    assert ff is f1*f0
    assert ff.key not in cat.cache
    assert ff.key in cat.homcat.cache # this Pair lives in the homcat

    assert A1.identity*f0 == f0*A0.identity
    assert (f2*f1)*f0 == f2*(f1*f0)

    gf = g0<<f0
    assert gf is g0<<f0
    assert gf.key in cat.cache

    assert gf.src == B0<<A0
    #assert gf * BA.identity == gf # FAIL


def test_globular():

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

    test_category()
    test_bicategory()
    test_globular()

    print("OK\n")



