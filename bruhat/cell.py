#!/usr/bin/env python3

"""
Render 2-morphisms in a bicategory using sheet/string diagrams.

previous version: huygens/cell.py

"""

import copy
from random import random, randint, shuffle
import operator
from functools import reduce
from math import pi, sin, cos
from time import sleep

from huygens.front import color
#from huygens.sat import System, Listener, Variable
#from huygens.back import Compound, Deco, Path, Transform, Scale
#from huygens.front import path, style, canvas, color, Canvas
#from huygens.back import Visitor
#from huygens.argv import argv
#
#from huygens import pov
#from huygens.pov import View, Mat

EPSILON = 1e-6

black = (0,0,0,1)

def conv(a, b, alpha=0.5):
    return (1.-alpha)*a + alpha*b


"""

The three directions:
    H : horizontal : operator <<  : x coordinate : width property  : left, right
    D : depth      : operator @   : y coordinate : depth property  : front, back
    V : vertical   : operator *   : z coordinate : height property : top, bot

atributes:
              pip_x  pip_y  pip_z
    -ve dir:  left   back   bot
    +ve dir:  right  front  top

"""


class Shape(object):
    def __init__(self, name, weight=1., **kw):
        self.name = name
        self.weight = weight

    def __str__(self):
        return self.name


class Compound(object):
    def __len__(self):
        return len(self.cells)

    def __getitem__(self, idx):
        return self.cells[idx]

    def __contains__(self, other):
        for cell in self.cells:
            if cell is other:
                return True
        return False

    def index(self, other):
        for idx, cell in enumerate(self.cells):
            if cell is other:
                return idx
        assert 0, "not found"

    def _associate(self, items):
        cls = self.__class__
        itemss = [(item.cells if isinstance(item, cls) else [item])
            for item in items]
        items = reduce(operator.add, itemss, [])
        return items

    def visit(self, callback, **kw):
        for child in self.cells:
            child.visit(callback, **kw)
        Shape.visit(self, callback, **kw)

    def deepstr(self, depth=0):
        lines = [Shape.deepstr(self, depth)]
        lines += [cell.deepstr(depth+1) for cell in self.cells]
        return "\n".join(lines)

def setop(cls, opname, parent):
    def meth(left, right):
        return parent([left, right])
    setattr(cls, opname, meth)


# -------------------------------------------------------


class Cell0(Shape):
    "These are the 0-cells"

    color = None # fill
    stroke = black # outline
    show_pip = False

    def __init__(self, name, **kw):
        Shape.__init__(self, name, **kw)
        if "color" in kw:
            self.color = kw["color"]
        if "stroke" in kw:
            self.stroke = kw["stroke"]

    def __len__(self):
        return 1

    def __getitem__(self, idx):
        return [self][idx]



class DCell0(Compound, Cell0):
    def __init__(self, cells, alias=False, **kw):
        cells = self._associate(cells)
        name = "@".join(cell.name for cell in cells) or "ii"
        Cell0.__init__(self, name, **kw)
        self.cells = cells


setop(Cell0, "__matmul__", DCell0)

Cell0.ii = DCell0([])



# -------------------------------------------------------


class Cell1(Shape):
    """
        These are the 1-cells.
    """

    color = black
    show_pip = True
    pip_radius = 0.3

    def __init__(self, tgt, src, name=None, weight=1.0, alias=False, **kw):
        assert isinstance(tgt, Cell0)
        assert isinstance(src, Cell0)
        if name is None:
            name = "(%s<---%s)"%(tgt, src)
        Shape.__init__(self, name, weight, **kw) # will update kw's on __dict__
        self.tgt = tgt
        self.src = src
        self.hom = (self.tgt, self.src)

    def extrude(self, show_pip=False, **kw):
        return Cell2(self, self, show_pip=show_pip, **kw)

    def reassoc(self):
        yield [self]
        

class DCell1(Compound, Cell1):
    def __init__(self, cells, alias=False, **kw):
        cells = self._associate(cells)
        tgt = DCell0([cell.tgt for cell in cells])
        src = DCell0([cell.src for cell in cells])
        name = "(" + "@".join(cell.name for cell in cells) + ")"
        Cell1.__init__(self, tgt, src, name, alias=True, **kw)
        self.cells = cells

    def extrude(self, show_pip=False, **kw):
        cells = [cell.extrude(show_pip=show_pip, **kw) for cell in self.cells]
        return DCell2(cells, show_pip=show_pip, **kw)

    #def reassoc(self):

        

class HCell1(Compound, Cell1):
    def __init__(self, cells, alias=False, **kw):
        cells = self._associate(cells)
        tgt = cells[0].tgt
        src = cells[-1].src
        i = 0
        while i+1 < len(cells):
            if cells[i].src.name != cells[i+1].tgt.name:
                msg = ("can't compose %s and %s"%(cells[i], cells[i+1]))
                raise TypeError(msg)
            i += 1
        name = "(" + "<<".join(cell.name for cell in cells) + ")"
        Cell1.__init__(self, tgt, src, name, alias=True, **kw)
        self.cells = cells

    def extrude(self, show_pip=False, **kw):
        cells = [cell.extrude(show_pip=show_pip, **kw) for cell in self.cells]
        return HCell2(cells, show_pip=show_pip, **kw)


setop(Cell1, "__matmul__", DCell1)
setop(Cell1, "__lshift__", HCell1)

# -------------------------------------------------------



class Cell2(Shape):
    "These are the 2-cells"

    DEBUG = False
    color = black
    show_pip = True
    pip_radius = 0.5
    cone = 0.6 # closer to 1. is more cone-like
    pip_cvs = None

    def __init__(self, tgt, src, name=None, alias=False, **kw):
        assert isinstance(tgt, Cell1)
        assert isinstance(src, Cell1)
        assert tgt.src.name == src.src.name, "%s != %s" % (tgt.src, src.src)
        assert tgt.tgt.name == src.tgt.name, "%s != %s" % (tgt.tgt, tgt.src)
        if name is None:
            name = "(%s<===%s)"%(tgt, src)
        Shape.__init__(self, name, **kw) # update's kw's on __dict__
        self.tgt = tgt
        self.src = src
        self.hom = (self.tgt, self.src)

    @property
    def center(self):
        return self.pip_x, self.pip_y, self.pip_z

    @property
    def rect(self):
        return (
            self.pip_x-self.left,  self.pip_y-self.front, self.pip_z-self.bot,
            self.pip_x+self.right, self.pip_y+self.back, self.pip_z+self.top,
        )

    @property
    def hunits(self):
        return 1

    @property
    def dunits(self):
        return 1

    @property
    def vunits(self):
        return 1

    def vflip(self, alias=False):
        tgt, src = self.src, self.tgt
        return Cell2(tgt, src, alias=alias)


class DCell2(Compound, Cell2):
    def __init__(self, cells, alias=False, **kw):
        cells = self._associate(cells)
        tgt = DCell1([cell.tgt for cell in cells])
        src = DCell1([cell.src for cell in cells])
        name = "(" + "@".join(cell.name for cell in cells) + ")"
        Cell2.__init__(self, tgt, src, name, alias=True, **kw)
        self.cells = cells

    @property
    def hunits(self):
        return max(cell.hunits for cell in self.cells)

    @property
    def dunits(self):
        return sum(cell.dunits for cell in self.cells)

    @property
    def vunits(self):
        return max(cell.vunits for cell in self.cells)

    def vflip(self, alias=False):
        cells = [cell.vflip(alias) for cell in self.cells]
        return DCell2(cells, alias)


class HCell2(Compound, Cell2):
    def __init__(self, cells, alias=False, **kw):
        cells = self._associate(cells)
        tgt = HCell1([cell.tgt for cell in cells])
        src = HCell1([cell.src for cell in cells])
        name = "(" + "<<".join(cell.name for cell in cells) + ")"
        Cell2.__init__(self, tgt, src, name, alias=True, **kw)
        self.cells = cells

    @property
    def hunits(self):
        return sum(cell.hunits for cell in self.cells)

    @property
    def dunits(self):
        return max(cell.dunits for cell in self.cells)

    @property
    def vunits(self):
        return max(cell.vunits for cell in self.cells)

    def vflip(self, alias=False):
        cells = [cell.vflip(alias) for cell in self.cells]
        return HCell2(cells, alias)


class VCell2(Compound, Cell2):
    def __init__(self, cells, alias=False, **kw):
        cells = self._associate(cells)
        tgt = cells[0].tgt
        src = cells[-1].src
        i = 0
        while i+1 < len(cells):
            if cells[i].src.name != cells[i+1].tgt.name:
                msg = ("can't compose\n%s and\n%s"%(cells[i], cells[i+1]))
                #raise TypeError(msg)
                #print("VCell2.__init__: WARNING", msg)
            i += 1
        name = "(" + "*".join(cell.name for cell in cells) + ")"
        Cell2.__init__(self, tgt, src, name, alias=True, **kw)
        self.cells = cells

    @property
    def hunits(self):
        return max(cell.hunits for cell in self.cells)

    @property
    def dunits(self):
        return max(cell.dunits for cell in self.cells)

    @property
    def vunits(self):
        return sum(cell.vunits for cell in self.cells)

    def vflip(self, alias=False):
        cells = [cell.vflip(alias) for cell in reversed(self.cells)]
        return VCell2(cells, alias)


setop(Cell2, "__matmul__", DCell2)
setop(Cell2, "__lshift__", HCell2)
setop(Cell2, "__mul__", VCell2)

# -------------------------------------------------------


def test():
    l = Cell0("l")
    m = Cell0("m")
    n = Cell0("n")
    o = Cell0("o")
    p = Cell0("p")

    assert str(m@n) == "m@n", str(m@n)
    #assert m != n
    #assert m@m == m@m

    ii = Cell0.ii
    mm = m@m
    mmm = m@m@m
    mmmm = m@m@m@m
    A = Cell1(mm, mm)
    B = Cell1(m, m) @ Cell1(m, m)
    AA = A<<A

    f = Cell2(B, B)

    #cell = Cell2(B, AA) * Cell2(AA, A)
    #cell = cell @ (f*f*f)
    #cell = Cell2(B, A<<B)

    cell = Cell1(mm,mm)<<((Cell1(m,mmm)) @ Cell1(m,m)) << Cell1(mmmm,m)
    cell = Cell2(cell, Cell1(mm,m))

    mm_ = Cell1(mm, ii)
    m_m = Cell1(m, m)
    mm_m = Cell1(mm, m)
    _mm = Cell1(ii, mm)
    _m = Cell1(ii, m)
    cell = Cell2(_mm, _mm) << Cell2(mm_, mm_)
    cell = Cell2(m_m, (m_m @ _m) << mm_m)

    A, B = Cell1(m@n, l), Cell1(l, p)
    A1, B1 = Cell1(m@n, o@o), Cell1(o@o, p)
    A2, B2 = Cell1(m@n, p@l), Cell1(p@l, p)

    AB = A << B
    str(A@A)
    str(AB)

    f = Cell2(A, A)
    g = Cell2(B, B)
    str(f)
    str(f@f)
    str(f << g)

    f = Cell2(A<<B, A1<<B1)
    g = Cell2(A1<<B1, A2<<B2)
    str(f*g)


def more_test():

    scheme = "ff5e5b-d8d8d8-ffffea-00cecb-ffed66"
    scheme = scheme.split("-")
    scheme = [color.rgbhex(rgb).alpha(0.5) for rgb in scheme]

    names = 'lmnop'
    l, m, n, o, p = [
        Cell0(name, color=scheme[i%len(scheme)], address=name) 
        for i, name in enumerate('lmnop')]
    i0 = Cell0("i", color=None)

    I_l = Cell1(l, l, show_pip=False, color=None)
    I_m = Cell1(m, m, show_pip=False, color=None)
    I_n = Cell1(n, n, show_pip=False, color=None)
    I_o = Cell1(o, o, show_pip=False, color=None)

    cell = Cell1(m, m@m) << Cell1(m@m, n)
    cell = cell @ cell

    # i don't think we can handle bubbles... 
    ii = Cell0.ii
    bubble = Cell1(ii, m@m) << Cell1(m@m, ii)

    cell = bubble

    cell = bubble<<bubble

    cell = bubble@bubble

    cell = Cell1(m, m) @ bubble


    l_mn = Cell1(l, m@n)
    mn_l = Cell1(m@n, l)
    o_mn = Cell1(o, m@n)
    mn_o = Cell1(m@n, o)
    l_l = Cell1(l, l)
    o_o = Cell1(o, o)
    left = (o_mn << mn_o << o_o) @ (l_mn << mn_l << l_l)
    right = (o_mn @ l_mn) << (mn_o @ mn_l) << (o_o @ l_l)

    top = left.extrude()
    #bot = right.extrude()
    #cell = top * bot
    cell = top

    o_mn = Cell1(o, m@n)
    cell = Cell2(I_o<<o_mn, o_mn<<(I_m@I_n), cone=1.)


    swap = lambda a,b : Cell1(a@b, b@a)
    # Yang-Baxter
    def yang_baxter(n, m, l, reversed=False, **kw):
        I_l = Cell1(l, l, show_pip=False, color=None)
        I_m = Cell1(m, m, show_pip=False, color=None)
        I_n = Cell1(n, n, show_pip=False, color=None)
        tgt = (I_n @ swap(m, l)) << (swap(n, l) @ I_m) << (I_l @ swap(n, m))
        src = (swap(n, m) @ I_l) << (I_m @ swap(n, l)) << (swap(m, l) @ I_n)
        if reversed:
            tgt, src = src, tgt
        morph = Cell2(tgt, src, cone=1.0, show_pip=False, **kw)
        return morph
    
    #tgt = (I_n @ I_m ) << S_nm
    #src = S_nm << (I_m @ I_n)
    #morph = Cell2(tgt, src)
    
    # a part of the Zamolodchikov Tetrahedron Equation
    lhs = (I_o @ ((I_n@swap(m,l))<<(swap(n,l)@I_m)) ).extrude(show_pip=False)
    rhs = (I_l @ ((swap(o,m)@I_n)<<(I_m@swap(o,n)) )).extrude(show_pip=False)
    
    tgt = swap(n,m) << (I_m @ I_n)
    src = (I_n @ I_m) << swap(n,m)
    back = Cell2(tgt, src, show_pip=False, cone=1.0)
    tgt = (I_o @ I_l) << swap(o,l)
    src = swap(o,l) << (I_l @ I_o)
    front = Cell2(tgt, src, show_pip=False, cone=1.0)
    
    morph_0 = lhs << (front @ back) << rhs
    
    rhs = I_l.extrude(show_pip=False) @ yang_baxter(o, n, m, reversed=True)
    lhs = (I_o @ I_n @ swap(m,l)) << (I_o @ swap(n,l) @ I_m) << (swap(o,l) @ I_n @ I_m)
    lhs = lhs.extrude()
    morph_1 = lhs << rhs
    morph_1 = morph_1.vflip()
    
            


if __name__ == "__main__":
    print("\n")
    test()
    more_test()

    print("OK")




