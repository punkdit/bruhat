#!/usr/bin/env python3

from random import random, choice, seed

from math import pi


from bruhat.render.sat import Expr, Variable, System
from bruhat.render.back import RGBA, Canvas, Scale
from bruhat.render.back import path, style
from bruhat.argv import argv


class Box(object):

    DEBUG = False
    did_layout = False

    @classmethod
    def promote(cls, item):
        if isinstance(item, Box):
            return item
        if isinstance(item, str):
            return TextBox(item)
        if isinstance(item, (tuple, list)):
            return HBox(item)
        raise TypeError(repr(item))

    def on_layout(self, cvs, system):
        assert not self.did_layout, "already called on_layout"
        if self.DEBUG:
            print("%s.on_layout" % (self.__class__.__name__,))
        for attr in 'x y left right top bot'.split():
            if attr in self.__dict__:
                continue
            stem = self.__class__.__name__ + '.' + attr

            # We don't try to minimize the absolute coordinate values.
            weight = 1.0 if attr not in 'xy' else 0.0
            vmin = None if attr in 'xy' else 0.
            v = system.get_var(stem, weight, vmin=vmin)
            setattr(self, attr, v)
        self.did_layout = True

    def on_render(self, cvs, system):
#        x = system[self.x]
#        y = system[self.y]
#        left = system[self.left]
#        right = system[self.right]
#        top = system[self.top]
#        bot = system[self.bot]
        attrs = list(self.__dict__.keys())
        for attr in attrs:
            value = getattr(self, attr)
            if not isinstance(value, Expr):
                continue
            value = system[value]
            setattr(self, attr, value)
        if not self.DEBUG:
            return
        x = self.x
        y = self.y
        left = self.left
        right = self.right
        top = self.top
        bot = self.bot
        #cvs.set_line_width(0.5)
        cl = RGBA(1., 0., 0., 0.5)
        r = 0.1
        cvs.stroke(path.line(x-r, y-r, x+r, y+r), [cl]) #, style.linewidth.Thick])
        cvs.stroke(path.line(x+r, y-r, x-r, y+r), [cl])
        #bg = RGBA(0.5*random(), 0.5*random(), 0.5*random(), 0.5)
        bg = RGBA(0.5, 0.5, 0., 0.1)
        cvs.fill(path.rect(x-left, y-bot, left+right, top+bot), [bg])
        cvs.stroke(path.rect(x-left, y-bot, left+right, top+bot), [cl])

    @property
    def width(self):
        return self.left + self.right

    @property
    def height(self):
        return self.top + self.bot

    @property
    def llx(self):
        return self.x - self.left

    @property
    def lly(self):
        return self.y - self.bot

    @property
    def urx(self):
        return self.x + self.right

    @property
    def ury(self):
        return self.y + self.top

    @property
    def bound(self):
        return self.llx, self.lly, self.urx, self.ury

    def get_align(self, align):
        llx, lly, urx, ury = self.bound
        xmid = 0.5*(llx + urx)
        ymid = 0.5*(lly + ury)
        if align == "center":
            x, y = xmid, ymid
        elif align == "north":
            x, y = xmid, ury
        elif align == "south":
            x, y = xmid, lly
        elif align == "east":
            x, y = urx, ymid
        elif align == "west":
            x, y = llx, ymid
        elif align == "northeast":
            x, y = urx, ury
        elif align == "northwest":
            x, y = llx, ury
        elif align == "southeast":
            x, y = urx, lly
        elif align == "southwest":
            x, y = llx, lly
        else:
            assert 0, "alignment %r not understood" % align
    
        return x, y


    def render(self, cvs, x=0, y=0):
        system = System()
        self.on_layout(cvs, system)
        system.add(self.x == x)
        system.add(self.y == y)
        system.solve()
        self.on_render(cvs, system)


class EmptyBox(Box):
    def __init__(self, top=0., bot=0., left=0., right=0.):
        self.top = top
        self.bot = bot
        self.left = left
        self.right = right


class EmptyBox(Box):
    def __init__(self, top=None, bot=None, left=None, right=None):
        if top is not None:
            self.top = top
        if bot is not None:
            self.bot = bot
        if left is not None:
            self.left = left
        if right is not None:
            self.right = right


    
class StrokeBox(Box):
    def __init__(self, width, height, rgba=(0., 0., 0., 1.)):
        self.top = height
        self.bot = 0.
        self.left = 0.
        self.right = width
        self.rgba = rgba

    def on_render(self, cvs, system):
        Box.on_render(self, cvs, system)
        x = system[self.x]
        y = system[self.y]
        cvs.stroke(path.rect(x, y, self.width, self.height), [RGBA(*self.rgba)])


class FillBox(StrokeBox):
    def on_render(self, cvs, system):
        Box.on_render(self, cvs, system)
        x = system[self.x]
        y = system[self.y]
        cvs.fill(path.rect(x, y, self.width, self.height), [RGBA(*self.rgba)])


class TextBox(Box):
    def __init__(self, text):
        self.text = text

    def on_layout(self, cvs, system):
        Box.on_layout(self, cvs, system)
        extents = cvs.text_extents(self.text)
        dx, dy, width, height = extents
        system.add(self.left + self.right == width+dx)
        system.add(self.top + self.bot == height, 1.0)
        system.add(self.left == 0)
        assert dy >= 0., dy
        system.add(self.top == dy)

    def on_render(self, cvs, system):
        Box.on_render(self, cvs, system)
        x = system[self.x]
        y = system[self.y]
        cvs.text(x, y, self.text)


class ChildBox(Box):
    def __init__(self, child):
        self.child = Box.promote(child)

    def on_render(self, cvs, system):
        Box.on_render(self, cvs, system)
        self.child.on_render(cvs, system)


class MarginBox(ChildBox):
    def __init__(self, child, margin):
        ChildBox.__init__(self, child)
        self.margin = margin

    def SLOW_on_layout(self, cvs, system):
        Box.on_layout(self, cvs, system)
        child = self.child
        child.on_layout(cvs, system)
        system.add(self.x == child.x)
        system.add(self.y == child.y)
        margin = self.margin
        system.add(self.left == child.left + margin)
        system.add(self.right == child.right + margin)
        system.add(self.top == child.top + margin)
        system.add(self.bot == child.bot + margin)

    def on_layout(self, cvs, system):
        child = self.child
        child.on_layout(cvs, system)
        self.x = child.x
        self.y = child.y
        margin = self.margin
        self.left = child.left + margin
        self.right = child.right + margin
        self.top = child.top + margin
        self.bot = child.bot + margin
        Box.on_layout(self, cvs, system)


class AlignBox(ChildBox):
    def __init__(self, child, align):
        ChildBox.__init__(self, child)
        self.align = align
        assert align in ("center north south east west northeast"+
            "northwest southeast southwest").split(), "align %r not understood"%align

    def on_layout(self, cvs, system):
        child = self.child
        child.on_layout(cvs, system)
        x, y = child.get_align(self.align)
        self.x = x
        self.y = y
        self.left = x - child.llx
        self.right = child.urx - x
        self.bot = y - child.lly
        self.top = child.ury - y
        Box.on_layout(self, cvs, system)


class CompoundBox(Box):
    def __init__(self, boxs, weight=None):
        assert len(boxs)
        self.boxs = [Box.promote(box) for box in boxs]
        self.weight = weight

    def on_layout(self, cvs, system):
        Box.on_layout(self, cvs, system)
        for box in self.boxs:
            box.on_layout(cvs, system)

    def on_render(self, cvs, system):
        Box.on_render(self, cvs, system)
        for box in self.boxs:
            box.on_render(cvs, system)


class OBox(CompoundBox):
    "Overlay boxes on top of each other, with matching anchors"
    strict = False
    def on_layout(self, cvs, system):
        CompoundBox.on_layout(self, cvs, system)
        boxs = self.boxs
        for box in boxs:
            system.add(self.x == box.x) # align
            system.add(self.y == box.y) # align
            if self.strict:
                system.add(box.left == self.left, self.weight)
                system.add(box.right == self.right, self.weight)
                system.add(box.top == self.top, self.weight)
                system.add(box.bot == self.bot, self.weight)
            else:
                system.add(box.left <= self.left)
                system.add(box.right <= self.right)
                system.add(box.top <= self.top)
                system.add(box.bot <= self.bot)


class HBox(CompoundBox):
    "horizontal compound box: anchor left"
    strict = False
    def on_layout(self, cvs, system):
        CompoundBox.on_layout(self, cvs, system)
        boxs = self.boxs
        system.add(self.left == 0.) # left anchor
        left = self.x
        for box in boxs:
            system.add(self.y == box.y) # align
            system.add(box.x - box.left == left)
            left += box.width
            if self.strict:
                system.add(box.top == self.top, self.weight)
                system.add(box.bot == self.bot, self.weight)
            else:
                system.add(box.top <= self.top)
                system.add(box.bot <= self.bot)
        system.add(self.x + self.width == left)

class StrictHBox(HBox):
    strict = True

class VBox(CompoundBox):
    "vertical compound box: anchor top"
    strict = False
    def on_layout(self, cvs, system):
        CompoundBox.on_layout(self, cvs, system)
        boxs = self.boxs
        system.add(self.top == 0.) # top anchor
        y = self.y
        for box in boxs:
            system.add(self.x == box.x) # align
            system.add(box.y + box.top == y)
            y -= box.height
            if self.strict:
                system.add(box.left == self.left, self.weight)
                system.add(box.right == self.right, self.weight)
            else:
                system.add(box.left <= self.left)
                system.add(box.right <= self.right)
        system.add(self.y - self.bot == y)

class StrictVBox(VBox):
    strict = True


class TableBox(CompoundBox):
    def __init__(self, rows, grid=False):
        boxs = []
        m = len(rows) # rows
        n = len(rows[0]) # cols
        assert n>0
        for row in rows:
            assert len(row) == n
            for box in row:
                box = Box.promote(box)
                boxs.append(box)
        self.rows = [list(row) for row in rows]
        self.shape = m, n
        # anchor is top left
        self.top = 0.
        self.left = 0.
        self.grid = grid
        CompoundBox.__init__(self, boxs)

    def on_layout(self, cvs, system):
        CompoundBox.on_layout(self, cvs, system)
        m, n = self.shape
        rows = self.rows
        xs, ys = {}, {}
        ws, hs = {}, {} # width's, height's
        for i in range(m): # rows
            ys[i] = system.get_var("TableBox.row(%d)"%i, weight=0.)
            hs[i] = system.get_var("TableBox.height(%d)"%i, weight=1.) # minimize
        for j in range(n): # cols
            xs[j] = system.get_var("TableBox.col(%d)"%j, weight=0.)
            ws[j] = system.get_var("TableBox.width(%d)"%j, weight=1.) # minimize
        for i in range(m): # rows
            for j in range(n): # cols
                box = rows[i][j]
                system.add(box.y == ys[i]) # align
                system.add(box.x == xs[j]) # align

        for i in range(m): # rows
            x = self.x
            for j in range(n): # cols
                box = rows[i][j]
                system.add(box.x - box.left >= x)
                width = ws[j] # width of this col
                x += width
                system.add(box.x + box.right <= x)
            system.add(self.x + self.width >= x)

        for j in range(n): # cols
            y = self.y
            for i in range(m): # rows
                box = rows[i][j]
                system.add(box.y + box.top <= y)
                height = hs[i]
                y -= height
                system.add(box.y - box.bot >= y)
            system.add(self.y - self.height <= y)
        self.vs = xs, ys, ws, hs

    def on_render(self, cvs, system):
        CompoundBox.on_render(self, cvs, system)
        if not self.grid:
            return
        m, n = self.shape
        xs, ys, ws, hs = self.vs
        width = system[self.width]
        height = system[self.height]
        x = system[self.x]
        y = system[self.y]
        for j in range(n):
            cvs.stroke(path.line(x, y, x, y-height))
            x += system[ws[j]]
        #cvs.stroke(path.line(x, y, x, y-height))
        x = system[self.x]
        y = system[self.y]
        for i in range(m):
            cvs.stroke(path.line(x, y, x+width, y))
            y -= system[hs[i]]
        x = system[self.x]
        y = system[self.y]
        cvs.stroke(path.rect(x, y-height, width, height))



def test_build():

    def rnd(a, b):
        return (b-a)*random() + a


    box = EmptyBox(1., 1., 1., 1.)
    yield box, "empty"

    box = TextBox("Hey there!")
    #box = TextBox(".")
    yield box, "text"

    box = HBox("geghh xxde xyeey".split())
    yield box, "hbox-text"

    box = VBox("geghh xxde xyeey".split())
    yield box, "vbox-text"


    box = OBox([
        EmptyBox(.4, .1, 0., 2.2),
        EmptyBox(.3, 0., .5, 2.5),
        EmptyBox(1., .5, .5, .5),
        FillBox(.2, .2),
    ])
    yield box, "obox"


    box = HBox([
        VBox([TextBox(text) for text in "xxx1 ggg2 xxx3 xx4".split()]),
        VBox([TextBox(text) for text in "123 xfdl sdal".split()]),
    ])
    yield box, "hbox-vbox"


    box = TableBox([
        [EmptyBox(.4, .1, 0.2, 2.2), EmptyBox(.3, 1.2, .5, 2.5),],
        [EmptyBox(.8, .1, 0.4, 1.2), EmptyBox(.5, 0.4, .5, 1.5),]
    ])
    yield box, "table"


    a, b = 0.2, 1.0
    rows = []
    for row in range(3):
        row = []
        for col in range(3):
            box = EmptyBox(rnd(a,b), rnd(a,b), rnd(a,b), rnd(a,b))
            row.append(box)
        rows.append(row)

    box = TableBox(rows)
    yield box, "table-2"


    rows = []
    for i in range(3):
        row = []
        for j in range(3):
            box = TextBox(choice("xbcgef")*(i+1)*(j+1))
            box = MarginBox(box, 0.1)
            box = AlignBox(box, "north")
            row.append(box)
        row.append(EmptyBox(bot=1.))
        rows.append(row)
    box = TableBox(rows)
    yield box, "table-3"


    a, b = 0.2, 2.0
    rows = []
    for row in range(2):
        boxs = []
        for col in range(2):
            top = rnd(a, b)
            bot = rnd(a, b)
            left = rnd(a, b)
            right = rnd(a, b)
            box = EmptyBox(top, bot, left, right)
            boxs.append(box)
        boxs.append(TextBox("Hig%d !"%row))
        #boxs.append(TextBox(r"$\to$"))
        box = HBox(boxs)
        rows.append(box)
    box = VBox(rows)
    yield box, "table-4"


def test():

    #Box.DEBUG = argv.get("debug") or argv.get("DEBUG")
    Box.DEBUG = True
    EmptyBox.DEBUG = True

    for (box, name) in test_build():
    
        try:
            print("rendering", name)
            cvs = Canvas()
            cvs.append(Scale(2.0))
            box.render(cvs)
            #cvs = Canvas([Scale(2.0, 2.0)]+cvs.items)
            cvs.writePDFfile("doc/pic-%s.pdf" % name)
            cvs.writeSVGfile("doc/pic-%s.svg" % name)
            print()
        except:
            print("render failed for %r"%name)
            raise


if __name__ == "__main__":

    seed(0)
    test()

    print("OK\n")


