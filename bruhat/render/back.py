#!/usr/bin/env python3

"""
Expose a _drawing api with multiple (hopefully) backends.

The api uses the first quadrant coordinate system:

        y 
        ^
        |
        |
        | positive x and y
        |    
        |
      --+-------------------> x
        |
    
"""

from math import pi, sqrt, sin, cos, sqrt, floor

from bruhat.render.base import SCALE_CM_TO_POINT, Base, Matrix
from bruhat.render.text import make_text


# ----------------------------------------------------------------------------
# 
#

def n_min(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return min(a, b)

def n_max(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return max(a, b)


class Bound(Base):
    def __init__(self, llx=None, lly=None, urx=None, ury=None):
        assert urx is None or llx is None or urx-llx >= 0., (llx, lly, urx, ury)
        assert ury is None or lly is None or ury-lly >= 0., (llx, lly, urx, ury)
        self.llx = llx
        self.lly = lly
        self.urx = urx
        self.ury = ury

    def scale_point_to_cm(self):
        llx = self.llx / SCALE_CM_TO_POINT
        lly = self.lly / SCALE_CM_TO_POINT
        urx = self.urx / SCALE_CM_TO_POINT
        ury = self.ury / SCALE_CM_TO_POINT
        return Bound(llx, lly, urx, ury)

    def union(self, other):
        "union"
        llx = n_min(self.llx, other.llx)
        lly = n_min(self.lly, other.lly)
        urx = n_max(self.urx, other.urx)
        ury = n_max(self.ury, other.ury)
        return Bound(llx, lly, urx, ury)
    __add__ = union

    def update(self, other):
        "union update"
        self.llx = n_min(self.llx, other.llx)
        self.lly = n_min(self.lly, other.lly)
        self.urx = n_max(self.urx, other.urx)
        self.ury = n_max(self.ury, other.ury)

    def nonempty(self):
        return (self.llx is not None or self.lly is not None or 
            self.urx is not None or self.ury is not None)

    def is_empty(self):
        return (self.llx is None and self.lly is None and 
            self.urx is None and self.ury is None)

    def __getitem__(self, idx):
        return (self.llx, self.lly, self.urx, self.ury)[idx]

    @property
    def width(self):
        return self.urx - self.llx

    @property
    def height(self):
        return self.ury - self.lly




# ----------------------------------------------------------------------------
# 
#


class Visitor(object):
    def on_item(self, item):
        pass

class DumpVisitor(Visitor):
    def on_item(self, item):
        print('\t%s'%item)


# XXX this probably should be a full fledged Context XXX
class BoundVisitor(Visitor):
    def __init__(self):
        self.pos = None
        self.lw = _defaultlinewidth*SCALE_CM_TO_POINT
        self.bound = Bound()

    def on_item(self, item):
        if isinstance(item, (Scale, Rotate, Translate_Pt)):
            assert 0, "trafo %s not implemented"%item
        elif isinstance(item, Compound):
            assert 0, "%s: save restore not implemented"%item
        elif isinstance(item, MoveTo_Pt):
            self.pos = item.x, item.y
        elif isinstance(item, (Stroke, Fill)):
            self.pos = None
        elif isinstance(item, LineWidth_Pt):
            self.lw = item.lw
        elif isinstance(item, (LineTo_Pt, CurveTo_Pt)):
            b = item.get_bound()
            if b.is_empty():
                assert 0
            # add a generous linewidth... does not work very well... XXX
            r = 0.5*self.lw
            b.llx -= r
            b.lly -= r
            b.urx += r
            b.ury += r
            self.bound.update(b)
            if self.pos is not None:
                x, y = self.pos
                self.bound.update(Bound(x, y, x, y))
        else:
            pass
        #print("BoundVisitor.on_item", self.bound)



# ----------------------------------------------------------------------------
# 
#


class Item(Base):

    DEBUG = False

    def get_bound(self):
        return Bound()

    def process_cairo(self, cxt):
        pass

    def visit(self, visitor):
        visitor.on_item(self)

    def dump(self):
        visitor = DumpVisitor()
        self.visit(visitor)



class ClosePath(Item):
    def process_cairo(self, cxt):
        cxt.close_path()


class MoveTo_Pt(Item):
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def get_bound(self):
        return Bound(self.x, self.y, self.x, self.y)

    def process_cairo(self, cxt):
        if self.DEBUG:
            print("ctx.move_to", self.x, self.y)
        cxt.move_to(self.x, -self.y)


class MoveTo(MoveTo_Pt):
    def __init__(self, x, y):
        x, y = float(x), float(y)
        self.x = SCALE_CM_TO_POINT*x
        self.y = SCALE_CM_TO_POINT*y


class LineTo_Pt(Item):
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def get_bound(self):
        return Bound(self.x, self.y, self.x, self.y)

    def process_cairo(self, cxt):
        if self.DEBUG:
            print("ctx.line_to", self.x, self.y)
        cxt.line_to(self.x, -self.y)


class LineTo(LineTo_Pt):
    def __init__(self, x, y):
        x, y = float(x), float(y)
        self.x = SCALE_CM_TO_POINT*x
        self.y = SCALE_CM_TO_POINT*y


class CurveTo_Pt(Item):
    def __init__(self, x0, y0, x1, y1, x2, y2):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2

    def get_bound(self):
        x0 = min([self.x0, self.x1, self.x2])
        x1 = max([self.x0, self.x1, self.x2])
        y0 = min([self.y0, self.y1, self.y2])
        y1 = max([self.y0, self.y1, self.y2])
        return Bound(x0, y0, x1, y1)

    def process_cairo(self, cxt):
        if self.DEBUG:
            print("ctx.curve_to", self.x0, self.y0, self.x1, self.y1, self.x2, self.y2)
        cxt.curve_to(self.x0, -self.y0, self.x1, -self.y1, self.x2, -self.y2)


class CurveTo(CurveTo_Pt):
    def __init__(self, x0, y0, x1, y1, x2, y2):
        self.x0 = SCALE_CM_TO_POINT*x0
        self.y0 = SCALE_CM_TO_POINT*y0
        self.x1 = SCALE_CM_TO_POINT*x1
        self.y1 = SCALE_CM_TO_POINT*y1
        self.x2 = SCALE_CM_TO_POINT*x2
        self.y2 = SCALE_CM_TO_POINT*y2


class Arc_Pt(Item):
    def __init__(self, x, y, r, angle1, angle2):
        "angle in radians"
        self.x = x
        self.y = y
        self.r = r
        self.angle1 = angle1
        self.angle2 = angle2

    def get_bound(self):
        r = self.r
        return Bound(self.x-r, self.y-r, self.x+r, self.y+r) # XXX TODO XXX

    def process_cairo(self, cxt):
        cxt.arc_negative(self.x, -self.y, self.r, 2*pi-self.angle1, 2*pi-self.angle2)


class Arc(Arc_Pt):
    def __init__(self, x, y, r, angle1, angle2):
        "angle in radians"
        self.x = SCALE_CM_TO_POINT*x
        self.y = SCALE_CM_TO_POINT*y
        self.r = SCALE_CM_TO_POINT*r
        self.angle1 = angle1
        self.angle2 = angle2


class _ArcnMixin(object):
    def process_cairo(self, cxt):
        cxt.arc(self.x, -self.y, self.r, 2*pi-self.angle1, 2*pi-self.angle2)


class Arcn(Arc, _ArcnMixin):
    pass


class Arcn_Pt(Arc_Pt, _ArcnMixin):
    pass


# ----------------------------------------------------------------------------
# 
#

class Polygon(Item):
    def __init__(self, pts, fill=None, stroke=None):
        Item.__init__(self)
        assert len(pts)>1
        self.pts = [(x*SCALE_CM_TO_POINT, y*SCALE_CM_TO_POINT) for (x, y) in pts]
        self.fill = fill
        self.stroke = stroke

    def process_cairo(self, cxt):
        pts = self.pts
        cxt.save() # <--------- save
        #cxt.set_line_width(2.0)

        fill = self.fill #or (0., 0., 0., 1.)
        if fill is not None:
            cxt.set_source_rgba(*fill)
            x, y = pts[0]
            cxt.move_to(x, -y)
            for (x, y) in pts[1:]:
                cxt.line_to(x, -y)
            cxt.close_path()
            cxt.fill()
    
        stroke = self.stroke # or (1., 1., 1., 1.)
        if stroke is not None:
            cxt.set_source_rgba(*stroke)
            x, y = pts[0]
            cxt.move_to(x, -y)
            for (x, y) in pts[1:]:
                cxt.line_to(x, -y)
            cxt.close_path()
            cxt.stroke()

        cxt.restore() # <------- restore


class Polymesh(Item):
    def __init__(self, pts, fills):
        Item.__init__(self)
        assert len(pts)>=3
        assert len(pts) == len(fills)
        self.pts = [(x*SCALE_CM_TO_POINT, y*SCALE_CM_TO_POINT) for (x, y) in pts]
        self.fills = fills

    def process_cairo(self, cxt):
        import cairo
        pts = self.pts
        fills = self.fills
        cxt.save() # <--------- save

        pts = list(self.pts)
        if len(pts)==3:
            pts.append(pts[-1])
            fills.append(fills[-1])
        assert len(pts)==4, len(pts)
        m = cairo.MeshPattern()
        m.begin_patch()
        x, y = pts[0]
        m.move_to(x, -y)
        for (x,y) in pts[1:]:
            m.line_to(x, -y)
        for i in range(4):
            m.set_corner_color_rgba(i, *fills[i])
        m.end_patch()

        cxt.set_source(m)
        x, y = pts[0]
        cxt.move_to(x, -y)
        for (x, y) in pts[1:]:
            cxt.line_to(x, -y)
        cxt.close_path()
        cxt.fill()

        cxt.restore() # <------- restore


class Ball(Item):
    def __init__(self, x, y, radius, rgb0=None, rgb1=None):
        Item.__init__(self)
        self.x = x*SCALE_CM_TO_POINT
        self.y = y*SCALE_CM_TO_POINT
        self.radius = radius*SCALE_CM_TO_POINT
        self.rgb0 = rgb0
        self.rgb1 = rgb1

    def process_cairo(self, cxt):
        x, y, radius = self.x, -self.y, self.radius
        cxt.save()
        cxt.set_line_width(0.5)
        cxt.arc(x, y, radius, 0., pi*2)
        cxt.stroke()
        cx0, cy0 = cx1, cy1 = x+0.8*radius, y+0.8*radius
        radius0 = 0.2*radius
        radius1 = 2.0*radius
        p = cairo.RadialGradient(cx0, cy0, radius0, cx1, cy1, radius1)
        if self.rgb0 is not None:
            p.add_color_stop_rgba(0, *self.rgb0)
        else:
            p.add_color_stop_rgba(0, 0.9, 0.9, 0.9, 1)
        if self.rgb1 is not None:
            p.add_color_stop_rgba(1, *self.rgb1)
        else:
            p.add_color_stop_rgba(1, 0.6, 0.6, 0.6, 1.)
        cxt.set_source(p)
        cxt.arc(x, y, radius, 0., pi*2)
        cxt.fill()
        cxt.restore()


# ----------------------------------------------------------------------------
#
#

class Compound(Item):
    def __init__(self, *args):
        items = []
        for arg in args:
            if type(arg) is list:
                items += arg
            else:
                items.append(arg)
        self.items = items

    def append(self, item):
        items = self.items
        items.append(item)

    def __len__(self):
        return len(self.items)

    def __getitem__(self, idx):
        return self.items[idx]

    def visit(self, visitor):
        for item in self.items:
            item.visit(visitor)

    def get_bound(self):
        visitor = BoundVisitor()
        self.visit(visitor)
        return visitor.bound

    def get_bound_cairo(self):
        # Not as tight as it could be.... ?
        import cairo
        surface = cairo.RecordingSurface(cairo.Content.COLOR_ALPHA, None)
        cxt = cairo.Context(surface)
        self.process_cairo(cxt)
        extents = surface.ink_extents()
        (ulx, uly, width, height) = extents
        llx, lly = ulx, -uly-height
        urx, ury = llx+width, lly+height
        return Bound(llx, lly, urx, ury)

    def process_cairo(self, cxt):
        cxt.save()
        for item in self.items:
            item.process_cairo(cxt)
        cxt.restore()


class Path(Compound):
    def reversed(self, pos=None):
        if len(self)==0:
            return self
        #print("Path.reversed")
        #for item in self:
        #    print('\t', item)
        n = len(self)
        assert n%2 == 0, "not implemented"
        items = []
        idx = n-2
        while idx>=0:
            item = self[idx]
            assert isinstance(item, LineTo_Pt), "not implemented"
            x, y = item.x, item.y
            item = self[idx+1]
            assert isinstance(item, CurveTo_Pt), "not implemented"
            x0, y0 = item.x0, item.y0
            x1, y1 = item.x1, item.y1
            x2, y2 = item.x2, item.y2
            items.append(LineTo_Pt(x2, y2))
            items.append(CurveTo_Pt(x1, y1, x0, y0, x, y))
            idx -= 2
        p = Path(items)
        return p



class Line(Path):
    def __init__(self, x0, y0, x1, y1):
        Path.__init__(self, [
            MoveTo(x0, y0), 
            LineTo(x1, y1)])


class Curve(Path):
    def __init__(self, x0, y0, x1, y1, x2, y2, x3, y3):
        Path.__init__(self, [
            MoveTo(x0, y0), 
            CurveTo(x1, y1, x2, y2, x3, y3)])
    

class Rect(Path):
    def __init__(self, x, y, width, height):
        Path.__init__(self, [
            MoveTo(x, y), 
            LineTo(x+width, y),
            LineTo(x+width, y+height),
            LineTo(x, y+height),
            ClosePath()])


class Circle(Path):
    def __init__(self, x, y, r):
        Path.__init__(self, [
            MoveTo(x+r, y),
            Arc(x, y, r, 0, 2*pi),
            ClosePath()])
        


# ----------------------------------------------------------------------------
# Deco
#


class Deco(Item):
    pass


class CompoundDeco(Deco):
    def __init__(self, decos):
        self.decos = list(decos)

    def process_cairo(self, cxt):
        for deco in self.decos:
            deco.process_cairo(cxt)


class Stroke(Deco):
    def process_cairo(self, cxt):
        cxt.stroke()


class Fill(Deco):
    def process_cairo(self, cxt):
        cxt.fill()


class FillPreserve(Deco):
    def process_cairo(self, cxt):
        cxt.fill_preserve()


class Clip(Deco):
    def process_cairo(self, cxt):
        cxt.clip()


class RGBA(Deco):
    def __init__(self, r, g, b, a=1.0):
        self.cl = (r, g, b, a)

    def process_cairo(self, cxt):
        cxt.set_source_rgba(*self.cl)


_defaultlinewidth = 0.02 # cm


class LineWidth_Pt(Deco):
    def __init__(self, lw):
        self.lw = lw

    def process_cairo(self, cxt):
        cxt.set_line_width(self.lw)


class LineWidth(LineWidth_Pt):
    def __init__(self, lw):
        self.lw = lw*SCALE_CM_TO_POINT


# cairo constants:
#    cairo.LINE_CAP_BUTT           
#    cairo.LINE_CAP_ROUND         
#    cairo.LINE_CAP_SQUARE       
#    cairo.LINE_JOIN_BEVEL      
#    cairo.LINE_JOIN_MITER     
#    cairo.LINE_JOIN_ROUND    


class LineCap(Deco):
    def __init__(self, desc):
        assert type(desc) is str
        assert desc.lower() in "butt round square".split()
        self.desc = desc.upper()

    def process_cairo(self, cxt):
        import cairo
        cap = getattr(cairo, "LINE_CAP_%s"%(self.desc,))
        cxt.set_line_cap(cap)


class LineJoin(Deco):
    def __init__(self, desc):
        assert type(desc) is str
        assert desc.lower() in "bevel miter round".split()
        self.desc = desc.upper()

    def process_cairo(self, cxt):
        import cairo
        cap = getattr(cairo, "LINE_JOIN_%s"%(self.desc,))
        cxt.set_line_join(cap)


class LineDash_Pt(Deco):
    def __init__(self, dashes, offset=0):
        self.dashes = dashes
        self.offset = offset

    def process_cairo(self, cxt):
        import cairo
        cxt.set_dash(self.dashes, self.offset)


#class LineDash(LineDash_Pt):
#    def __init__(self, dashes, offset=0):
#        dashes = [sz*SCALE_CM_TO_POINT for sz in dashes]
#        LineDash_Pt.__init__(self, dashes, offset)


class LineDash(Deco):
    def __init__(self, dashes, offset=0):
        self.dashes = dashes
        self.offset = offset

    def process_cairo(self, cxt):
        import cairo
        lw = cxt.get_line_width()
        scale = lw / _defaultlinewidth / SCALE_CM_TO_POINT
        dashes = [sz * scale for sz in self.dashes]
        cxt.set_dash(dashes, self.offset)


class TextSize(Deco):
    def __init__(self, size):
        #self.size = float(size) * SCALE_CM_TO_POINT
        self.size = size # ???

    def process_cairo(self, cxt):
        pass
        #cxt.set_font_size(self.size)


class TextAlign(Deco):
    def __init__(self, desc):
        self.desc = desc



class Translate_Pt(Deco):
    def __init__(self, dx, dy):
        self.dx = dx
        self.dy = dy

    def process_cairo(self, cxt):
        cxt.translate(self.dx, -self.dy)


class Translate(Translate_Pt):
    def __init__(self, dx, dy):
        self.dx = dx*SCALE_CM_TO_POINT
        self.dy = dy*SCALE_CM_TO_POINT


class Scale(Deco):
    def __init__(self, sx, sy=None):
        if sy is None:
            sy = sx
        self.sx = float(sx)
        self.sy = float(sy)

    def process_cairo(self, cxt):
        cxt.scale(self.sx, self.sy)


class Rotate(Deco):
    def __init__(self, angle):
        self.angle = angle

    def process_cairo(self, cxt):
        cxt.rotate(self.angle)




# ----------------------------------------------------------------------------
# 
#


def text_extents_cairo(text):
    import cairo
    #surface = cairo.PDFSurface("/dev/null", 0, 0)
    # only in cairo 1.11.0:
    surface = cairo.RecordingSurface(cairo.Content.COLOR_ALPHA, None)
    cxt = cairo.Context(surface)
    #cxt.move_to(0., 0.)
    #cxt.show_text(text)
    #extents = surface.ink_extents()
    #(ulx, uly, width, height) = extents
    #print("get_bound", extents) # same as text_extents ...
    ex = cxt.text_extents(text)
    return ex


class Text(object):
    def __new__(cls, x, y, text):
        ob = object.__new__(the_text_cls)
        return ob


class CairoText(Item, Text):
    def __init__(self, x, y, text):
        self.x = SCALE_CM_TO_POINT*x
        self.y = SCALE_CM_TO_POINT*y
        self.text = text

    def get_bound(self):
        extents = text_extents_cairo(self.text)
        (dx, dy, width, height, x1, y1) = extents
        #print("CairoText.get_bound", repr(self.text), dx, dy, width, height, x1, y1)
        #assert dx>=0, dx
        x, y = self.x, self.y
        llx, lly = x, y-dy-height
        urx, ury = llx+dx+width, lly+height
        b = Bound(llx, lly, urx, ury)
        #print("CairoText.get_bound", b)
        return b

    def process_cairo(self, cxt):
        cxt.save()
        cxt.move_to(self.x, -self.y)
        cxt.show_text(self.text)
        cxt.restore()


class MkText(Compound, Text):

    tex_engine = "pdftex"

    _baseline = None
    @classmethod
    def _get_baseline(cls):
        if cls._baseline is None:
            item = make_text("x") # measure this font using "x"
            bound = item.get_bound()
            #print("_get_baseline", bound)
            cls._baseline = bound.lly
        return cls._baseline

    def __init__(self, x, y, text):
        assert text
        self.x = x = SCALE_CM_TO_POINT*x
        self.y = y = SCALE_CM_TO_POINT*y
        self.text = text
        item = make_text(text, self.tex_engine)
        bound = item.get_bound()
        assert not bound.is_empty(), bound
        llx, lly = bound.llx, bound.lly
        #print("Text.__init__", x, y, bound)
        y0 = self._get_baseline()
        self.bot = y0 - bound.lly
        trafo = Translate_Pt(x-llx, y-lly-self.bot)
        #print("Text.__init__ trafo:", trafo)
        self.bound = Bound(x, y-self.bot, x+bound.width, y+bound.height-self.bot)
        items = list(item.items)
        items = [trafo] + items
        Compound.__init__(self, items)

    def get_bound(self):
        return self.bound


the_text_cls = CairoText

#def Text(*args, **kw):
#    return the_text_cls(*args, **kw)



# ----------------------------------------------------------------------------
# Some good code copied from PyX
#


def arc_to_curve_pt(x_pt, y_pt, r_pt, angle1, angle2):
    dangle = angle2-angle1

    if dangle==0:
        return None

    x0_pt, y0_pt = x_pt+r_pt*cos(angle1), y_pt+r_pt*sin(angle1)
    x3_pt, y3_pt = x_pt+r_pt*cos(angle2), y_pt+r_pt*sin(angle2)

    l = r_pt*4*(1-cos(dangle/2))/(3*sin(dangle/2))

    x1_pt, y1_pt = x0_pt-l*sin(angle1), y0_pt+l*cos(angle1)
    x2_pt, y2_pt = x3_pt+l*sin(angle2), y3_pt-l*cos(angle2)

    items = [
        LineTo_Pt(x0_pt, y0_pt), 
        CurveTo_Pt(x1_pt, y1_pt, x2_pt, y2_pt, x3_pt, y3_pt)]
    return items


def arc_to_bezier(x, y, r, angle1, angle2, danglemax=0.5*pi):
    x_pt = x*SCALE_CM_TO_POINT
    y_pt = y*SCALE_CM_TO_POINT
    r_pt = r*SCALE_CM_TO_POINT
    p = arc_to_bezier_pt(x_pt, y_pt, r_pt, angle1, angle2, danglemax)
    return p


def arc_to_bezier_pt(x_pt, y_pt, r_pt, angle1, angle2, danglemax=0.5*pi):
    if angle2<angle1:
        angle2 = angle2 + (floor((angle1-angle2)/(2*pi))+1)*2*pi
    elif angle2>angle1+2*pi:
        angle2 = angle2 - (floor((angle2-angle1)/(2*pi))-1)*2*pi

    if r_pt == 0 or angle1-angle2 == 0:
        return []

    subdivisions = int((angle2-angle1)/danglemax)+1

    dangle = (angle2-angle1)/subdivisions

    items = []
    for i in range(subdivisions):
        items += arc_to_curve_pt(x_pt, y_pt, r_pt, angle1+i*dangle, angle1+(i+1)*dangle)

    p = Path(items)
    return p


# ----------------------------------------------------------------------------
#
#


def test():
    pass


if __name__ == "__main__":
    test()





