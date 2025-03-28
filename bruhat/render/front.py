#!/usr/bin/env python3

"""
Expose a _drawing api similar to the pyx api.

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

from bruhat.render.base import EPSILON, NS, SCALE_CM_TO_POINT
from bruhat.render.back import * # XXX
from bruhat.render.back import _defaultlinewidth
from bruhat.render.flatten import Flatten


# ----------------------------------------------------------------------------
# 
#



path = NS(
    line=Line, curve=Curve, rect=Rect, circle=Circle, path=Path,
    arc=Arc, arcn=Arcn, moveto=MoveTo, lineto=LineTo, 
    curveto=CurveTo, closepath=ClosePath)



RGB = RGBA
RGB.red = RGB(1., 0., 0.)
RGB.green = RGB(0., 1., 0.)
RGB.blue = RGB(0., 0., 1.)
RGB.white = RGB(1., 1., 1.)
RGB.grey = RGB.gray = RGB(0.8, 0.8, 0.8)
RGB.black = RGB(0., 0., 0.)

color = NS(rgb=RGBA)


style = NS(
    linecap = NS(
        butt = LineCap("butt"),
        round = LineCap("round"), 
        square = LineCap("square")),
    linejoin = NS(
        bevel = LineJoin("bevel"),
        miter = LineJoin("miter"),
        round = LineJoin("round")),
    linewidth = NS(
        THIN = LineWidth(_defaultlinewidth/sqrt(32)),
        THIn = LineWidth(_defaultlinewidth/sqrt(16)),
        THin = LineWidth(_defaultlinewidth/sqrt(8)),
        Thin = LineWidth(_defaultlinewidth/sqrt(4)),
        thin = LineWidth(_defaultlinewidth/sqrt(2)),
        normal = LineWidth(_defaultlinewidth),
        thick = LineWidth(_defaultlinewidth*sqrt(2)),
        Thick = LineWidth(_defaultlinewidth*sqrt(4)),
        THick = LineWidth(_defaultlinewidth*sqrt(8)),
        THIck = LineWidth(_defaultlinewidth*sqrt(16)),
        THICk = LineWidth(_defaultlinewidth*sqrt(32)),
        THICK = LineWidth(_defaultlinewidth*sqrt(64)),
    ))


# XXX TODO
text = NS(
    size = NS(
        tiny = TextSize(-4),
        script = TextSize(-3),
        footnote = TextSize(-2),
        small = TextSize(-1),
        normal = TextSize(0),
        large = TextSize(1),
        Large = TextSize(2),
        LARGE = TextSize(3),
        huge = TextSize(4),
        Huge = TextSize(5)),
    halign = NS(
        left = TextAlign("left"),
        center = TextAlign("center"),
        right = TextAlign("right"),
        clear = TextAlign("clear"),
        boxleft = TextAlign("boxleft"),
        boxcenter = TextAlign("boxcenter"),
        boxright = TextAlign("boxright"),
        flushleft = TextAlign("flushleft"),
        flushcenter = TextAlign("flushcenter"),
        flushright = TextAlign("flushright")),
    valign = NS(
        top = TextAlign("top"),
        middle = TextAlign("middle"),
        bottom = TextAlign("bottom")))




linestyle = NS(
    solid = CompoundDeco([style.linecap.butt, LineDash([])]),
    dashed = CompoundDeco([style.linecap.butt, LineDash([2])]),
    dotted = CompoundDeco([style.linecap.round, LineDash([0, 2])]),
    dashdotted = CompoundDeco([style.linecap.round, LineDash([0, 2, 2, 2])]))


trafo = NS(translate = Translate, scale = Scale, rotate = Rotate)



class Canvas(Compound):

    def stroke(self, path, decos=[]):
        assert type(decos) is list
        assert isinstance(path, Item), repr(path)
        item = Compound(decos, path, Stroke())
        self.append(item)

    def fill(self, path, decos=[]):
        #assert type(decos) is list
        assert isinstance(path, Item), repr(path)
        item = Compound(decos, path, Fill())
        self.append(item)

    def clip(self, path):
        self.append(path)
        self.append(Clip())

#    def text_extents(self, text):
#        dx, dy, width, height, _, _ = text_extents_cairo(text)
#        return (dx/SCALE_CM_TO_POINT, -dy/SCALE_CM_TO_POINT,  # <-- sign flip
#            width/SCALE_CM_TO_POINT, height/SCALE_CM_TO_POINT)

    def text_extents(self, text):
        item = Text(0., 0., text)
        bound = item.get_bound()
        llx, lly, urx, ury = bound
        llx /= SCALE_CM_TO_POINT
        lly /= SCALE_CM_TO_POINT
        urx /= SCALE_CM_TO_POINT
        ury /= SCALE_CM_TO_POINT
        return (0., ury, urx-llx, ury-lly)

    def text(self, x, y, text, decos=[]):
        assert type(decos) is list
        color = None
        for deco in decos:
            if isinstance(deco, RGBA):
                color = deco
        item = Compound(decos, Text(x, y, text, color))
        self.append(item)

    def _write_cairo(self, method, name):

        if 0:
            #self.dump()
            #bound = self.get_bound()
            #print("_write_cairo: self.get_bound()", bound)
    
            cxt = Flatten()
            self.process_cairo(cxt)
            item = Compound(cxt.paths)
            #print("Flatten:")
            #item.dump()
            bound = item.get_bound()
            #print("_write_cairo: item.get_bound()", bound)
            assert not bound.is_empty()

        else:
            bound = self.get_bound_cairo()

        import cairo

        W = bound.width
        H = bound.height
        surface = method(name, W, H)

        dx = 0 - bound.llx
        dy = H + bound.lly
        surface.set_device_offset(dx, dy)

        cxt = cairo.Context(surface)
        cxt.set_line_width(_defaultlinewidth * SCALE_CM_TO_POINT)
        self.process_cairo(cxt)
        #item.process_cairo(cxt)
        return surface

    def writePDFfile(self, name):
        assert name.endswith(".pdf") or name=="/dev/null" 
        import cairo
        method = cairo.PDFSurface
        surface = self._write_cairo(method, name)
        surface.finish()

    def writeSVGfile(self, name):
        assert name.endswith(".svg") or name=="/dev/null" 
        import cairo
        method = cairo.SVGSurface
        surface = self._write_cairo(method, name)
        surface.finish()

    def writePNGfile(self, name):
        assert name.endswith(".png") or name=="/dev/null" 
        import cairo
        def method(name, W, H):
            W = int(round(W))
            H = int(round(H))
            surface = cairo.ImageSurface(cairo.Format.RGB24, W, H)
            return surface
        surface = self._write_cairo(method, name)
        surface.write_to_png(name)
        surface.finish()


canvas = NS(canvas=Canvas)



# ----------------------------------------------------------------------------
#
#


def test():

    cvs = canvas.canvas()

    def cross(x, y):
        r = 0.1
        st = [color.rgb.blue, style.linewidth.THick, style.linecap.round]
        cvs.stroke(path.line(x-r, y-r, x+r, y+r), st)
        cvs.stroke(path.line(x-r, y+r, x+r, y-r), st)

    p = path.path([
        path.moveto(0., 0.),
        path.arc(0., 0., 1., 0., 0.5*pi),
        path.lineto(-1., 1.),
        path.arc(-1., 0., 1., 0.5*pi, 1.0*pi),
        path.arc(-1.5, 0., 0.5, 1.0*pi, 2.0*pi),
        path.closepath()
    ])

    items = (
    [ 
        path.moveto(0., 0.),
        path.arc(0., 0., 1., 0., 0.5*pi),
        path.lineto(-1., 1.), path.arc(-1., 0., 1., 0.5*pi, 1.0*pi),
        path.arc(-1.5, 0., 0.5, 1.0*pi, 2.0*pi), path.closepath() ])
    p = path.path(items)

    cvs.fill(p, [color.rgb.red, trafo.scale(0.8, 0.8)])
    cvs.stroke(p, [color.rgb.black, style.linewidth.THick])

    cross(0., 0.)
    cross(-1.2, 1.2)

    if 0:
        x, y, r, angle1, angle2 = 0., 0., 1., 0., 0.5*pi
        p = arc_to_bezier(x, y, r, angle1, angle2, danglemax=pi/2.)
        cvs.stroke(p, [color.rgb.white])
    
        x, y, r, angle1, angle2 = 0., 0., 1., -0.5*pi, 0.
        p = arc_to_bezier(x, y, r, angle1, angle2, danglemax=pi/2.)
        cvs.stroke(p, [color.rgb.red])

    cvs.writePDFfile("output.pdf")

    print("OK")


def test():

    cvs = canvas.canvas()

    def cross(x, y):
        r = 0.1
        st = [color.rgb.blue, style.linewidth.normal, style.linecap.round]
        cvs.stroke(path.line(x-r, y-r, x+r, y+r), st)
        cvs.stroke(path.line(x-r, y+r, x+r, y-r), st)

    #cvs.append(Translate(1., 1.))
    cross(0., 0.)

    cvs.text(0., 0., "hey there!")

    cvs.writePDFfile("output.pdf")

    print("OK\n")


if __name__ == "__main__":
    test()





