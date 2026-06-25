#!/usr/bin/env python3

"""

Render circles and lines under Mobius group actions.

See:
"An Introduction to the McKay Correspondence
Master Thesis in Physics"
Max Lindh


"""

import sys
from random import randint, seed, choice, random, shuffle
from functools import reduce, cache
from operator import mul
from cmath import pi, exp
from math import atan

import numpy

from huygens.namespace import (
    Canvas, Scale, Rotate, MoveTo, LineTo,
    red, green, blue, grey, black, white, orange,
    path, st_THick, st_center, thin)
from huygens.back import arc_to_bezier, _arc_to_bezier

if __name__=="__main__":
    from huygens import config
    config(text="pdflatex", latex_header=r"\usepackage{amsmath}\usepackage{amssymb}")


def save(cvs, name):
    print("save", name)
    assert "images" not in name
    name = "images/"+name
    cvs.writePDFfile(name)

from bruhat.util import cross
from bruhat.smap import SMap
from bruhat.todd_coxeter import Schreier
from bruhat.mobius import Mobius, mulclose, EPSILON, d_poincare
from bruhat.disc import get_pos_angle, get_angle
from bruhat.argv import argv



rnd = lambda r=1: (1-2*random())*r + (1-2*random())*r*1j
metric = lambda w,z=0. : ((w-z) * complex(w-z).conjugate()).real**0.5
eq = lambda z0, z1: abs(complex(z0-z1)) < EPSILON

infty = None
I = Mobius()

rotate = lambda theta : Mobius.rotate(theta/2) # double cover .... ?!?
translate = lambda x : Mobius(1, x, 0, 1)

#_det = lambda a,b,c,d:a*d-b*c
#det = lambda z,w:_det(z.real, z.imag, w.real, w.imag)


def zpromote(z):
    return complex(z) if z is not None else z

def zeq(z0, z1):
    if z0==z1:
        return True
    if z0 is None or z1 is None:
        return False
    return abs(z0-z1)<EPSILON


def get_BT(refl=False):
    # Binary tetrahedral group (see Lindh p32)
    r2 = 2**0.5
    r3 = 3**0.5
    i = 1j
    a = Mobius((1+i*r3)/2, 0, 0, (1-i*r3)/2)
    assert a.order() == 6
    s = 2*r3
    b = Mobius((r3-i)/s, (-i*2*r2)/s, (-i*2*r2)/s, (r3+i)/s)
    assert b.order() == 6

    gen = [a,b]
    if refl:
        c = Mobius(1, 0, 0, 1, True) # conjugate
        gen = [a, b, c]

    G = mulclose(gen)
    assert len(G) == 48 if refl else 24
    #print("get_BT:", len(G))
    return G


def get_BO(refl=False):
    # Binary octahedral group (Lindh p34)
    r2 = 2**0.5
    i = 1j
    a = Mobius((1+i)/r2, 0, 0, (1-i)/r2)
    b = Mobius(1/r2, -i/r2, -i/r2, 1/r2)
    gen = [a,b]
    if refl:
        c = Mobius(1, 0, 0, 1, True) # conjugate
        gen = [a, b, c]
    G = mulclose(gen)
    assert len(G) == 96 if refl else 48 
    #print("get_BO:", len(G))
    return G


def get_BD(refl=False):
    # Binary dodecahedral group
    i = 1j
    e = exp(i*pi/5)
    r5 = 5**0.5

    a = Mobius(e, 0, 0, e**9)
    b = Mobius(
        (5*(e-e**4)-r5*(e+e**4))/10,
        -2*r5*(e+e**4)/10, 
        -2*r5*(e+e**4)/10, 
        (5*(e-e**4)+r5*(e+e**4))/10,
    )

    gen = [a,b]
    if refl:
        c = Mobius(1, 0, 0, 1, True) # conjugate
        gen = [a, b, c]
    G = mulclose(gen)
    assert len(G) == 240 if refl else 120
    #print("get_BD:", len(G))
    return G


def get_orbit(G, item):
    found = []
    for g in G:
        gtem = g*item
        if gtem not in found:
            found.append(gtem)
    return found

def unique(items):
    found = []
    for item in items:
        if item not in found:
            found.append(item)
    return found


stack = []
def push(): # campaign to save globals from extinction
    #print("push")
    save = {}
    for name,cls in globals().items():
        if (isinstance(cls, object.__class__) 
            and issubclass(cls, (Feature,Geometry))):
            ns = cls.__dict__
            for (attr, value) in ns.items():
                #print(k, list(ns.keys()))
                if attr.startswith("_"):
                    continue
                if isinstance(value, (list,tuple,type(None),int,float,complex)):
                    #print("\t", name, attr, value)
                    save[(cls, attr)] = value
    stack.append(save)

def pop():
    #print("pop")
    save = stack.pop()
    for key,value in save.items():
        cls, attr = key
        #print("\t", cls.__name__, attr, value)
        setattr(cls, attr, value)

    
# ----------------------------------------------------------------------------
#
#     Feature's live in Geometry's
#
# ----------------------------------------------------------------------------

class Feature:
    scale = 1.0

    def __init__(self, zs, **kw):
        if zs is None or isinstance(zs, (float, complex, int)):
            zs = [zs]
        zs = [zpromote(z) for z in zs]
        self.zs = zs
        self.kw = kw
        self.__dict__.update(kw)

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.zs,)
    __repr__ = __str__

    def __rmul__(self, g):
        assert isinstance(g, Mobius)
        zs = [g(z) for z in self.zs]
        return self.__class__(zs, **self.kw)

    def __eq__(self, other):
        for z,w in zip(self.zs, other.zs):
            if not zeq(z,w):
                return False
        return True

    def render(self, cvs):
        pass


class Point(Feature):

    scale = 1.0
    radius = 0.06
    st = [black]

    def __eq__(self, other):
        assert isinstance(other, Point)
        z0, z1 = self.zs[0], other.zs[0]
        return zeq(z0, z1)
        
    def render(self, cvs):
        st = self.st
        radius = self.radius
        z = self.zs[0]
        if z is None:
            cvs.stroke(path.circle(0, 0, Spherical.RADIUS), st+st_THick)
        else:
            cvs.fill(path.circle(z.real, z.imag, self.scale*radius), st)

    def line(self, other):
        z0 = self.zs[0]
        z2 = other.zs[0]
        z1 = (z0+z2)/2
        return Arc([z0, z1, z2])


class Circle(Feature):
    """ A circle (or line) in the complex plane 
        defined by (passing through) three given points.
    """ 
    def __init__(self, zs, **kw):
        assert len(zs) == 3, zs
        assert zs.count(None) <= 1, zs
        Feature.__init__(self, zs, **kw)
        a, b, c = self.zs
        assert not zeq(a, b), zs
        assert not zeq(a, c), zs
        assert not zeq(b, c), zs

    def __eq__(self, other):
        assert isinstance(other, Circle)
        l, r = self.is_colinear(), other.is_colinear()
        if l and r:
            a, b, c = self.zs
            d, e, f = other.zs
            if a is None:
                a, b, c = b, c, a # rotate points
            if d is None:
                d, e, f = e, f, d # rotate points
            ab = a-b
            de = d-e
            if abs((ab/de).imag) > EPSILON:
                return False
            # parallel lines
            if abs(a-e) < EPSILON:
                return True
            ae = a-e
            return abs((ab/ae).imag) < EPSILON
        elif not l and not r:
            z0 = self.get_center()
            z1 = other.get_center()
            if abs(z0-z1) > EPSILON:
                return False
            r0 = self.get_radius()
            r1 = other.get_radius()
            return abs(r0-r1) < EPSILON
        else:
            return False

    def is_colinear(self):
        if None in self.zs:
            return True
        a, b, c = self.zs
        a, b = a-c, b-c
        ab = a/b
        return abs(ab.imag) < EPSILON

    def get_center(self): # XXX cache
        assert not self.is_colinear()
        a, b, c = self.zs
        ax, ay = a.real, a.imag
        bx, by = b.real, b.imag
        cx, cy = c.real, c.imag
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        #assert abs(d)>EPSILON, "is_colinear"
        assert abs(d)>1e-14, "yikes"
    
        ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
        uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    
        return ux + 1j*uy

    def get_radius(self):
        assert not self.is_colinear()
        z0 = self.get_center()
        return metric(z0, self.zs[0])

    st_fill = None
    st_stroke = []

    #def render(self, cvs, st_stroke=[], st_fill=None):
    def render(self, cvs):
        if self.is_colinear():
            a, b, c = self.zs
            dz = b-a
            dz = dz/metric(dz)
            #z0, z1 = a-100*dz, a+100*dz
            # XXX check we go through zero XXX
            z0, z1 = -Spherical.RADIUS*dz, Spherical.RADIUS*dz
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)
        
        else:
            z0 = self.get_center()
            r = self.get_radius()
            p = path.circle(z0.real, z0.imag, r)
        if self.st_fill is not None:
            cvs.fill(p, self.st_fill)
        if self.st_stroke is not None:
            cvs.stroke(p, self.st_stroke)


class Arc(Circle):
    "directed arc of a circle defined by three points: start, (random) middle, end"
    def __eq__(self, other):
        assert isinstance(other, Arc)
        if not Circle.__eq__(self, other):
            return False
        # FIX FIX TODO TODO FIX FIX
        return (zeq(self.zs[0], other.zs[0]) and zeq(self.zs[2], other.zs[2])
        or zeq(self.zs[0], other.zs[2]) and zeq(self.zs[2], other.zs[0])) # no ...

    @property
    def src(self):
        return self.zs[0]

    @property
    def tgt(self):
        return self.zs[2]

    def get_path(self):
        zs = self.zs
        if None in zs:
            # arc to infty
            assert zs.count(None)==1, zs
            assert zs[1] is not None, zs # ouch: "middle" should not be infty !?!
            z0, z1, z2 = zs
            z1 = z1*Spherical.RADIUS/abs(z1)
            if z0 is None:
                p = path.line(z1.real, z1.imag, z2.real, z2.imag)
            elif z2 is None:
                p = path.line(z0.real, z0.imag, z1.real, z1.imag)

        elif self.is_colinear():
            z0, _, z1 = zs
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)

        else:
            zc = self.get_center()
            r = self.get_radius()
            assert r>EPSILON
            xc, yc = zc.real, zc.imag
            z0, z1, z2 = zs
            dz0, dz1, dz2 = [z-zc for z in self.zs]
            dz0, dz1, dz2 = [1, dz1/dz0, dz2/dz0]
            swap = get_pos_angle(dz1) > get_pos_angle(dz2)
            if swap:
                z0, z2 = z2, z0
            theta0 = get_pos_angle(z0-zc)
            theta1 = get_pos_angle(z1-zc)
            theta2 = get_pos_angle(z2-zc)
            #print(theta0, theta1, theta2)
            p = arc_to_bezier(xc, yc, r, theta0, theta2, relative=False)
            #items = _arc_to_bezier(xc, yc, r, theta0, theta2, relative=False)
            #p = path.path(items)
            # HACK THIS ?!?
            #items = [item for item in p.items if not isinstance(item, LineTo)]
            #p = path.path(items)
            if swap:
                p = p.backwards()
                #print("backwards")

        return p
        
    st_fill = None
    st_stroke = []

    def render(self, cvs):
        p = self.get_path()
        st = [self.scale * thin] + self.st_stroke
        cvs.stroke(p, st)


class Bigon(Arc):

    line_width = 0.01
    st = []
    def render(self, cvs):

        assert None not in self.zs, "not implemented"
        a, b, c = self.zs
        if self.is_colinear():
            dz = (c-a)*1j

        else:
            z0 = self.get_center()
            r = self.get_radius()
            dz = b-z0

        dz = dz / abs(dz)
        line_width = self.line_width * self.scale
        b0 = b - line_width*dz
        b1 = b + line_width*dz
        larc = Arc([a, b0, c])
        rarc = Arc([c, b1, a])
        p = larc.get_path() + rarc.get_path()
        cvs.fill(p, self.st)



class Triangle(Feature):
    def __init__(self, zs, **kw):
        Feature.__init__(self, zs, **kw)
        assert len(zs) == 6
        p0, m0, p1, m1, p2, m2 = zs
        arcs = [
            Arc([p0, m0, p1]),
            Arc([p1, m1, p2]),
            Arc([p2, m2, p0])]
        self.arcs = arcs

    def get_center(self):
        assert None not in self.zs, "TODO"
        p0, m0, p1, m1, p2, m2 = self.zs
        z = (p0+p1+p2)/3
        return z

    def is_clockwise(self): # or is it anti-clockwise ?
        assert None not in self.zs, "TODO"
        zc = self.get_center()
        z0, z1 = self.zs[:2]
        return get_angle((z1-zc)/(z0-zc)) > 0.

    #def render(self, cvs, st_stroke=None, st_fill=None, lbl=None, err=0.00, debug=False):
    st_stroke = None
    st_fill = [orange.alpha(0.5)]
    lbl = None
    err = 0.00
    def render(self, cvs):
        zs = self.zs
        arcs = list(self.arcs)
        if None in zs:
            # this is an unbounded triangle, add an extra Arc at the RADIUS
            assert zs.count(None)==1
            assert zs[4] is None, "i'm expecting the green here"
            a, b, c = arcs
            x0, y0 = b.get_path().getat(1)
            x2, y2 = c.get_path().getat(0)
            r = (x0**2 + y0**2)**0.5
            x1, y1 = (x0+x2)*0.5, (y0+y2)*0.5
            r1 = (x1**2 + y1**2)**0.5
            x1 = r*x1/r1
            y1 = r*y1/r1
            arc = Arc([x0+y0*1j, x1+y1*1j, x2+y2*1j])
            arcs.insert(2, arc)

        ps = [arc.get_path() for arc in arcs]
        items = list(ps[0].items)
        for p in ps[1:]:
            jtems = list(p.items)
            if isinstance(jtems[0], MoveTo):
                jtems.pop(0)
            items += jtems
        p = path.path(items)
        if self.lbl is not None:
            zc = self.get_center()
            cvs.text(zc.real+self.err*random(), zc.imag+self.err*random(), 
                self.lbl, [Scale(0.5)]+st_center)

        if self.st_fill is not None:
            cvs.fill(p, self.st_fill)
        if self.st_stroke is not None:
            cvs.stroke(p, self.st_stroke)


# ----------------------------------------------------------------------------
#
#     Geometry's
#
# ----------------------------------------------------------------------------



class Geometry:
    scale = 3.0

    def __init__(self, G, features, **kw):
        self.G = G
        self.features = list(features)
        self.__dict__.update(kw)
        self.kw = kw

    @classmethod
    def from_arcs(cls, G, rb_arc, bg_arc, gr_arc):
        z_red = rb_arc.src
        z_blue = rb_arc.tgt

        # make a Triangle:
        for g in G:
            item = g*bg_arc
            if zeq(item.src, z_blue):
                bg_arc = item
                break
        else:
            assert 0
        z_green = bg_arc.tgt

        for g in G:
            item = g*gr_arc
            if zeq(item.src, z_green) and zeq(item.tgt, z_red):
                gr_arc = item
                break
        else:
            assert 0
        assert zeq(gr_arc.tgt, z_red), (gr_arc, z_red)
        p_red = Point([z_red], st=[red])
        p_blue = Point([z_blue], st=[blue])
        p_green = Point([z_green], st=[green])
        triangle = Triangle([
            rb_arc.zs[0], rb_arc.zs[1], 
            bg_arc.zs[0], bg_arc.zs[1], 
            gr_arc.zs[0], gr_arc.zs[1],
        ], st_fill=[orange.alpha(0.3)])

        features = [triangle, rb_arc, bg_arc, gr_arc, p_red, p_blue, p_green]
        geometry = cls(G, features)
        return geometry

    def get_oriented(self):
        G = [g for g in self.G if not g.conjugate]
        return self.__class__(G, self.features, **self.kw)

    def get_orbit(self, item):
        G = self.G
        found = []
        for g in G:
            gtem = g*item
            if gtem not in found:
                found.append(gtem)
        return found

    def get_cvs(self):
        return Canvas([Scale(self.scale)])

    def render(self, cvs=None):
        G = self.G
        cvs = cvs if cvs is not None else self.get_cvs()
        for item in self.features:
            #found = unique([g*item for g in G])
            found = self.get_orbit(item)
            for item in found:
                item.render(cvs)
        return cvs
    
    @staticmethod
    def get_BT():
        z_red = 0.0
        z_blue = (2-3**0.5)**0.5
        z_green = (-1/2**0.5)
    
        rb_arc = Arc([z_red, z_blue/2, z_blue])
        gr_arc = Arc([z_green, z_green/2, z_red])
        z_blue = z_green - (z_blue-z_green)
        bg_arc = Arc([z_blue, 2*z_blue, infty])
    
        G = get_BT(True)
    
        geometry = Spherical.from_arcs(G, rb_arc, bg_arc, gr_arc)
        return geometry
    
    @staticmethod
    def get_BO():
        G = get_BO(True)
    
        z_green = 0.0
        z_red = 0.3660254037844387+0.3660254037844386j
        z_blue = 2.4142135623730
        bg_arc = Arc([z_blue, 2*z_blue, infty])
        gr_arc = Arc([z_green, (z_red+z_green)/2, z_red])
    
        for g in G:
            p = g*Point([z_blue])
            z = p.zs[0] # arghhh...
            if z.real > EPSILON and z.imag > EPSILON:
                z_blue = z
                break
        else:
            assert 0
    
        rb_arc = Arc([z_red, (z_red+z_blue)/2, z_blue])
    
        geometry = Spherical.from_arcs(G, rb_arc, bg_arc, gr_arc)
        return geometry
    
    @staticmethod
    def get_BD():
        G = get_BD(True)
    
        #z_green = 0.0
        z_red = -2.956295201467611
        z_blue = 3.52014702134020
        bg_arc = Arc([z_blue, 2*z_blue, infty])
        gr_arc = Arc([infty, 2*z_red, z_red])
    
        for g in G:
            p = g*Point([z_blue])
            z = p.zs[0] # arghhh...
            if abs(z.imag) < EPSILON and z.real < -1.:
                z_blue = z
                break
        else:
            assert 0
    
        rb_arc = Arc([z_red, (z_red+z_blue)/2, z_blue])
    
        geometry = Spherical.from_arcs(G, rb_arc, bg_arc, gr_arc)
        return geometry


class Spherical(Geometry):
    scale = 3.0
    RADIUS = 4.0 
    def get_cvs(self):

        R = Spherical.RADIUS
        cvs = Canvas().scale(self.scale)
        p = path.circle(0, 0, R)
        #cvs.stroke(p, [green])
        cvs.stroke(path.circle(0,0,1.1*R), [white])
        cvs.clip(p)
        return cvs


class Affine(Geometry):
    R = 1.7
    scale = 10.0
    radius = 0.04

    @classmethod
    def clip(cls, items):
        R = cls.R
        jtems = []
        for item in items:
            for z in item.zs:
                if z.real > -R and z.real < +R and z.imag > -R and z.imag < +R:
                    jtems.append(item)
                    break
        return jtems

    @classmethod
    def accept(cls, g):
        R = cls.R
        z = g(0) * 0.7
        return z.real > -R and z.real < +R and z.imag > -R and z.imag < +R

    @classmethod
    def generate(cls, gens, maxsize=2000):
        G = mulclose(gens, maxsize=maxsize, accept=cls.accept)
        return G

    @classmethod
    def get_cvs(cls):
        R = cls.R
        scale = cls.scale
        cvs = Canvas([Scale(scale)])
        RR = 1.1*R
        p = path.rect(-RR, -RR, 2*RR, 2*RR)
        cvs.stroke(p, [white])
        p = path.rect(-R, -R, 2*R, 2*R)
        cvs.stroke(p, [grey]+st_THick)
        cvs.fill(p, [white])
        cvs.clip(p)
        return cvs



class Hyperbolic(Geometry):

    def get_orbit(self, item):
        orbit = []
        for item in Geometry.get_orbit(self, item):
            zs = [z for z in item.zs if z is not None]
            z = sum(zs)/len(zs)
            item.scale = 1/d_poincare(z)
            orbit.append(item)
        return orbit
        
        


# ----------------------------------------------------------------------------
#
#     drivers
#
# ----------------------------------------------------------------------------


def test():

    a = Mobius.rotate(2*pi/6)
    assert a.order() == 6

    m = Mobius.send(1, 1j, -1, 1, -1, 1j)
    assert eq(m(1), 1)
    assert eq(m(1j), -1)
    assert eq(m(-1), 1j)

    c = Circle([1., 1j, -1])
    assert eq(c.get_center(), 0)
    assert eq(c.get_radius(), 1)
    d = a*c

    assert c==c
    assert c==a*c
    assert c != Circle([1.1, 1j, -1])

    assert not c.is_colinear()
    z = 1+2.3j
    assert Circle([0+z, 1+z, 2+z]).is_colinear()

    for _ in range(10):
        c = Circle([rnd() for i in range(3)])

        z0 = c.get_center()

        r = None
        for z in c.zs:
            r1 = metric(z, z0)
            assert r is None or abs(r-r1) < EPSILON


def test_circles():

    scale = 3.0
    radius = 0.06

    cl = 0.6*white
    Circle.st_stroke = [cl]+st_THick

    G = get_BT()

    R = Spherical.RADIUS
    cvs = Canvas().scale(scale)
    p = path.circle(0, 0, R)
    #cvs.stroke(p, [green])
    cvs.stroke(path.circle(0,0,1.1*R), [white])
    #cvs.clip(p)
    
    #c = Circle([rnd() for i in range(3)])
    c = Circle([0, 1, 2])
    orbit = get_orbit(G, c)
    for item in orbit:
        item.render(cvs)

    for c in orbit:
        if c.is_colinear():
            continue
        z0 = c.get_center()
        if abs(z0.imag) < EPSILON:
            r = c.get_radius()
            z_blue = z0+r
    z_green = infty
    z_red = 0.

    for point in [Point(z_green, st=[green]), Point(z_red, st=[red]), Point(z_blue, st=[blue])]:
        orbit = get_orbit(G, point)
        #print("orbit:", len(orbit))
        for point in orbit:
            point.render(cvs)

    #for g in G:
    #    z = g(0.2 + 0.1j)
    #    cvs.text(z.real, z.imag, r"$\star$", st_center)

    # -----------------------------------------------------

    cvs = Canvas([cvs]).show_page().scale(scale)

    G = get_BO()
    p = path.circle(0, 0, R)
    #cvs.stroke(p, [green])
    cvs.stroke(path.circle(0,0,1.1*R), [white])
    #cvs.clip(p)
    
    diag = lambda x:x+x*1j

    c0 = Circle([diag(x) for x in [0.5, 1.1, 2]])
    c1 = Circle([0.5, 1.1, 2])

    orbit = get_orbit(G, c0) + get_orbit(G, c1)
    #print("orbit:", len(orbit))
    cl = 0.6*white
    for item in orbit:
        item.render(cvs)
        if item.is_colinear():
            continue
        z = item.get_center()
        if abs(z.imag) < EPSILON and z.real < -EPSILON:
            r = item.get_radius()
            z_blue = z-r
        #if abs(z.real) < EPSILON and z.imag > EPSILON:
        #    item.render(cvs)

    for g in G:
        pts = g.fixed()
        a, b = pts
        if abs(a.real-a.imag) < EPSILON and EPSILON<a.real<0.5:
            #Point([a]).render(cvs)
            z_red = a
    z_green = 0.

    for point in [Point(z_green, st=[green]), Point(z_red, st=[red]), Point(z_blue, st=[blue])]:
        orbit = get_orbit(G, point)
        for point in orbit:
            point.render(cvs)

    # -----------------------------------------------------
    # -----------------------------------------------------

    cvs = Canvas([cvs]).show_page().scale(scale)

    G = get_BD()
    p = path.circle(0, 0, R)
    #cvs.stroke(p, [green])
    cvs.stroke(path.circle(0,0,1.1*R), [white])
    #cvs.clip(p)
    
    diag = lambda x:x+x*1j

    c = Circle([0.5, 1.1, 2])

    orbit = get_orbit(G, c)
    #print("orbit:", len(orbit))
    cl = 0.6*white
    for item in orbit:
        item.render(cvs)

    items = []
    for g in G:
        #(g*p).render(cvs, radius, [black])
        pts = g.fixed()
        a, b = pts
        items.append(Point([a]))
        items.append(Point([b]))
        #Point([a]).render(cvs, radius, [green])
        #Point([b]).render(cvs, radius, [green])

    z_green = 0.
    found = unique(items)

    orbit = get_orbit(G, Point([0]))
    for item in orbit:
        #item.render(cvs, radius, [green])
        if item in found:
            found.remove(item)

    #print(len(found))
    for item in found:
        z = item.zs[0]
        if z.real>2 and abs(z.imag)<EPSILON:
            #item.render(cvs, radius, [red])
            z_blue = z
        if z.real< -2 and abs(z.imag)<EPSILON:
            #item.render(cvs, radius, [red])
            z_red = z

    for point in [Point(z_green, st=[green]), Point(z_red, st=[red]), Point(z_blue, st=[blue])]:
        orbit = get_orbit(G, point)
        for point in orbit:
            point.render(cvs)

    save(cvs, "stereo_circles.pdf")


def test_stereo():

    cvs = Canvas()
    
    for geometry in [Geometry.get_BT(), Geometry.get_BO(), Geometry.get_BD()]:
        geometry = geometry.get_oriented()
        fg = geometry.render()
        cvs.append(fg)
        cvs.show_page()
    
    save(cvs, "stereo.pdf")


def test_worksheet():

    cvs = Canvas()

    geometry = Geometry.get_BT()
    fg = geometry.render()
    cvs.insert(0, 0, fg)
    #bb = fg.get_bound_box()
    #cvs.insert(0, -bb.height, fg)

    save(cvs, "stereo_worksheet.pdf")



def test_wallpaper():

    push() # because we love globals don't we
    z_red = 0.
    z_blue = 0.5
    z_green = (0 + 1 + exp(2*1j*pi/6))/3.

    r_rot = rotate(2*pi/6)
    b = translate(1)
    g_refl = Mobius.conjugate()
    assert zeq(g_refl(z_blue), z_blue)

    G = Affine.generate([r_rot, b])

#    G = mulclose([r_refl, b_refl, g_refl], maxsize=200)

    Point.radius = 0.04
    p_red = Point([z_red], st=[red])
    p_blue = Point([z_blue], st=[blue])
    p_green = Point([z_green], st=[green])
    
    rb_arc = p_red.line(p_blue)
    gr_arc = p_green.line(p_red)
    bg_arc = p_blue.line(p_green)

    geometry = Affine(G, [rb_arc, bg_arc, gr_arc, p_red, p_blue, p_green])
    cvs = geometry.render()
    save(cvs, "stereo_wallpaper_refl.pdf")
    
    geometry = Affine.from_arcs(G, rb_arc, bg_arc, gr_arc)
    cvs = geometry.render()

    save(cvs, "stereo_wallpaper.pdf")

    pop()


def test_hyperbolic():
    push()
    Point.radius = 0.04
    #Geometry.scale = 4.0

    from bruhat.disc import mktriangle
    (l, m, n, maxsize) = (7,2,3,400)

    # build the rotation group generators
    a, b = [g.todisc() for g in mktriangle(l, m, n)]

    print(a)

    def accept(g):
        z = g(0.)
        return abs(z) < 0.99

    G = mulclose([a, b], accept=accept)
    print(len(G))

    cvs = Canvas()

    z_red = 0.
    p_red = Point([z_red], st=[red])

    #for item in get_orbit(G, p_red):
    #    print(item)

#    for g in G:
#        z,w = g.fixed()
#        #if EPSILON < abs(z) < 0.5 and abs(z.imag)<EPSILON: # z_blue
#        if z.imag > EPSILON and z.real > EPSILON and abs(z) < 0.5:
#            print(z)

    z_green = 0.270959736741934+0.13048733193368528j
    p_green = Point([z_green], st=[green])

    z_blue = 0.2660772452600881
    p_blue = Point([z_blue], st=[blue])

    rb_arc = p_red.line(p_blue)

    #arc = Arc([-1, 0, 1], st_stroke=[0.5*thin])
    arc = Bigon([-1, 0, 1])
    #G = list(G)[:100]

    #geometry = Hyperbolic(G, [rb_arc, p_red, p_blue, p_green])
    geometry = Hyperbolic(G, [arc, p_red, p_blue, p_green])

    cvs = Canvas([Scale(10.0)])
    cvs.stroke(path.circle(0,0,1.2), [white])
    cvs.clip(path.circle(0,0,1.))
    cvs.stroke(path.circle(0,0,1.), [0.95*white,2*thin])

    geometry.render(cvs)

    save(cvs, "stereo_hyperbolic.pdf")

    pop()


def test_logo():
    push()

    Point.radius = 0.04
    #Geometry.scale = 4.0

    from bruhat.disc import mktriangle
    (l, m, n, maxsize) = (7,2,3,400)

    # build the rotation group generators
    a, b = [g.todisc() for g in mktriangle(l, m, n)]

    print(a)

    def accept(g):
        z = g(0.)
        return abs(z) < 0.95

    G = mulclose([a, b], accept=accept)
    print(len(G))
    G = list(G)
    a, b = [g.todisc() for g in mktriangle(8, 2, 4)]
    g = a*b*a*Mobius(1,1,0,1).todisc()
    G = [g*h for h in G]
    G1 = G[:100]

    cvs = Canvas()

    z_red = 0.
    p_red = Point([z_red], st=[red])

    z_green = 0.270959736741934+0.13048733193368528j
    p_green = Point([z_green], st=[green])

    z_blue = 0.2660772452600881
    p_blue = Point([z_blue], st=[blue])

    rb_arc = p_red.line(p_blue)

    #arc = Arc([-1, 0, 1], st_stroke=[0.5*thin])
    arc = Bigon([-1, 0, 1], st=[0.5*white+0.1*blue])
    arc0 = Bigon([-1, 0, 1], st=[white], line_width=0.003)
    #G = list(G)[:100]


    cvs = Canvas([Scale(10.0)])
    r = 1.5
    p = path.rect(-r, -r, 2*r, 2*r)
    cvs.fill(p, [black])
    cvs.clip(path.circle(0,0,1.))
    p = path.circle(0,0,1)
    cvs.fill(p, [0.1*white])
    cvs.stroke(p, [0.55*white+0.2*blue,1*thin])

    cvs.rotate(1.2*2*pi/4)
    geometry = Hyperbolic(G1, [arc, arc0]).render(cvs)
    geometry = Hyperbolic(G, [p_red, p_blue, p_green]).render(cvs)

    save(cvs, "stereo_logo.pdf")
    pop()


def test_render():
    print("\ntest_render:")
    test()
    test_circles()
    test_stereo()
    test_worksheet()
    test_wallpaper()

    test_hyperbolic()
    

if __name__ == "__main__":

    from time import sleep, time
    start_time = time()
    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))

        
