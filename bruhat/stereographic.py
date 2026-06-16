#!/usr/bin/env python3

"""

Render circles and lines under Mobius group actions.

See:
"An Introduction to the McKay Correspondence
Master Thesis in Physics"
Max Lindh

TODO: 
shade regions bounded by arc and lines.

"""

import sys
from random import randint, seed, choice, random
from functools import reduce, cache
from operator import mul
from cmath import pi, exp
from math import atan

import numpy

from bruhat.util import cross
from bruhat.smap import SMap
from bruhat.argv import argv

from bruhat.todd_coxeter import Schreier
from bruhat.mobius import Mobius, mulclose, EPSILON
from bruhat.disc import get_pos_angle


try:
    from huygens.namespace import (
        Canvas, red, green, blue, grey, black, white, path,
        st_THick, st_center, 
        arc_to_bezier,)
    
    from huygens import config
    config(text="pdflatex", latex_header=r"""
    \usepackage{amsmath}
    \usepackage{amssymb}
    """)
except:
    print("\n\nWARNING: huygens not found\n\n")


rnd = lambda r=1: (1-2*random())*r + (1-2*random())*r*1j
metric = lambda w,z=0. : ((w-z) * complex(w-z).conjugate()).real**0.5
eq = lambda z0, z1: abs(complex(z0-z1)) < EPSILON

infty = None
I = Mobius()

def zpromote(z):
    return complex(z) if z is not None else z

def zeq(z0, z1):
    if z0==z1:
        return True
    if z0 is None or z1 is None:
        return False
    return abs(z0-z1)<EPSILON


class Feature:
    RADIUS = 4.0 

    def __init__(self, zs):
        if zs is None or isinstance(zs, (float, complex, int)):
            zs = [zs]
        zs = [zpromote(z) for z in zs]
        self.zs = zs

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.zs,)

    def __rmul__(self, g):
        assert isinstance(g, Mobius)
        zs = [g(z) for z in self.zs]
        return self.__class__(zs)


class Point(Feature):
    def __eq__(self, other):
        assert isinstance(other, Point)
        z0, z1 = self.zs[0], other.zs[0]
        return zeq(z0, z1)
        
    def render(self, cvs, radius, st):
        z = self.zs[0]
        if z is None:
            pass # ?
        else:
            cvs.fill(path.circle(z.real, z.imag, radius), st)


class Circle(Feature):
    """ A circle (or line) in the complex plane defined by three points.
    """ 
    def __init__(self, zs):
        assert len(zs) == 3
        Feature.__init__(self, zs)
        a, b, c = self.zs
        assert not zeq(a, b)
        assert not zeq(a, c)
        assert not zeq(b, c)

    def __eq__(self, other):
        assert isinstance(other, Circle)
        l, r = self.is_colinear(), other.is_colinear()
        if l and r:
            a, b, c = self.zs
            d, e, f = other.zs
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
        assert abs(d)>EPSILON, "is_colinear"
    
        ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
        uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    
        return ux + 1j*uy

    def get_radius(self):
        assert not self.is_colinear()
        z0 = self.get_center()
        return metric(z0, self.zs[0])

    def render(self, cvs, st_stroke=[], st_fill=None):
        if self.is_colinear():
            a, b, c = self.zs
            dz = b-a
            dz = dz/metric(dz)
            #z0, z1 = a-100*dz, a+100*dz
            # XXX check we go through zero XXX
            z0, z1 = -self.RADIUS*dz, self.RADIUS*dz
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)
        
        else:
            z0 = self.get_center()
            r = self.get_radius()
            p = path.circle(z0.real, z0.imag, r)
        if st_fill is not None:
            cvs.fill(p, st_fill)
        if st_stroke is not None:
            cvs.stroke(p, st_stroke)


_det = lambda a,b,c,d:a*d-b*c
det = lambda z,w:_det(z.real, z.imag, w.real, w.imag)

class Arc(Circle):
    "directed arc of a circle defined by three points: start, (random) middle, end"
    def __eq__(self, other):
        assert isinstance(other, Arc)
        if not Circle.__eq__(self, other):
            return False
        # FIX FIX TODO TODO FIX FIX
        return (zeq(self.zs[0], other.zs[0]) and zeq(self.zs[2], other.zs[2])
        or zeq(self.zs[0], other.zs[2]) and zeq(self.zs[2], other.zs[0])) # no ...
        
    def render(self, cvs, st_stroke=[], st_fill=None, debug=False):
        zs = self.zs
        if None in zs:
            # arc to infty
            assert zs.count(None)==1, zs
            assert zs[1] is not None, zs # ouch: "middle" should not be infty !?!
            z0 = zs[0] if zs[0] is not None else zs[2] # start here
            z1 = z0*Feature.RADIUS/abs(z0) # infty is here
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)
            cvs.stroke(p, st_stroke)
            return # <---------------- return

        if self.is_colinear():
            z0, _, z1 = self.zs
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)
            cvs.stroke(p, st_stroke)
            return # <---------------- return

        zc = self.get_center()
        r = self.get_radius()
        assert r>EPSILON
        xc, yc = zc.real, zc.imag
        z0, z1, z2 = self.zs
        dz0, dz1, dz2 = [z-zc for z in self.zs]
        dz0, dz1, dz2 = [1, dz1/dz0, dz2/dz0]
        if get_pos_angle(dz1) > get_pos_angle(dz2):
            #print("swap")
            z0, z2 = z2, z0
        theta0 = get_pos_angle(z0-zc)
        theta1 = get_pos_angle(z1-zc)
        theta2 = get_pos_angle(z2-zc)
        #print(theta0, theta1, theta2)
        p = arc_to_bezier(xc, yc, r, theta0, theta2)
        cvs.stroke(p, st_stroke)

        if debug:
            for z,c in zip(self.zs, [red, blue, green]):
                cvs.fill(path.circle(z.real, z.imag, 0.05), [c])



def get_BT():
    # Binary tetrahedral group (see Lindh p32)
    r2 = 2**0.5
    r3 = 3**0.5
    i = 1j
    a = Mobius((1+i*r3)/2, 0, 0, (1-i*r3)/2)
    assert a.order() == 6
    s = 2*r3
    b = Mobius((r3-i)/s, (-i*2*r2)/s, (-i*2*r2)/s, (r3+i)/s)
    assert b.order() == 6

    G = mulclose([a,b])
    assert len(G) == 24
    print("get_BT:", len(G))
    return G


def get_BO():
    # Binary octahedral group (Lindh p34)
    r2 = 2**0.5
    i = 1j
    a = Mobius((1+i)/r2, 0, 0, (1-i)/r2)
    b = Mobius(1/r2, -i/r2, -i/r2, 1/r2)
    G = mulclose([a,b])
    assert len(G) == 48
    print("get_BO:", len(G))
    return G


def get_BD():
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

    G = mulclose([a,b])
    assert len(G) == 120
    print("get_BD:", len(G))
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

    for trial in range(10):
        c = Circle([rnd() for i in range(3)])

        z0 = c.get_center()

        r = None
        for z in c.zs:
            r1 = metric(z, z0)
            assert r is None or abs(r-r1) < EPSILON


def render():

    scale = 3.0
    radius = 0.06

    G = get_BT()

    R = Feature.RADIUS
    cvs = Canvas().scale(scale)
    p = path.circle(0, 0, R)
    cvs.stroke(p, [green])
    cvs.stroke(path.circle(0,0,1.1*R), [white])
    cvs.clip(p)
    
    cl = 0.6*white
    #c = Circle([rnd() for i in range(3)])
    c = Circle([0, 1, 2])
    orbit = get_orbit(G, c)
    for item in orbit:
        item.render(cvs, [cl]+st_THick)

    for c in orbit:
        if c.is_colinear():
            continue
        z0 = c.get_center()
        if abs(z0.imag) < EPSILON:
            r = c.get_radius()
            z_blue = z0+r

    #print("z_blue:", z_blue.real)
    #return

#    z = Point(infty)
#    for g in G:
#        z1 = (g*z).zs[0]
#        if z1 is None:
#            continue
#        if abs(z1.imag) < EPSILON and z1.real < 0:
#            print("green:", z1.real**2)
#    return

    for point,cl in [(Point(infty),green), (Point(0), red), (Point(z_blue), blue)]:
        orbit = get_orbit(G, point)
        print("orbit:", len(orbit))
        for point in orbit:
            point.render(cvs, radius, [cl])

    #for g in G:
    #    z = g(0.2 + 0.1j)
    #    cvs.text(z.real, z.imag, r"$\star$", st_center)

    # -----------------------------------------------------

    cvs = Canvas([cvs]).show_page().scale(scale)

    G = get_BO()
    p = path.circle(0, 0, R)
    cvs.stroke(p, [green])
    cvs.stroke(path.circle(0,0,1.1*R), [white])
    cvs.clip(p)
    
    diag = lambda x:x+x*1j

    c0 = Circle([diag(x) for x in [0.5, 1.1, 2]])
    c1 = Circle([0.5, 1.1, 2])

    orbit = get_orbit(G, c0) + get_orbit(G, c1)
    print("orbit:", len(orbit))
    cl = 0.6*white
    for item in orbit:
        item.render(cvs, [cl]+st_THick)
        if item.is_colinear():
            continue
        z = item.get_center()
        if abs(z.imag) < EPSILON and z.real < -EPSILON:
            r = item.get_radius()
            z_blue = z-r
        #if abs(z.real) < EPSILON and z.imag > EPSILON:
        #    item.render(cvs, [red]+st_THick)

    for g in G:
        pts = g.fixed()
        a, b = pts
        if abs(a.real-a.imag) < EPSILON and EPSILON<a.real<0.5:
            #Point([a]).render(cvs, radius, [red])
            z_red = a

    orbit = get_orbit(G, Point([0]))
    for item in orbit:
        item.render(cvs, radius, [green])

    orbit = get_orbit(G, Point([z_blue]))
    for item in orbit:
        item.render(cvs, radius, [blue])

    orbit = get_orbit(G, Point([z_red]))
    for item in orbit:
        item.render(cvs, radius, [red])

    # -----------------------------------------------------
    # -----------------------------------------------------

    cvs = Canvas([cvs]).show_page().scale(scale)

    G = get_BD()
    p = path.circle(0, 0, R)
    cvs.stroke(p, [green])
    cvs.stroke(path.circle(0,0,1.1*R), [white])
    cvs.clip(p)
    
    diag = lambda x:x+x*1j

    c = Circle([0.5, 1.1, 2])

    orbit = get_orbit(G, c)
    print("orbit:", len(orbit))
    cl = 0.6*white
    for item in orbit:
        item.render(cvs, [cl]+st_THick)

    items = []
    for g in G:
        #(g*p).render(cvs, radius, [black])
        pts = g.fixed()
        a, b = pts
        items.append(Point([a]))
        items.append(Point([b]))
        #Point([a]).render(cvs, radius, [green])
        #Point([b]).render(cvs, radius, [green])

    found = unique(items)

    orbit = get_orbit(G, Point([0]))
    for item in orbit:
        item.render(cvs, radius, [green])
        if item in found:
            found.remove(item)

    print(len(found))
    for item in found:
        z = item.zs[0]
        if z.real>2 and abs(z.imag)<EPSILON:
            #item.render(cvs, radius, [red])
            z_blue = z
        if z.real< -2 and abs(z.imag)<EPSILON:
            #item.render(cvs, radius, [red])
            z_red = z

    orbit = get_orbit(G, Point([z_blue]))
    for item in orbit:
        item.render(cvs, radius, [blue])

    orbit = get_orbit(G, Point([z_red]))
    for item in orbit:
        item.render(cvs, radius, [red])

    cvs.writePDFfile("Circles.pdf")


def test_arc():

    scale = 3.0
    radius = 0.06

    R = Feature.RADIUS
    cvs = Canvas().scale(scale)
    p = path.circle(0, 0, R)
    cvs.stroke(p, [green])
    cvs.stroke(path.circle(0,0,1.1*R), [white])
    cvs.clip(p)

    i = 1j
    item = Arc([0, i, 1+i])
    #item = Arc([i+1, i, 0])
    #item = Arc([1, -1, -i])
    #item.render(cvs, [black], debug=True)

    z_red = 0.0
    z_blue = 0.5176380902050413
    #z_green = (-1/2**0.5)*exp(4*pi*i/3)
    z_green = (-1/2**0.5)

    rb_item = Arc([z_red, z_blue/2, z_blue])
    rg_item = Arc([z_red, z_green/2, z_green])
    z_blue = z_green - (z_blue-z_green)
    bg_item = Arc([z_blue, 2*z_blue, infty])

    G = get_BT()

    for item in [rb_item, rg_item, bg_item]:
        found = unique([g*item for g in G])
        print("found:", len(found))
        for item in found:
            item.render(cvs, [grey]+st_THick, debug=False)

    cvs.writePDFfile("test_arc.pdf")

    

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

        
