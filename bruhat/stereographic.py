#!/usr/bin/env python3

"""

"""

import sys
from random import randint, seed, choice, random
from functools import reduce, cache
from operator import mul
from math import pi

import numpy

from bruhat.util import cross
from bruhat.smap import SMap
from bruhat.argv import argv

from bruhat.todd_coxeter import Schreier
from bruhat.mobius import Mobius, mulclose, EPSILON


try:
    from huygens.namespace import *
    
    #from huygens import config
    #config(text="pdflatex", latex_header=r"""
    #\usepackage{amsmath}
    #\usepackage{amssymb}
    #""")
except:
    print("huygens not found")


rnd = lambda r=1: (1-2*random())*r + (1-2*random())*r*1j
metric = lambda w,z=0. : ((w-z) * complex(w-z).conjugate()).real**0.5
eq = lambda z0, z1: abs(complex(z0-z1)) < EPSILON

infty = None
I = Mobius()

class Circle:
    """ A circle (or line) in the complex plane defined by three points.
    """ 
    def __init__(self, zs):
        assert len(zs) == 3
        zs = [complex(z) for z in zs]
        a, b, c = zs
        assert abs(a-b) > EPSILON, zs
        assert abs(a-c) > EPSILON, zs
        assert abs(b-c) > EPSILON, zs
        self.zs = zs

    def __str__(self):
        return "Circle(%s)"%(self.zs,)

    def __rmul__(self, g):
        assert isinstance(g, Mobius)
        zs = [g(z) for z in self.zs]
        return Circle(zs)

#    def gop(self, other):
#        assert isinstance(other, Circle)
#        m = Mobius.send(*self.zs, *other.zs)
#        a = Mobius.send(*self.zs, 1, 1j, -1)
#        g = a*m*(~a)
#        return g

    def failed_eq(self, other):
        assert isinstance(other, Circle)
        m = Mobius.send(*self.zs, *other.zs)
        a = Mobius.send(*self.zs, 1, 1j, -1)
        g = a*m*(~a)
        return g.is_su() # nope

    def __eq__(self, other):
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
            z0, z1 = a-100*dz, a+100*dz
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)
        
        else:
            z0 = self.get_center()
            r = self.get_radius()
            p = path.circle(z0.real, z0.imag, r)
        if st_fill is not None:
            cvs.fill(p, st_fill)
        if st_stroke is not None:
            cvs.stroke(p, st_stroke)



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

    # Binary octahedral group (see Lindh p32)
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

    R = 4.
    cvs = Canvas()
    p = path.circle(0, 0, R)
    cvs.stroke(p, [green])
    cvs.clip(p)
    
    #c = Circle([rnd() for i in range(3)])
    c = Circle([0, 1, 2])
    found = []
    for g in G:
        pts = g.fixed()
        gc = g*c
        if gc in found:
            continue
        #d = Circle([z+0.04*random() for z in gc.zs])
        gc.render(cvs, [black.alpha(0.2)]+st_THick)
        if not gc.is_colinear():
            z0 = gc.get_center()
            r0 = gc.get_radius()
            for other in found:
                if other.is_colinear():
                    continue
                z1 = other.get_center()
                r1 = other.get_radius()
                if abs(z0-z1)<EPSILON and abs(r0-r1)<EPSILON:
                    print(gc)
                    print(other)
                    gop = gc.gop(other)
                    print(gop, gop.is_su())
                    gop.dump_su()
                    assert 0
            
        found.append(gc)
    print(len(found), len(G))

    for c in found:
        if c.is_colinear():
            continue
        z0 = c.get_center()
        if abs(z0.imag) < EPSILON:
            r = c.get_radius()
            z_blue = z0+r

    for z0,cl in [(infty,green), (0, red), (z_blue, blue)]:
        for g in G:
            z = g(z0)
            if z is not infty:
                cvs.fill(path.circle(z.real, z.imag, 0.05), [cl])
    print()

    cvs.writePDFfile("Circles.pdf")



    

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

        
