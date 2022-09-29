#!/usr/bin/env python3

"""
Visualize Poincare disc model of hyperbolic space.

copied from transversal2018.hyperbolic


"""

import sys
from math import pi, hypot, cos, sin, tan, acosh, atan
from functools import reduce
from operator import mul
import cmath

from huygens.namespace import *

from bruhat.mobius import Mobius, mktriangle, mulclose, d_poincare

def get_angle(z):
    z = complex(z)
    r = (z*z.conjugate())**0.5
    assert abs(r) > EPSILON
    z /= r
    theta = cmath.log(z).imag
    return theta


class Geodesic(object):
    @classmethod
    def construct(cls, z0, z1):
        if abs(z0-z1)<EPSILON or \
            abs(z0) < EPSILON or \
            abs(z1) < EPSILON or \
            (abs(z0.real)<EPSILON and abs(z1.real) < EPSILON):
            gamma = LineGeodesic(z0, z1)
        elif abs(z0.real)>EPSILON and abs(z1.real)>EPSILON and \
            abs(z0.imag/z0.real - z1.imag/z1.real) < EPSILON:
            gamma = LineGeodesic(z0, z1)
        elif (abs(z0*z0.conjugate() - 1.) > EPSILON or abs(z1*z1.conjugate() - 1.) > EPSILON):
            gamma = cls._construct_short(z0, z1)
        else:
            gamma = cls._construct_long(z0, z1)
        return gamma

    @classmethod
    def _construct_long(cls, z0, z1):
        assert abs(z0*z0.conjugate() - 1.) < EPSILON
        assert abs(z1*z1.conjugate() - 1.) < EPSILON

        theta0 = get_angle(z0)
        theta1 = get_angle(z1)
        r_cos = 1 / cos(0.5*abs(theta0 - theta1))
        radius = tan(0.5*abs(theta0 - theta1))
        z_center = r_cos * cmath.exp(1j * 0.5*(theta0+theta1))
        #print(theta0*360 / (2*pi))
        #print(theta1*360 / (2*pi))

        radius = abs(radius)
        t0 = (theta0 + 0.5*pi) % (2*pi)
        t1 = (theta1 - 0.5*pi) % (2*pi)
        gamma = CircleGeodesic(z_center, radius, t0, t1, path.arc)

        return gamma

    @classmethod
    def _construct_short(cls, z0, z1):

        #print("_construct_short", z0, z1)
        if abs(z0)>1-EPSILON:
            #z0 *= 0.999 # argh, stupid hack
            g = Mobius.get_translate(0.999*z0, z1)
        elif abs(z1)>1-EPSILON:
            #z1 *= 0.999 # argh, stupid hack
            g = Mobius.get_translate(z0, 0.999*z1)
        else:
            g = Mobius.get_translate(z0, z1)

        l, r = g.fixed()
        src, tgt = l, r

        theta0 = get_angle(src)
        theta1 = get_angle(tgt)
        r_cos = 1 / cos(0.5*abs(theta0 - theta1))
        radius = tan(0.5*abs(theta0 - theta1))
        z_center = r_cos * cmath.exp(1j * 0.5*(theta0+theta1))

        if radius < 0:
            radius = -radius

        theta0 = get_angle(z0-z_center) % (2*pi)
        theta1 = get_angle(z1-z_center) % (2*pi)
        theta2 = 0.5*(theta0 + theta1)
        z2 = z_center + radius*cmath.exp(1j*theta2)
        meths = [path.arc, path.arcn]
        idx = 0
        if theta0 > theta1:
            idx = 1
        #print(theta1-theta0)
        if abs(theta1-theta0) > pi:
            theta1 = theta1 - 2*pi
            idx = 1-idx

        gamma = CircleGeodesic(z_center, radius, theta0, theta1, meths[idx])
        assert abs(gamma.z0 - z0) < 1e-4, (gamma.z0, z0)
        assert abs(gamma.z1 - z1) < 1e-4, (gamma.z1, z1)

        return gamma

    def __init__(self, z0, z1, z2):
        self.z0 = complex(z0)
        self.z1 = complex(z1)
        self.z2 = complex(z2) # midpoint

    def get_refl(self):
        z0, z1, z2 = self.z0, self.z1, self.z2
        g = Mobius.send(z0, z1, z2, -1, 0., 1.)
        C = Mobius(1, 0, 0, 1, True) # conjugate
        return (~g)*C*g


class CircleGeodesic(Geodesic):
    def __init__(self, z_center, radius, theta0, theta1, method=path.arc):
        self.z_center = complex(z_center)
        #assert abs(self.z_center) > 1+EPSILON, " not a geodesic: %s "%(abs(self.z_center),)
        assert abs(self.z_center) > 1, " not a geodesic: %s "%(abs(self.z_center),)
        self.radius = radius
        self.theta0 = theta0
        self.theta1 = theta1
        self.p = method(z_center.real, z_center.imag, radius, theta0, theta1)
        z0 = radius*(cos(theta0) + 1j*sin(theta0)) + z_center
        z1 = radius*(cos(theta1) + 1j*sin(theta1)) + z_center
        theta2 = 0.5*(theta0 + theta1)
        z2 = radius*(cos(theta2) + 1j*sin(theta2)) + z_center
        z3 = -radius*(cos(theta2) + 1j*sin(theta2)) + z_center
        lhs = abs(z0-z2) + abs(z1-z2)
        rhs = abs(z0-z3) + abs(z1-z3)
        if lhs > rhs:
            z2 = z3
        Geodesic.__init__(self, z0, z1, z2)



class LineGeodesic(Geodesic):
    def __init__(self, z0, z1):
        self.z0 = z0 = complex(z0)
        self.z1 = z1 = complex(z1)
        self.p = path.line(z0.real, z0.imag, z1.real, z1.imag)
        z2 = 0.5*(z0 + z1)
        Geodesic.__init__(self, z0, z1, z2)
        if abs(z0) < EPSILON or abs(z1)<EPSILON:
            return
        t0 = get_angle(z0)
        t1 = get_angle(z1)
        if abs(t0-t1) > EPSILON:
            print("LineGeodesic fail:", z0, z1, t0, t1)


class Disc(object):

    def __init__(self, cvs=None, z_center=1j, scale=5.0, bg=white):
        self.g_center = Mobius.cayley(z_center) * ~Mobius.cayley()
        if cvs is None:
            cvs = Canvas()
        cvs.append(Scale(scale))
        if bg is not None:
            cvs.fill(path.rect(-1, -1, 2, 2), [bg])
        self.cvs = cvs

    def show_fixed(self, g):
        #if g.trace.real < 2.-EPSILON and abs(g.c)>EPSILON:
        cl = red
        prev = 999.
        for z in g.fixed():
            if abs(z) > 1+EPSILON:
                continue
            if abs(z-prev) < EPSILON:
                continue
            prev = z
            #z = g.inner_fixed()
            w = g(z)
            assert abs(w-z) < EPSILON
            z = self.g_center(z)
            x, y = z.real, z.imag
            self.cvs.stroke(path.circle(x, y, 0.02), [cl])
            cl = green

    def d_poincare(self, z):
        g = self.g_center
        return d_poincare(g(z))

    def show_point(self, z, fill_attrs=None, stroke_attrs=None, radius=0.05):
        z = self.g_center(z)
        r = 2./d_poincare(z)
        p = path.circle(z.real, z.imag, radius*r)
        if fill_attrs is not None:
            self.cvs.fill(p, fill_attrs)
        if stroke_attrs is not None:
            self.cvs.stroke(p, stroke_attrs + [normal*r])

    def get_geodesic(self, z0, z1):
        z0, z1 = self.g_center(z0), self.g_center(z1)
        #print("show_short_geodesic", z0, z1)
        gamma = Geodesic.construct(z0, z1)
        return gamma

    def show_polygon(self, pts, st_stroke=[], st_fill=None): # XXX broken XXX
        cvs = self.cvs
        n = len(pts)
        z = sum(pts) / n
        # sometimes these geodesics are reversed XXX
        items = [self.get_geodesic(pts[i], pts[(i+1)%n]) for i in range(n)]
#        for item in items:
#            if isinstance(item, LineGeodesic):
#                #print(item.z0, item.z1)
#                #st_fill = [grey]
#                cvs.stroke(item.p, [thin])
#            else:
#                cvs.stroke(item.p, [thin])
#        return
        items = [item.p for item in items]
        items.append(path.closepath())
        p = path.path(items)
        if st_fill is not None:
            self.cvs.fill(p, st_fill)
        r = 2./d_poincare(z)
        #if isinstance(p[0], path.line):
        cvs.stroke(p, st_stroke+[r*normal])

    def show_geodesic(self, z0, z1, attrs=[]):
        "show_geodesic between two points"
        z0, z1 = self.g_center(z0), self.g_center(z1)
        gamma = Geodesic.construct(z0, z1)
        r = 2./d_poincare(gamma.z2)
        self.cvs.stroke(gamma.p, attrs+[r*normal])

    def show_tiles(self, G):
        # -----------------------------------------------------------
        # draw pentagons & qubits
        N = 5
        #edge = [(conv(j/(N-1), 0., z3.real), 0.) for j in range(N)]
        edge = [(conv(j/(N-1), 0., z2.real), 0.) for j in range(N)]
        for g in G:
            #print("trace:", g.trace, "fixed:", g.fixed())
            #assert g.is_sl()
            assert g.is_su()
            #self.cvs.fill(path.circle(x1, y1, 0.01), [green])

            g2 = self.g_center*g

            edge1 = [g2.trafo(*p) for p in edge]
            for idx in range(len(edge1)-1):
                p, q = edge1[idx:idx+2]
                lw = 2./d_poincare(p[0] + 1.j*p[1])
                self.cvs.stroke(path.line(*p, *q), [lw*thin]+st_round)

            z = g2(z2)
            ds = 2 / (1 - z*z.conjugate()).real
            self.show_qubit(z.real, z.imag, 1.5/ds)

    def show_qubit(self, z, scale=None):
        cvs = self.cvs
        x, y = z.real, z.imag
        if scale is None:
            scale = 1.5/d_poincare(z)
        p = path.circle(x, y, 0.05*scale)
        cvs.fill(p, [white])
        cvs.fill(p, [black.alpha(0.1)])
        cvs.stroke(p, [black, normal*scale])

    def show_label(self, z, label, scale=1.0):
        cvs = self.cvs
        x, y = z.real, z.imag
        scale = scale*2.5/d_poincare(z)
        if scale > 1e-5:
            cvs.text(x, y, label, [Scale(scale)]+st_center)

    def fini(self):
        cvs = self.cvs
        p = path.circle(0, 0, 1)
        cvs.clip(p)
        cvs.stroke(p, [1.5*thin]+[black.alpha(0.5)])
        cvs.stroke(p, [1.0*thin])
    
    def save(self, name):
        print("save(%r)"%(name,))
        cvs = self.cvs
        cvs.writePDFfile("i"+"mages/"+name+".pdf")
        cvs.writeSVGfile("i"+"mages/"+name+".svg")
        cvs = Canvas([Scale(3.), cvs])
        cvs.writePNGfile("i"+"mages/"+name+".png")
    

# -------------------------------------------------------------
#


def main_bring_1():
    # build the rotation group generators
    a, b = [g.todisc() for g in mktriangle(5, 2, 5)]
    
    disc = Disc()

    z_face = 0j
    z_vert = (a*b).inner_fixed()

    gamma = Geodesic.construct(z_vert, (~a)(z_vert))
    z_edge = gamma.z2 # midpoint
    g_face = gamma.get_refl()
    gamma = Geodesic.construct(z_face, z_vert)
    g_edge = gamma.get_refl()
    g_vert = Mobius.conjugate()

    gens = [g_face, g_edge, g_vert]
    gens = gens + [~g for g in gens]
    G = mulclose(gens, verbose=True, maxsize=8000) # 8000 is enough

    faces, edges, verts = [], [], []
    for g in G:
        #disc.show_geodesic(g(z_face), g(z_vert), attrs=[grey])
        #disc.show_geodesic(g(z_face), g(z_edge), attrs=[grey])
        faces.append(g(z_face))
        edges.append(g(z_edge))
        verts.append(g(z_vert))

    def key1(z):
        s = "(%.6f,%.6f)"%(z.real, z.imag)
        s = s.replace("-0.000", "0.000")
        return s
    def key(item):
        a, b = item
        k = [key1(a), key1(b)]
        return k[0]+k[1]
    
    def quantize(items):
        items = dict((key(item), item) for item in items)
        return list(items.values())

    items = [(g(z_vert), g(z_edge)) for g in G]
    print(len(items))
    items = quantize(items)
    print(len(items))

#    for g in G:
#        disc.show_geodesic(g(z_vert), g(z_edge), attrs=st_round)
#        break
    for item in items:
        disc.show_geodesic(*item, attrs=st_round)

    #disc.show_polygon([z_face, z_edge, z_vert], st_fill=[grey])
    #disc.show_point(z_vert, radius=0.1)
    #disc.show_point((~a)(z_vert), radius=0.1)
    #disc.show_point(z_edge, radius=0.1)

#    z1, z2 = z_vert, a(z_vert)
#    faces.append(0.j)
#    for i in range(5):
#        disc.show_geodesic(z1, z2)
#        gamma = Geodesic.construct(z1, z2)
#        disc.show_geodesic(0, z1)
#        disc.show_geodesic(0, gamma.z2)
#        edges.append(gamma.z2)
#        verts.append(gamma.z0)
#        #g = gamma.get_refl()
#        #disc.show_point(g(0))
#        z1, z2 = z2, a(z2)

    #for z in edges:
    #    disc.show_qubit(z)
    for z in verts:
        disc.show_qubit(z)

    #for [cl, zs] in ([green, faces], [blue, edges], [red, verts]):
    #    for z in zs:
    #        disc.show_point(z, [cl])

    I = Mobius()
    z = 0.j
    pts0 = []
    for (g, label) in [
        (I, "1"),
        (b, "6"),
        (a*b, "2"),
        (a*a*b, "3"),
        (a*a*a*b, "4"),
        (a*a*a*a*b, "5"),
        (b*a*b, "8"),
        (a*b*a*b, "10"),
        (a*a*b*a*b, "7"),
        (a*a*a*b*a*b, "9"),
        (a*a*a*a*b*a*b, "11"),
        (b*a*b*a*b, "7"),
        (a*b*a*b*a*b, "9"),
        (a*a*b*a*b*a*b, "11"),
        (a*a*a*b*a*b*a*b, "8"),
        (a*a*a*a*b*a*b*a*b, "10"),
        (a*b*a*a*b, "5"),
        (a*a*b*a*a*b, "6"),
        (a*a*a*b*a*a*b, "2"),
        (a*a*a*a*b*a*a*b, "3"),
        (a*a*a*a*a*b*a*a*b, "4"),
        (b*a*a*a*b, "3"),
        (a*b*a*a*a*b, "4"),
        (a*a*b*a*a*a*b, "5"),
        (a*a*a*b*a*a*a*b, "6"),
        (a*a*a*a*b*a*a*a*b, "2"),
        (b*a*a*a*a*b*a*b*a*b, "7"),
        (a*b*a*a*a*a*b*a*b*a*b, "9"),
        (a*a*b*a*a*a*a*b*a*b*a*b, "11"),
        (a*a*a*b*a*a*a*a*b*a*b*a*b, "8"),
        (a*a*a*a*b*a*a*a*a*b*a*b*a*b, "10"),
        (b*a*a*a*a*b*a*b, "12"),
        (a*b*a*a*a*a*b*a*b, "12"),
        (a*a*b*a*a*a*a*b*a*b, "12"),
        (a*a*a*b*a*a*a*a*b*a*b, "12"),
        (a*a*a*a*b*a*a*a*a*b*a*b, "12"),
        (b*a*b*a*a*a*a*b, "12"),
        (a*b*a*b*a*a*a*a*b, "12"),
        (a*a*b*a*b*a*a*a*a*b, "12"),
        (a*a*a*b*a*b*a*a*a*a*b, "12"),
        (a*a*a*a*b*a*b*a*a*a*a*b, "12"),
        (b*a*a*b*a*b, "10"),
        (a*b*a*a*b*a*b, "7"),
        (a*a*b*a*a*b*a*b, "9"),
        (a*a*a*b*a*a*b*a*b, "11"),
        (a*a*a*a*b*a*a*b*a*b, "8"),
    ]:
        z1 = g(z)
        z1 = 0.98*z1
        disc.show_label(z1, label)
        #pts0.append(z1)
        #disc.show_point(z1, [red], radius=0.02)

    disc.fini()
    disc.save("brings-curve-1")
    

def main():

    w, h = 5, 3
    dx = w/3.

    for (l, m, n,gens ) in [(5,2,5,"abc"), (5,2,4,"def")]:

        print(l, m, n)

        cvs = Canvas()

        radius = 0.14
        x, y = 0., 0.

        xs = [dx*i for i in range(3)]

        cvs.stroke(path.line(xs[0], y, xs[2], y), st_Thick)

        for x,cl,label in zip(xs, [red, blue, green], "vertex edge face".split()):
            cvs.fill(path.circle(x, y, radius), [cl])
            cvs.text(x, y-2*radius, label, st_north)

        for x,gen in zip(xs, gens):
            cvs.text(x, y+2*radius, gen, st_south)

        cvs.text(0.5*(xs[0] + xs[1]), y+2*radius, "%d"%l, st_center)
        cvs.text(0.5*(xs[1] + xs[2]), y+2*radius, "%d"%n, st_center)

        cvs.writePDFfile("images/coxeter-%d%d%d"%(l,m,n))




def main_poincare():

    scale = argv.get("scale", 1.)
    for (l, m, n, maxsize) in [(5,2,5,8000), (5,2,4,12000)]:

        # build the rotation group generators
        a, b = [g.todisc() for g in mktriangle(l, m, n)]
        
        cvs = Canvas([Scale(scale)])
        if n==4:
            cvs.append(Rotate(pi))
        disc = Disc(cvs)
    
        z_face = 0j
        z_vert = (a*b).inner_fixed()
    
        gamma = Geodesic.construct(z_vert, (~a)(z_vert))
        z_edge = gamma.z2 # midpoint
        g_face = gamma.get_refl()
        gamma = Geodesic.construct(z_face, z_vert)
        g_edge = gamma.get_refl()
        g_vert = Mobius.conjugate()
    
        gens = [g_face, g_edge, g_vert]
        gens = gens + [~g for g in gens]
        G = mulclose(gens, verbose=True, maxsize=maxsize)
    
        faces, edges, verts = [], [], []
        for g in G:
            faces.append(g(z_face))
            edges.append(g(z_edge))
            verts.append(g(z_vert))
    
        for g in G:
            disc.show_geodesic(g(z_vert), g(z_edge), attrs=st_round)
            disc.show_geodesic(g(z_face), g(z_edge), attrs=st_round+[grey])
            disc.show_geodesic(g(z_face), g(z_vert), attrs=st_round+[grey])
    
    
        for [cl, zs] in ([green, faces], [blue, edges], [red, verts]):
            for z in zs:
                disc.show_point(z, [cl])

        disc.fini()
        disc.save("poincare-disc-%d%d%d"%(l,m,n))
    

def main_poincare_1():

    for (l, m, n, maxsize) in [(5,2,4,4000)]:

        # build the rotation group generators
        a, b = [g.todisc() for g in mktriangle(l, m, n)]
        
        cvs = Canvas()
        cvs.append(Scale(2.))
        disc = Disc(cvs)
    
        z_face = 0j
        z_vert = (a*b).inner_fixed()
    
        gamma = Geodesic.construct(z_vert, (~a)(z_vert))
        z_edge = gamma.z2 # midpoint
        g_face = gamma.get_refl()
        gamma = Geodesic.construct(z_face, z_vert)
        g_edge = gamma.get_refl()
        g_vert = Mobius.conjugate()
    
        gens = [g_face, g_edge, g_vert]
        gens = gens + [~g for g in gens]
        G = mulclose(gens, verbose=True, maxsize=maxsize)
    
        faces, edges, verts = [], [], []
        for g in G:
            faces.append(g(z_face))
            edges.append(g(z_edge))
            verts.append(g(z_vert))
    
        for g in G:
            disc.show_geodesic(g(z_face), g(z_vert), attrs=st_round+[grey.alpha(0.1)])
        for g in G:
            disc.show_geodesic(g(z_face), g(z_edge), attrs=st_round+[grey])
        for g in G:
            disc.show_geodesic(g(z_vert), g(z_edge), attrs=st_round)

    
#        for [cl, zs] in ([green, faces], [blue, edges], [red, verts]):
        for [cl, zs] in ([red, verts],):
            for z in zs:
                disc.show_point(z, [cl])
    
        disc.fini()
        disc.save("poincare-rotation-%d%d%d"%(l,m,n))
    

def render_group(l, m, n, words=[], rels=[], labels={}, name="output", maxsize=1000):

    # build the rotation group generators
    a, b = [g.todisc() for g in mktriangle(l, m, n)]
    assert (a.order()) == 10 # SL(2) not PSL(2)  
    assert (b.order()) == 4  # SL(2) not PSL(2)  
    c = (~b)*(~a)
    
    cvs = Canvas()
    cvs.append(Scale(4.))
    disc = Disc(cvs)

    z_face = 0j
    z_vert = (a*b).inner_fixed()

    gamma = Geodesic.construct(z_vert, (~a)(z_vert))
    z_edge = gamma.z2 # midpoint
    #z_tile = (1/3)*(z_face + z_edge + z_vert)
    z_tile = (1/2)*(z_face + z_vert)
    g_face = gamma.get_refl()
    gamma = Geodesic.construct(z_face, z_vert)
    g_edge = gamma.get_refl()
    g_vert = Mobius.conjugate()

    gens = [g_vert, ~g_vert, g_edge, ~g_edge, g_face, ~g_face]
    G = mulclose(gens, verbose=True, maxsize=maxsize)

    faces, edges, verts = [], [], []
    for g in G:
        faces.append(g(z_face))
        edges.append(g(z_edge))
        verts.append(g(z_vert))

    #for g in G:
    #    disc.show_geodesic(g(z_face), g(z_vert), attrs=st_round+[grey.alpha(0.1)])
    for g in G:
        disc.show_geodesic(g(z_face), g(z_edge), attrs=st_round+[grey])
    for g in G:
        disc.show_geodesic(g(z_vert), g(z_edge), attrs=st_round)

#    for [cl, zs] in ([green, faces], [blue, edges], [red, verts]):
    for [cl, zs] in ([red, verts],):
        for z in zs:
            disc.show_point(z, fill_attrs=[white], stroke_attrs=[black])

    # ------------------------------------------------------

    I = Mobius()
    gens = [a, ~a, b, ~b, c, ~c]
    def get(word):
        g = reduce(mul, [gens[i] for i in word], I)
        return g

    for word in words:
        g = get(word)
        disc.show_point(g(z_tile), [blue.alpha(0.5)])

    scale = 0.4
    if labels:
        for (word, label) in labels.items():
            g = reduce(mul, [gens[w] for w in reversed(word)], I)
            disc.show_label(g(z_tile), label, scale)
            
    else:
        for (g, label) in [
            (I, r"$\star$"), 
            (a, r"$a$"),
            (b, r"$b$"),
            (c, r"$c$"),
            #(c**2, r"$c^2$"),
            #(c**3, r"$c^3$"),
        ]:
            disc.show_label(g(z_tile), label, scale)

    for rel in rels:
        g = get(reversed(rel))
        disc.show_label(g(z_tile), "I", scale)

    disc.fini()
    disc.save(name)
    



# -------------------------------------------------------------
#

def main_bring():
    global cvs
    
    gens = [Mobius(1, 1, 0, 1), Mobius(0, 1, -1, 0)]
    gens = [Mobius(1, 2, 0, 1), Mobius(1, 0, -2, 1)]
    gens = mktriangle(5, 2, 5)
    
    a, b = gens
    z1 = (a*b).inner_fixed() # pentagon center
    gens = [g.todisc() for g in gens] # --------------- todisc --------------
    
    I = Mobius()
    
    cvs = Canvas([Scale(5)])
    disc = Disc(cvs, z1)
    
    a, b = gens
    ab = a*b
    c = ab
    ci = ~c
    assert (ci*ci*a).order() == 5
    assert ab.order() == 5
    
    z0 = a.inner_fixed()
    z1 = c.inner_fixed() # pentagon center
    z2 = b.inner_fixed() # edge center
    z3 = b(z0) # other end of an edge
    
    #G = mulclose(gens, maxsize=2000)
    #G = [Mobius(), a, b, a*b, b*a]
    G = [ab**idx for idx in range(5)]
    
    # cvs.paint([white]) # broken ?!?
    
    aspect = 16/9
    h = 2.2
    w = aspect*h
    #cvs.fill(path.rect(-0.5*w, -0.5*h, w, h), [white])
    
    p = path.circle(0, 0, 1)
    #cvs.fill(p, [grey])
    #cvs.clip(p)
    #cvs.stroke(p)
    
    #G = mulclose([a, b, ~a, ~b], g0=~disc.g_center, maxsize=2000)
    G = mulclose([a, b], maxsize=800) # 8000 is good
    op = ~disc.g_center
    #G = [op * g for g in G]
    
    bring = [I]
    #disc.show_fixed(a)
    edges = []
    for conj in [I, c, c*c, c*c*c, c*c*c*c]:
        b1 = conj*b*~conj
        a1 = conj*a*~conj
        #disc.show_fixed(a1) # _vertices
        #disc.show_fixed(b1) # edge midpoints
        z = b1.inner_fixed()
        edges.append(z)
        tx0, tx1 = (a1*a1*b1)**3, (a1*b1*a1)**3
        #bring += [tx0, ~tx0, tx1, ~tx1]
    
    def shrink(pts, center=None, alpha=0.07):
        if center is None:
            center = sum(pts)/len(pts)
        pts = [conv(alpha, z, center) for z in pts]
        return pts
    
    star = []
    for i in range(5):
        g = a**i
        star.append(g(edges[0]))
    #star = shrink(star, z0)
    #disc.show_polygon(star, [green]+st_round, [green.alpha(0.5)])
    
    edges = shrink(edges, z1)
    star = shrink(star, z0)
    
    def get_translates(w):
        found = []
        for g in G:
            if not g.is_translate():
                continue
            for z in found:
                if abs(z - g(w)) < EPSILON:
                    break
            else:
                yield g
                z = g(w)
                found.append(z)
    
    def show_sigma():
        g = a*ci*ci*a
        src, tgt = g.fixed()
        for i in range(5):
            op = ci**i
            disc.show_geodesic(op(src), op(tgt), [orange.alpha(0.4)])
            break
    
        #print(src)
        src = disc.g_center(src)
        #disc.show_point(src)
        x, y = src.real, src.imag
        #cvs.stroke(path.circle(x, y, 0.1))
        z = src*1j
        dx, dy = 0.2*z.real, 0.2*z.imag
        cvs.text(x+0.2, y-0.1, r"$\sigma$", [Scale(0.5)]+st_center)
    
        r = 1.12
        cvs.stroke(path.line(x, y, r*x, r*y), [orange.alpha(0.4)])
        x, y = 1.1*x, 1.1*y
        st = [deco.earrow(size=0.1)]+[0.7*thin]
        cvs.stroke(path.line(x, y, x+dx, y+dy), st)
        cvs.stroke(path.line(x, y, x-dx, y-dy), st)
    
    
    def show_faces(min_lw = 0.00):
    
        for g in get_translates(z1):
            z = g(z1)
            lw = 4./disc.d_poincare(z)
            if lw < min_lw:
                break
            pts = [g(v) for v in edges]
            disc.show_polygon(pts, [lw*thin, green]+st_round, [green.alpha(0.3)])
    
        #disc.show_point(z0)
    
        for g in get_translates(z0):
            z = g(z0)
            lw = 4./disc.d_poincare(z)
            if lw < min_lw:
                break
            pts = [g(v) for v in star]
            disc.show_polygon(pts, [lw*thin, blue]+st_round, [blue.alpha(0.3)])
    
    #bring = [I]
    #G0 = mulclose([a, b], maxsize=20)
    #bring += [g*tx0*~g for g in G0]
    #bring += [g*tx1*~g for g in G0]
    
    #for op in bring[1:]:
    #    #disc.show_fixed(op)
    #    src, tgt = (op).fixed()
    #    disc.show_geodesic(src, tgt, [green.alpha(0.5)])
    
    #bring = [I, tx0*~tx1*tx0]
    #bring = mulclose([tx0, tx1, ~tx0, ~tx1], maxsize=20)
    #disc.show_fixed(bring[0])
    #bring = mulclose(bring, maxsize=200)
    
    
    def show_bring():
        tiles = [
            (1,  I, ),
            (2,  (c**3)*a, ),
            (3,  (c**2)*a, ),
            (4,  c*a, ),
            (5,  a, ),
            (6,  (c**4)*a, ),
            (7,  a**2, ),
            (8,  a**3, ),
            (9,  ci*a**2, ),
            (10, ci*a**3, ),
            (11, ci**2*a**2, ),
            (7,  ci**2*a**3),
            (8,  ci**3*a**2, ),
            (9,  ci**3*a**3),
            (10, ci**4*a**2, ),
            (11, ci**4*a**3),
        ]
        #g = ci*a**2
        #h = ci**3*a**3
        #disc.show_fixed(h*~g)
        #bring = [I, h*~g]
        g = a**3*b*a**-3
        a1 = g*a*~g
        #disc.show_fixed(a1*a1)
        #assert a1.order() == 5, a1.order()
        tiles.append((11, a1**2*a**3))
        tiles.append((9, a1**4*a**3))
    
        for (count, g) in tiles:
            #if count == 5:
            #    g = (tx0**-3)*g
            for h in bring:
                g2 = disc.g_center*h*g
                z = g2(z1)
                ds = 1./d_poincare(z) ** 0.7
                if ds > 0.05:
                    cvs.text(z.real, z.imag, str(count), [Scale(ds)]+st_center)
                #break
    
    
    
    #show_sigma()
    #show_faces()
    show_bring()
    
    disc.show_tiles(G) # broken XXX
    
    p = path.circle(0, 0, 1)
    cvs.clip(p)
    cvs.stroke(p, [1.0*thin])
    
    #save("hyperbolic-55")
    #save("fold-hyperbolic-55")
    save("brings-curve")
    
    print("OK")
    print()
    
    
    # -------------------------------------------------------------
    


if __name__ == "__main__":
    from huygens import config
    config(text="pdflatex", latex_header=r"""
    \usepackage{amsmath}
    \usepackage{amssymb}
    """)

    from bruhat.argv import argv
    name = argv.next() or "main"
    fn = eval(name)
    fn()

    print("OK\n")







