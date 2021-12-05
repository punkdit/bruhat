#!/usr/bin/env python3

"""
Visualize Poincare disc model of hyperbolic space.

copied from transversal2018.hyperbolic


"""

import sys
from math import pi, hypot, cos, sin, tan, acosh, atan
import cmath

from huygens import config
config(text="pdflatex", latex_header=r"""
\usepackage{amsmath}
\usepackage{amssymb}
""")

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
    pass


class CircleGeodesic(Geodesic):
    def __init__(self, z_center, radius, theta0, theta1):
        self.z_center = complex(z_center)
        assert abs(self.z_center) > 1+EPSILON, " not a geodesic "
        self.radius = radius
        self.theta0 = theta0
        self.theta1 = theta1


class LineGeodesic(Geodesic):
    def __init__(self, z0, z1):
        self.z0 = complex(z0)
        self.z1 = complex(z1)


class Disc(object):

    def __init__(self, z_center=1j):
        self.g_center = Mobius.cayley(z_center) * ~Mobius.cayley()

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
            cvs.stroke(path.circle(x, y, 0.02), [cl])
            cl = green

    def d_poincare(self, z):
        g = self.g_center
        return d_poincare(g(z))

    def show_point(self, z, attrs=[]):
        z = self.g_center(z)
        cvs.stroke(path.circle(z.real, z.imag, 0.04), attrs)

    def get_short_geodesic(self, z0, z1):
        z0, z1 = self.g_center(z0), self.g_center(z1)
        #print("show_short_geodesic", z0, z1)
        g = Mobius.get_translate(z0, z1)
        #self.show_fixed(g)
        l, r = g.fixed()
        #self.show_geodesic(l, r)
        src, tgt = l, r

        theta0 = get_angle(src)
        theta1 = get_angle(tgt)
        if abs(src + tgt) < EPSILON:
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)
            return p # <---------------- return

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
            #dbg = meths[idx](z_center.real, z_center.imag, radius,
            #    conv(0.9, theta0, theta1), theta1)
            #cvs.stroke(dbg, [black]+st_THick)
            #print("HERE")
        #if theta0 > theta1:
        #    theta0, theta1 = theta1, theta0
        #    idx = 1-idx
        #if abs(z2) > 1.:
        #    idx = 1-idx
        #    theta0, theta1 = theta1, theta0

        #print("show_short_geodesic")
        #print(z_center, radius, theta0, theta1)
        #print(abs(z2) < 1.)

        p = meths[idx](z_center.real, z_center.imag, radius, theta0, theta1)
        return p

    def show_short_geodesic(self, z0, z1, attrs=[]):
        #print("show_short_geodesic", z0, z1)
        p = self.get_short_geodesic(z0, z1)
        cvs.stroke(p, attrs)

    def show_polygon(self, pts, st_stroke=[], st_fill=None):
        n = len(pts)
        items = [self.get_short_geodesic(pts[i], pts[(i+1)%n])
            for i in range(n)]
        items.append(path.closepath())
        p = path.path(items)
        if st_fill is not None:
            cvs.fill(p, st_fill)
        cvs.stroke(p, st_stroke)

    def show_geodesic(self, z0, z1, attrs=[]):
        "show_geodesic between two points"

        if (abs(z0*z0.conjugate() - 1.) > EPSILON or
            abs(z1*z1.conjugate() - 1.) > EPSILON):
            self.show_short_geodesic(z0, z1, attrs)
            return # <<<<<<<<<<-- return

        assert abs(z0*z0.conjugate() - 1.) < EPSILON
        assert abs(z1*z1.conjugate() - 1.) < EPSILON
        z0, z1 = self.g_center(z0), self.g_center(z1)

        #from huygens.pov import get_angle
        #z = z0/z1
        #theta = get_angle(z.real, z.imag)
        #print(theta*360 / (2*pi))

        if abs(z0 + z1) < EPSILON:
            p = path.line(z0.real, z0.imag, z1.real, z1.imag)

        else:
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
            #cvs.stroke(path.line(z0.real, z0.imag, z1.real, z1.imag))
            #print("arc", radius, t0, t1)
            #cvs.stroke(path.line(0., 0., z_center.real, z_center.imag))
            p = path.arc(z_center.real, z_center.imag, radius, t0, t1)

        cvs.stroke(p, attrs)

    #def get_refl(self, z0, z1):

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
            #cvs.fill(path.circle(x1, y1, 0.01), [green])

            g2 = self.g_center*g

            edge1 = [g2.trafo(*p) for p in edge]
            for idx in range(len(edge1)-1):
                p, q = edge1[idx:idx+2]
                lw = 2./d_poincare(p[0] + 1.j*p[1])
                cvs.stroke(path.line(*p, *q), [lw*thin]+st_round)

            z = g2(z2)
            ds = 2 / (1 - z*z.conjugate()).real
            draw_qubit(z.real, z.imag, 1.5/ds)



# -------------------------------------------------------------
#

print()

def save(name):
    print("save(%r)"%(name,))
    cvs.writePDFfile("i"+"mages/"+name+".pdf")
    cvs.writeSVGfile("i"+"mages/"+name+".svg")
    cvs.writePNGfile("i"+"mages/"+name+".png")


def draw_qubit(x, y, scale=1.0):
    p = path.circle(x, y, 0.05*scale)
    cvs.fill(p, [white])
    cvs.fill(p, [black.alpha(0.1)])
    cvs.stroke(p, [black, normal*scale])


# -------------------------------------------------------------
#


def main():
    global cvs

    gens = mktriangle(5, 2, 5)
    
    gens = [g.todisc() for g in gens] # --------------- todisc --------------
    a, b = gens
    
    I = Mobius()
    
    disc = Disc()

    cvs = Canvas([Scale(5)])

    z0 = cos(0.1) + 1j*sin(0.1)
    z1 = -z0
    disc.show_geodesic(-1, 1, [green])
    disc.show_geodesic(z0, z1, [green])

    z0 = 0j
    z1 = (a*b).inner_fixed() # pentagon center

    for g in [I, a, a**2, a**3, a**4]:
        for h in [g, b*g, b*g*b]:
            disc.show_geodesic(h(z0), h(z1), [thin])
    
    p = path.circle(0, 0, 1)
    cvs.clip(p)
    cvs.stroke(p, [1.0*thin])
    
    save("brings-curve")
    



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
    
    disc = Disc(z1)
    
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
    
    cvs = Canvas([Scale(5)])
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
        #disc.show_fixed(a1) # vertices
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
    
    disc.show_tiles(G)
    
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
    main()







