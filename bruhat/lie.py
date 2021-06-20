#!/usr/bin/env python3

"""
Construct algebraic representations of Lie algebras.
"""

import os
from functools import reduce
import operator
from random import randint, shuffle, seed
from math import sin, cos, pi

import numpy
from scipy.spatial import ConvexHull

from bruhat.argv import argv
from bruhat.util import cross, factorial, choose, determinant
from bruhat.poly import Fraction, Q, Ring, Poly, grobner



def sage_grobner(gens):
    from sage.all_cmdline import PolynomialRing, QQ, ideal
    names = tuple('abcdefghuvwxyz')
    R = PolynomialRing(QQ, order='lex', names=names)
    #(a, b, c, d, u, v, w, x, y, z,) = R._first_ngens(10)
    vs = R._first_ngens(len(names))
    #ns = locals()
    ns = dict(zip(names, vs))
    ps = [eval(str(p).replace("^", "**"), {}, ns) for p in gens]
    if ps == [1]:
        return 1

    basis = ideal(ps)
    return len(basis.groebner_basis())


def all_monomials(vs, deg, ring):
    n = len(vs)
    assert n>0
    if n==1:
        v = vs[0]
        yield v**deg
        return

    items = list(range(deg+1))
    for idxs in cross((items,)*(n-1)):
        idxs = list(idxs)
        remain = deg - sum(idxs)
        if remain < 0:
            continue
        idxs.append(remain)
        p = Poly({():1}, ring)
        assert p==1
        for idx, v in zip(idxs, vs):
            p = p * v**idx
        yield p
            


def test_sl2():

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    a, b, c, d, e, f, g, h = [Poly(c, ring) for c in 'abcdefgh']

    def repr_2d(a, b, c, d):
        A = numpy.empty((2, 2), dtype=object)
        A[:] = [[a, b], [c, d]]
        #A = A.transpose()
        return A

    def repr_3d(a, b, c, d):
        A = numpy.empty((3, 3), dtype=object)
        A[:] = [[a*a, a*c, c*c], [2*a*b, a*d+b*c, 2*c*d], [b*b, b*d, d*d]]
        A = A.transpose()
        return A

    #a, b, c, d, e, f, g, h = [1, 1, 1, 2, 1, 2, 0, 1]

    dot = numpy.dot
    A2 = repr_2d(a, b, c, d)
    B2 = repr_2d(e, f, g, h)
    #print(A2)
    #print(B2)
    AB2 = dot(A2, B2)
    #print(AB2)
    a0, b0, c0, d0 = AB2.flat

    A3 = repr_3d(a, b, c, d)
    B3 = repr_3d(e, f, g, h)
    AB3 = dot(A3, B3)
    #print(A3)
    #print(B3)
    #print(AB3)
    #print(repr_3d(a0, b0, c0, d0))
    assert numpy.alltrue(AB3 == repr_3d(a0, b0, c0, d0))



dim_sl3 = lambda a, b : (a+1)*(b+1)*(a+b+2)//2

def test_hilbert_sl3():
    # ---------------------------------------------------------------
    # Here we work out the coefficients of a Hilbert polynomial
    # given by the rational function top/bot.
    # These coefficients give the dimension of irreps of SL(3).
    # See also test_graded_sl3 below.

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)

    top = one - x*y
    bot = ((one-x)**3) * ((one-y)**3)

    def diff(top, bot, var):
        top, bot = top.diff(var)*bot - bot.diff(var)*top, bot*bot
        return top, bot

    fracs = {}
    fracs[0,0] = (top, bot)

    N = 3
    for i in range(N):
      for j in range(N):
        top, bot = fracs[i, j]
        top, bot = diff(top, bot, 'x')
        fracs[i, j+1] = top, bot

        top, bot = fracs[i, j]
        top, bot = diff(top, bot, 'y')
        fracs[i+1, j] = top, bot
        print(".", end="", flush=True)
    print()

    for i in range(N+1):
      for j in range(N+1):
        if (i, j) not in fracs:
            continue
        top, bot = fracs[i, j]
        t = top.get_const()
        b = bot.get_const()
        val = t//b//factorial(i)//factorial(j)
        assert val == dim_sl3(i, j)
        print(val, end=" ")
      print()


def test_graded_sl3():
    # ---------------------------------------------------------------
    # slightly more explicit calculation than test_hilbert above

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)
    u = Poly("u", ring)
    v = Poly("v", ring)
    w = Poly("w", ring)

    rel = x*u + y*v + z*w
    rels = [rel]
    rels = grobner(rels)

    for a in range(4):
      for b in range(4):
        gens = []
        for p in all_monomials([x, y, z], a, ring):
          for q in all_monomials([u, v, w], b, ring):
            rem = p*q
            #gens.append(pq)
            for rel in rels:
                div, rem = rel.reduce(rem)
                #print(pq, div, rem)
            gens.append(rem)

        basis = grobner(gens)
        assert len(basis) == dim_sl3(a, b)
        print(len(basis), end=' ', flush=True)
      print()


def plot_sl3():
    a = argv.get("a", 2)
    b = argv.get("b", 6)
    weights = get_weights_sl3(a, b)

    weights = dict((k, str(v)) for (k,v) in weights.items())

    name = argv.get("name", "output.pdf")

    #plot_weights(weights)
    diagram = WeightDiagram.make_A2(weights)
    diagram.plot(name)


def plot_all_sl3():

    for a in range(4):
      for b in range(4):
        weights = get_weights_sl3(a, b)
    
        weights = dict((k, str(v)) for (k,v) in weights.items())
    
        name = "sl3_weights_%d%d"%(a, b)
        print(name)
        #plot_weights(weights)
        diagram = WeightDiagram.make_A2(weights)
        diagram.plot(name)


def find_grade(grades, p):
    grade = None
    for m in p.keys:
        g = (0, 0)
        for v, e in m:
            i, j = grades[v]
            g = g[0]+e*i, g[1]+e*j
        assert grade is None or grade==g,\
            "%s is non-homogeneous: %s != %s" % (p, grade, g)
        grade = g
    return grade


def get_weights_sl3(a, b):
    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)
    u = Poly("u", ring)
    v = Poly("v", ring)
    w = Poly("w", ring)

    # See: Fulton & Harris, page 183
    grades = {
        #      L_1, L_3
        'x' : ( 0,   1),
        'y' : (-1,  -1),
        'z' : ( 1,   0),
        'u' : ( 0,  -1),
        'v' : ( 1,   1),
        'w' : (-1,   0),
    }
    rel = x*u + y*v + z*w
    rels = [rel]
    rels = grobner(rels)

    gens = []
    for p in all_monomials([x, y, z], a, ring):
      for q in all_monomials([u, v, w], b, ring):
        rem = p*q
        for rel in rels:
            div, rem = rel.reduce(rem)
        gens.append(rem)

    basis = grobner(gens)

    size = sage_grobner(gens)
    assert size == len(basis)

    weights = {}
    for p in basis:
        g = find_grade(grades, p)
        weights[g] = weights.get(g, 0) + 1
    return weights


def show_weights_62():
    weights = {(-4, -4): 3, (5, 2): 2, (2, 2): 3, (0, 1):
    3, (1, 0): 3, (1, 3): 3, (-1, 5): 2, (0, 4): 3, (3, 1):
    3, (-3, -2): 3, (5, -1): 2, (3, -2): 2, (1, -3): 2, (-1,
    -4): 2, (-3, -5): 2, (-5, -6): 2, (7, 0): 1, (-7, -7):
    1, (4, -3): 1, (2, -4): 1, (0, -5): 1, (-2, -6): 1, (-4,
    -7): 1, (-6, -8): 1, (6, -2): 1, (-1, -1): 3, (-5, -3):
    2, (-2, 0): 3, (-4, -1): 2, (4, 3): 2, (-1, 2): 3, (-3,
    1): 2, (3, 4): 2, (-2, 3): 2, (2, 5): 2, (1, 6): 2, (0,
    7): 1, (4, 0): 3, (2, -1): 3, (0, -2): 3, (-2, -3): 3,
    (-6, -5): 2, (6, 1): 2, (7, 3): 1, (6, 4): 1, (5, 5):
    1, (4, 6): 1, (3, 7): 1, (2, 8): 1, (8, 2): 1, (-7, -4):
    1, (-6, -2): 1, (-5, 0): 1, (-4, 2): 1, (-3, 4): 1, (-2,
    6): 1, (-8, -6): 1}

    weights = dict((k, str(v)) for (k,v) in weights.items())
    diagram = WeightDiagram.make_A2(weights)
    diagram.plot("sl3_weights_62")


try:
    from huygens import config
    config(text="pdflatex")
    from huygens.front import Canvas, path, color, text, style
except ImportError:
    pass

class WeightDiagram(object):
    def __init__(self, basis, weights):
        self.basis = basis
        self.weights = dict(weights)
        pts = []
        for (i,j) in weights.keys():
            x, y = self.coord(i, j)
            pts.append((x, y))
        self.pts = pts
        if len(pts)>=3:
            self.hull = ConvexHull(pts)
        else:
            self.hull = None

    def coord(self, i, j):
        u, v = self.basis
        x, y = i*u[0] + j*v[0], i*u[1] + j*v[1]
        return x, y

    def in_hull(self, i, j):
        hull = self.hull
        if not hull:
            return i==0 and j==0
        x, y = self.coord(i, j)
        for row in hull.equations:
            x0, y0, a = row
            val = x*x0 + y*y0
            if val >= -a + 0.01:
                return False
        return True

    def mark_weight(self, cvs, i, j, value):
        st_center = [text.halign.boxcenter, text.valign.middle]
        x, y = self.coord(i, j)
        if type(value) is int:
            radius = 0.3
            for count in range(value):
                pth = path.circle(x, y, 0.2*radius)
                if count==0:
                    cvs.fill(pth)
                else:
                    cvs.stroke(pth)
                radius += 0.3
        elif type(value) is str:
            if value:
                cvs.fill(path.circle(x, y, 0.20), [color.rgb(0.8, 0.8, 0.8, 1.0)])
                cvs.text(x, y, value, st_center)
            else:
                cvs.fill(path.circle(x, y, 0.20), [color.rgb(0.8, 0.8, 0.8, 0.5)])
        else:
            assert 0, repr(value)

    def mark_hull(self, cvs):
        hull = self.hull
        pts = self.pts
        weights = self.weights
        lred = color.rgb(1.0, 0.3, 0.2, 0.7)

        fg = Canvas()

        # construct the boundary of the hull
        pairs = [(i,j) for (i, j) in hull.simplices]
        bdy = []
        i,j = pairs.pop()
        bdy.append(i)
        while pairs:
            for pair in pairs:
                if j in pair:
                    break
            else:
                assert 0
            pairs.remove(pair)
            bdy.append(j)
            if pair[0]==j:
                j = pair[1]
            else:
                j = pair[0]
        items = [path.moveto(pts[bdy[0]][0], pts[bdy[0]][1])]
        for i in bdy[1:]:
            items.append(path.lineto(pts[i][0], pts[i][1]))
        items.append(path.closepath())
        p = path.path(items)
        fg.stroke(p, [style.linewidth.Thick])

        fg.clip(p)

        imin, imax = 0, 0
        for (i, j) in weights.keys():
            imin = min(imin, min(i, j))
            imax = max(imax, max(i, j))

        for i in range(imin, imax+1):
          for j in range(imin, imax+1):
            if not self.in_hull(i, j):
                continue
            x, y = self.coord(i, j)
            for row in hull.equations:
                x0, y0, a = row
                dx, dy = 2*a*x0, 2*a*y0
                fg.stroke(path.line(x-dx, y-dy, x+dx, y+dy))
            #for (dx, dy) in self.basis:
            #    fg.stroke(path.line(x-dx, y-dy, x+dx, y+dy))

        #for (dx, dy) in self.basis:
        for row in hull.equations:
            x0, y0, a = row
            dx, dy = 2*a*x0, 2*a*y0
            fg.stroke(path.line(-dx, -dy, +dx, +dy), [lred, style.linewidth.Thick])

        for i in range(imin, imax+1):
          for j in range(imin, imax+1):
            if not self.in_hull(i, j):
                continue
            x, y = self.coord(i, j)
            #self.mark_weight(fg, i, j, "")
            fg.fill(path.circle(x, y, 0.05))
        cvs.append(fg)
    
    def plot(self, filename="output"):
        self.cvs = cvs = Canvas()

        if self.hull is not None:
            self.mark_hull(cvs)

        weights = self.weights
        for (key, value) in weights.items():
            i, j = key
            self.mark_weight(cvs, i, j, value)
        assert ".pdf" not in filename
        cvs.writePDFfile(filename+".pdf")
        cvs.writeSVGfile(filename+".svg")

    @classmethod
    def make_A2(cls, weights, scale=0.8):
        basis = [
            (1.*scale, 0*scale),
            (cos(4/3*pi)*scale, sin(4/3*pi)*scale)]
        return WeightDiagram(basis, weights)

    @classmethod
    def make_B2(cls, weights, scale=0.8):
        basis = [(-1.*scale/2, 1*scale/2), (0*scale, 1*scale)]
        return WeightDiagram(basis, weights)




def test_graded_sl4():
    # See: Miller & Sturmfels, p276

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    poly = lambda v : Poly(v, ring)

    p1 = poly("p1")
    p2 = poly("p2")
    p3 = poly("p3")
    p4 = poly("p4")
    p12 = poly("p12")
    p13 = poly("p13")
    p14 = poly("p14")
    p23 = poly("p23")
    p24 = poly("p24")
    p34 = poly("p34")
    p123 = poly("p123")
    p124 = poly("p124")
    p134 = poly("p134")
    p234 = poly("p234")

    rels = [
        p23*p1 - p13*p2 + p12*p3,       p24*p1 - p14*p2 + p12*p4,
        p34*p1 - p14*p3 + p13*p4,       p34*p2 - p24*p3 + p23*p4,
        p14*p23 - p13*p24 + p12*p34,    p234*p1 - p134*p2 + p124*p3 - p123*p4,
        p134*p12 - p124*p13 + p123*p14, p234*p12 - p124*p23 + p123*p24,
        p234*p13 - p134*p23 + p123*p34, p234*p14 - p134*p24 + p124*p34,
    ]
    rels = grobner(rels, verbose=True)
    print("rels:", rels)
    print()

    grades = [
        [p1, p2, p3, p4],
        [p12, p13, p14, p23, p24, p34],
        [p123, p124, p134, p234],
    ]
    multi = argv.get("multi")
    n = 5 if multi is None else sum(multi)+1
    n = argv.get("n", n)
    for g0 in range(n):
     for g1 in range(n):
      for g2 in range(n):
        if multi is not None and (g0, g1, g2)!=multi:
            #print(".  ", end='')
            continue
        elif g0+g1+g2 > n-1:
            print(".  ", end='')
            continue
        gens = []
        for m0 in all_monomials(grades[0], g0, ring):
         for m1 in all_monomials(grades[1], g1, ring):
          for m2 in all_monomials(grades[2], g2, ring):
            m = m0*m1*m2
            #for rel in rels:
            #    div, m = rel.reduce(m)
            m = reduce_many(rels, m)
            if m != 0:
                gens.append(m)
            
        print(len(gens), end=':', flush=True)
        basis = grobner(gens)
        lhs = len(basis)
        rhs = (g0+1)*(g1+1)*(g2+1)*(g0+g1+2)*(g1+g2+2)*(g0+g1+g2+3)//12
        assert lhs==rhs, ("%s != %s"%(lhs, rhs))
        print(len(basis), end=' ', flush=True)

#        basis.sort(key=str)
#        heads = {}
#        for p in basis:
#            print(p.head, p)
#            heads[p.head] = p
#        print(len(heads))
#        return

      print()
     print()



def test_plucker():
    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)

    rows, cols = argv.get("rows", 2), argv.get("cols", 4)
    U = numpy.empty((rows, cols), dtype=object)
    for i in range(rows):
      for j in range(cols):
        U[i, j] = Poly("x[%d,%d]"%(i, j), ring)

    print(U)
    COLS = list(range(cols))
    w = {} # the plucker coordinates
    for idxs in choose(COLS, rows):
        V = U[:, idxs]
        #print(V)
        a = determinant(V)
        w[idxs] = a
        #print(idxs, a)

    if (rows, cols) == (2, 4):
        assert w[0,1]*w[2,3]-w[0,2]*w[1,3]+w[0,3]*w[1,2] == 0

    for idxs in choose(COLS, rows-1):
      for jdxs in choose(COLS, rows+1):
        if len(idxs) and idxs[-1] >= jdxs[0]:
            continue
        #print(idxs, jdxs)
        sign = ring.one
        rel = ring.zero
        for l in range(rows+1):
            ldxs = idxs+(jdxs[l],)
            rdxs = jdxs[:l] + jdxs[l+1:]
            rel += sign*w[ldxs]*w[rdxs]
            sign *= -1
        assert rel==0


def test_plucker_flag():
    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)

    n = argv.get("n", 4)
    U = numpy.empty((n, n), dtype=object)
    for i in range(n):
      for j in range(n):
        U[i, j] = Poly("x[%d,%d]"%(i, j), ring)

    print(U)

    N = list(range(n))
    w = {} # the plucker coordinates
    for k in range(1, n):
      for idxs in choose(N, k):
        V = U[:k, idxs]
        #print(V)
        a = determinant(V)
        if k==1:
            w[idxs[0]] = a
        else:
            w[idxs] = a
        #print(idxs, a)

    assert n==4
    p1 = w[0]
    p2 = w[1]
    p3 = w[2]
    p4 = w[3]
    p12 = w[0,1]
    p13 = w[0,2]
    p14 = w[0,3]
    p23 = w[1,2]
    p24 = w[1,3]
    p34 = w[2,3]
    p123 = w[0,1,2]
    p124 = w[0,1,3]
    p134 = w[0,2,3]
    p234 = w[1,2,3]

    for rel in [
        p23*p1 - p13*p2 + p12*p3,       p24*p1 - p14*p2 + p12*p4,
        p34*p1 - p14*p3 + p13*p4,       p34*p2 - p24*p3 + p23*p4,
        p14*p23 - p13*p24 + p12*p34,    p234*p1 - p134*p2 + p124*p3 - p123*p4,
        p134*p12 - p124*p13 + p123*p14, p234*p12 - p124*p23 + p123*p24,
        p234*p13 - p134*p23 + p123*p34, p234*p14 - p134*p24 + p124*p34,
    ]:
        assert rel == 0

    return

    for idxs in choose(N, rows-1):
      for jdxs in choose(N, rows+1):
        if len(idxs) and idxs[-1] >= jdxs[0]:
            continue
        print(idxs, jdxs)
        sign = ring.one
        rel = ring.zero
        for l in range(rows+1):
            ldxs = idxs+(jdxs[l],)
            rdxs = jdxs[:l] + jdxs[l+1:]
            rel += sign*w[ldxs]*w[rdxs]
            sign *= -1
        print(rel)


def test_dimension():
    name = argv.next() or "A2"
    N = int(argv.next() or 4)
    data = os.popen("./sl.sage %s %s"%(name, N)).read()
    data = eval(data)
    print(data)
    data = numpy.array(data)

    series = name[0]
    n = int(name[1])

    def Adim(*xs):
        if len(xs)==1:
            x = xs[0]
            return x+1
        #print("Adim", xs)
        n = len(xs)
        value = Adim(*xs[:-1])
        #print("value =", value)
        for i in range(n):
            rhs = sum(xs[i:]) + n-i
            #print("   *= sum(%s) + %s"%(xs[i:], n-i))
            value *= rhs
        div = factorial(n)
        assert value%div == 0
        value //= factorial(n)
        #print("  =", value)
        return value

    def Cdim(*xs):
        assert len(xs)>1
        value = Adim(*xs)
        if len(xs)==2:
            return value * (xs[0] + 2*xs[1] + 3) // 3
        if len(xs)==3:
            a, b, c = xs
            value *= (b+2*c+3)*(a+b+2*c+4)*(a+2*b+2*c+5)
            div = 60
            assert value%div == 0
            value //= div
            return value
        if len(xs)==4:
            a, b, c, d = xs
            value *= (a+b+c+2*d+5)*(a+b+2*c+2*d+6)*(a+2*b+2*c+2*d+7)
            value *= (b+c+2*d+4)*(b+2*c+2*d+5)*(c+2*d+3)
            div = 5*6*7*4*5*3
            assert value%div == 0
            value //= div
            return value

    if name[0] == "A":
        fn = Adim
    elif name[0] == "B":
        fn = Bdim
    elif name[0] == "C":
        fn = Cdim
    else:
        assert 0, name

    idxs = tuple(range(N))
    idxss = tuple(idxs for i in range(n))
    result = numpy.zeros((N,)*n, dtype=int)
    for idx in numpy.ndindex(result.shape):
        value = fn(*idx)
        result[idx] = value
    print(result)
    assert numpy.alltrue(data == result)


def test_quaternion():
    ARing = type("ARing", (Ring,), {})
    ring = ARing()
    assert isinstance(ring, Ring)
    ring.zero = zero
    ring.one = one

    space = Space(ring, 4, name="U")

    def quaternion(a, b, c, d):
        # build matrix representation of quaternion
        A = numpy.empty((4, 4), dtype=object)
        A[:] = [
            [a, -b, -c, -d],
            [b, a, -d, c],
            [c, d, a, -b],
            [d, -c, b, a]]
        return Lin(space, space, A)
    
    e = quaternion(1, 0, 0, 0)
    i = quaternion(0, 1, 0, 0)
    j = quaternion(0, 0, 1, 0)
    k = quaternion(0, 0, 0, 1)

    l = quaternion(a, b, c, d)       # light quaternion
    s = quaternion(zero, x, y, z)    # space quaternion
    t = quaternion(zero, u, v, zero) # time quaternion

    # action of light quaternion by conjugation gives the 
    # isometry from time (2dim) to space (3dim):
    lhs = t*l
    rhs = l*s

    for i in range(4):
      for j in range(4):
        print(lhs[i, j], "=", rhs[i, j])

    result = """
    a*u + d*v = a*x + c*z - d*y
    a*v - d*u = a*y - b*z + d*x
    -b*u - c*v = -b*x - c*y - d*z
    b*v - c*u = -a*z - b*y + c*x
    """
    return


#class ARing(Ring):
#    pass

def test_so5():
    # Plucker embedding for SO(5) ~= SO(2,3)

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)

    a = Poly("a", ring)
    b = Poly("b", ring)
    c = Poly("c", ring)
    d = Poly("d", ring)

    # (2,3) space-time
    u = Poly("u", ring) # time coord
    v = Poly("v", ring) # time coord
    x = Poly("x", ring) # space coord
    y = Poly("y", ring) # space coord
    z = Poly("z", ring) # space coord

    # got from test_quaternion :
    rels = [
        u**2+v**2-x**2-y**2-z**2,
        a*u + d*v -(a*x + c*z - d*y),
        a*v - d*u -(a*y - b*z + d*x),
        -b*u - c*v -(-b*x - c*y - d*z),
        b*v - c*u -(-a*z - b*y + c*x),
    ]
    rels = grobner(rels)

    for idx in range(6):
      for jdx in range(6):
        gens = []
        for p in all_monomials([a, b, c, d], idx, ring):
          for q in all_monomials([u, v, x, y, z], jdx, ring):
            rem = p*q
            #for count in range(3):
            while 1:
                rem0 = rem
                for rel in rels:
                    div, rem = rel.reduce(rem)
                if rem == rem0:
                    break
            gens.append(rem)
    
        #print("gens:", len(gens))
        n = sage_grobner(gens)
        #basis = grobner(gens)
        print("%3d"%n, end=' ', flush=True)
      print()


def get_weights_so5(idx=0, jdx=0):
    # Plucker embedding for SO(5) ~= SO(2+3)

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)

    a = Poly("a", ring)
    b = Poly("b", ring)
    c = Poly("c", ring)
    d = Poly("d", ring)

    z = Poly("z", ring)

    e = Poly("e", ring)
    f = Poly("f", ring)
    g = Poly("g", ring)
    h = Poly("h", ring)

    half = one/2
    u = half*(e+g)
    v = half*(f+h)
    x = half*(e-g)
    y = half*(h-f)

    assert u**2+v**2-x**2-y**2-z**2 == e*g + f*h - z**2


    # got from test_quaternion :
    rels = [
        e*g + f*h - z**2,
        a*u + d*v -(a*x + c*z - d*y),
        a*v - d*u -(a*y - b*z + d*x),
        -b*u - c*v -(-b*x - c*y - d*z),
        b*v - c*u -(-a*z - b*y + c*x),
    ]
    print("rels:")
    for rel in rels:
        print(rel)
    print()
    rels = grobner(rels)

    gens = []
    for p in all_monomials([a, b, c, d], idx, ring):
      for q in all_monomials([e, f, g, h, z], jdx, ring):
        rem = p*q
        #for count in range(3):
        while 1:
            rem0 = rem
            for rel in rels:
                div, rem = rel.reduce(rem)
            if rem == rem0:
                break
        gens.append(rem)

    #print("gens:", len(gens))
    basis = grobner(gens)
    for p in basis:
        print(p)
    n = sage_grobner(gens)
    assert n == len(basis)
    print("%3d"%n)

    grades = {
        #       s    q
        'e' : ( 0,   1),
        'b' : (-1,   1),
        'f' : (-2,   1),
        'a' : ( 1,   0),
        'z' : ( 0,   0),
        'd' : (-1,   0),
        'h' : ( 2,  -1),
        'c' : ( 1,  -1),
        'g' : ( 0,  -1),
    }

    weights = {}
    for p in basis:
        g = find_grade(grades, p)
        weights[g] = weights.get(g, 0) + 1
    return weights


def plot_weights_so5():
    idx = argv.get("s", 0)
    jdx = argv.get("q", 0)
    weights = get_weights_so5(idx, jdx)
    name = argv.get("name", "weights_so5")
    weights = dict((k, str(v)) for (k,v) in weights.items())
    diagram = WeightDiagram.make_B2(weights)
    diagram.plot(name)


def plot_all_so5():
    for idx in range(3):
      for jdx in range(3):
        weights = get_weights_so5(idx, jdx)
        weights = dict((k, str(v)) for (k,v) in weights.items())
        name = "so5_weights_%d%d"%(idx, jdx)
        print(name)
        #plot_weights(weights)
        diagram = WeightDiagram.make_B2(weights)
        diagram.plot(name)


    


if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        seed(_seed)

    profile = argv.profile
    fn = argv.next()

    if profile:
        import cProfile as profile
        profile.run("test()")

    elif fn is None:
        test()

    else:
        fn = eval(fn)
        fn()

    print("OK")

