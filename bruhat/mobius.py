#!/usr/bin/env python3

"""
Mobius transforms, reflections, and triangle groups.

The automorphisms of the disc is SU(1,1), which
lives inside SL(2, C) as:
    [[u v]
     [v* u*]]
    such that u.u* - v.v* == 1.

Refs:
[1] "Indra's Pearls"
[2] https://arxiv.org/pdf/1311.2081v2.pdf
[3] https://www.math.cuhk.edu.hk/course_builder/2021/math4900e/MATH4900E_group3%20Oct%2012.pdf
"""


import sys
from random import shuffle, random
from math import pi, hypot, cos, sin, tan, acosh, atan
import cmath

import kdtree

from bruhat import equ
from bruhat.argv import argv

EPSILON = 1e-6


def euclid(p, q):
    return ((p[0]-q[0])**2 + (p[1]-q[1])**2)**0.5


def poincare(z, w):
    "poincare disc metric"
    z, w = z*1j, w*1j
    top = 2*(z-w)*(z-w).conjugate()
    bot = (1 - x*x.conjugate()) * (1 - y*y.conjugate())
    d = acosh(1 + (top/bot).real)
    return d

#d_poincare = lambda z : 2 / (1 - z*z.conjugate()).real
def d_poincare(z):
    z = complex(z)
    zz = z*z.conjugate()
    assert abs(zz) < 1
    ds = 2 / (1 - zz).real
    return ds


class Mobius(object):
    """ 
    Warning: these matrices are in GL(2) not PGL(2) so they
    act projectively on the Disc, etc.
    [[a b]
     [c d]]
    """
    def __init__(self, a=1., b=0., c=0., d=1., conjugate=False):
        #assert abs(det - 1.) < EPSILON, det
        self.a = complex(a)
        self.b = complex(b)
        self.c = complex(c)
        self.d = complex(d)
        self.conjugate = conjugate

    @property
    def det(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        return a*d - b*c

    @property
    def trace(self):
        return self.a + self.d

    def __eq__(g, h):
        return g.conjugate == h.conjugate and \
            abs(g.a-h.a)+abs(g.b-h.b)+ abs(g.c-h.c)+abs(g.d-h.d) < EPSILON

    def __str__(self):
        flds = (self.a, self.b, self.c, self.d)
        s = "[[%s %s]\n [%s %s]]"%flds
        if self.conjugate:
            s = "*"+s
        return s
    __repr__ = __str__

    def __mul__(g, h):
        """
        [[a b]  * [[a b]
         [c d]]    [c d]]
        """
        # See ref. [3], page 45
        conjugate = (g.conjugate != h.conjugate)
        if g.conjugate == False:
            a = g.a*h.a + g.b*h.c 
            b = g.a*h.b + g.b*h.d 
            c = g.c*h.a + g.d*h.c 
            d = g.c*h.b + g.d*h.d 
        else:
            a = g.a*h.a.conjugate() + g.b*h.c.conjugate() 
            b = g.a*h.b.conjugate() + g.b*h.d.conjugate() 
            c = g.c*h.a.conjugate() + g.d*h.c.conjugate() 
            d = g.c*h.b.conjugate() + g.d*h.d.conjugate() 
        return Mobius(a, b, c, d, conjugate)

    def __getitem__(self, idx):
        #return (self.a, self.b, self.c, self.d)[idx] # used by kdtree
        r = (
            self.a.real, self.b.real, self.c.real, self.d.real,
            self.a.imag, self.b.imag, self.c.imag, self.d.imag,
        )[idx] # used by kdtree
        assert isinstance(r, float)
        return r

    REAL_DIMS = 4
    def __len__(self):
        return self.REAL_DIMS # used by kdtree

    def fixed(self):
        assert self.conjugate == False, "not implemented"
        a, b, c, d = (self.a, self.b, self.c, self.d)
        if abs(c) < EPSILON and abs(b) < EPSILON:
            return (0.j, 0.j)
        elif abs(c) < EPSILON:
            return None # fixes point at infinity
        disc = (a+d)**2 - 4*(a*d - b*c)
        #assert disc >= 0.
        return (a-d+disc**0.5)/(2*c), (a-d-disc**0.5)/(2*c)

    def inner_fixed(self):
        assert self.conjugate == False, "not implemented"
        z, w = self.fixed()
        if abs(z) < 1+EPSILON:
            return z
        return w

    def is_translate(self):
        if self == I:
            return True
        if self.conjugate:
            return False
        z, w = self.fixed()
        if abs(z-w) < EPSILON:
            return False
        if abs(z*z.conjugate() - 1.) > EPSILON:
            return False
        if abs(w*w.conjugate() - 1.) > EPSILON:
            return False
        return True

    @classmethod
    def cayley(cls, z=1.j):
        "send upper half-plane to unit disc, and z to zero"
        z = complex(z)
        assert z.imag > EPSILON, "z must be in upper half-plane"
        return cls(1., -z, 1., -z.conjugate())

    @classmethod
    def disc_tozero(cls, z):
        K = ~cls.cayley()
        K1 = cls.cayley(K(z))
        g = K1*K
        assert abs(g(z)) < EPSILON
        return g

    def todisc(g, z=1.j):
        # send g in SL(2,R) to SU(1,1) automorphism of the disc
        K = Mobius.cayley(z)
        return K * g * ~K

    def inv(self):
        #assert self.conjugate == False, "not implemented"
        a, b, c, d = (self.a, self.b, self.c, self.d)
        det = a*d - b*c
        return Mobius(d/det, -b/det, -c/det, a/det, self.conjugate)
    __invert__ = inv

    @classmethod
    def rotate(cls, theta):
        a = cmath.exp(1j*theta)
        d = cmath.exp(-1j*theta)
        return Mobius(a, 0., 0., d)

    @classmethod
    def conjugate(cls):
        return Mobius(1, 0, 0, 1, True)

    def __call__(self, z):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        if z is None: # infinity
            top, bot = a, c
            if abs(bot) < EPSILON:
                return None # infinity
            return top / bot
        if self.conjugate:
            z = z.conjugate()
        top, bot = (a*z + b), (c*z + d)
        if abs(bot) < EPSILON:
            return None # infinity
        w = top / bot
        return w

    def trafo(self, x, y):
        z = x + 1.j*y
        w = self(z)
        assert w is not None
        return w.real, w.imag

    def is_sl(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        return abs(self.det - 1.) < EPSILON and not self.conjugate

    def is_su(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        # u, v, v*, u*
        return abs(complex(d).conjugate() - a) < EPSILON and \
            abs(complex(c).conjugate() - b) < EPSILON and not self.conjugate

    @classmethod
    def _send(cls, p, q, r):
        # [1] pg 93
        top, bot = q-r, q-p
        g = Mobius(top, -p*top, bot, -r*bot)
        return g

    @classmethod
    def send(cls, p0, q0, r0, p1, q1, r1):
        "Mobius transforms are triply transitive"
        g = cls._send(p0, q0, r0)
        h = cls._send(p1, q1, r1)
        return (~h)*g

    def __pow__(g, e):
        if e < 0:
            return (~g).__pow__(-e) # recurse
        op = Mobius()
        for i in range(e):
            op = g*op
        return op

    def order(g, maxorder=999):
        I = Mobius()
        h = g
        count = 1
        while h != I:
            h = g*h
            assert h != g
            count += 1
            if count >= maxorder:
                #assert 0, "infinite order?"
                return None
        return count

    def dist(g, h):
        assert g.conjugate == h.conjugate == False, "not implemented"
        a, b, c, d = (g.a-h.a, g.b-h.b, g.c-h.c, g.d-h.d)
        return (a*a.conjugate() + b*b.conjugate() 
            + c*c.conjugate() + d*d.conjugate()).real**0.5

    @classmethod
    def get_translate(cls, z, w):
        "find Mobius translation that sends z to w"
        # First send to the upper half-plane:
        #print("get_translate", z, "-->", w)
        K = Mobius.cayley()
        u, v = (~K)(z), (~K)(w)
        #print("\t", u, "-->", v)
        # Now find the geodesic that goes through u and v:
        if abs(v.real - u.real) < 1e-8:
            # Now send this geodesic to the +ve imaginary axis
            x = v.real
            L = Mobius(1, -x, 0, 1)
            assert abs(L(x)) < EPSILON
        else:
            bot = 2*(v.real - u.real)
            top = v.real**2 + v.imag**2 - u.real**2 - u.imag**2
            #print("%s / %s" % (top, bot))
            ax = top / bot
            radius = ((u.real - ax)**2 + u.imag**2)**0.5
            #print("ax:", ax, "radius:", radius)
            # The left and right endpoints of the geodesic:
            left, right = ax-radius, ax+radius
            # Now send this geodesic to the +ve imaginary axis
            L = Mobius(0., -1., 1., -right)
            L = Mobius(1., -L(left), 0., 1.) * L
            # L will send left -> 0, right -> inf
            assert abs(L(left)) < EPSILON
        # R will send L(u) -> L(v)
        Lv, Lu = L(v), L(u)
        assert Lv is not None, "v = %s"%v
        assert Lu is not None, "u = %s"%u
        #print("\t", Lu, Lv)
        R = Mobius(Lv.imag / Lu.imag, 0., 0., 1.)
        #print("\t\t", R(Lu))
        # The final translation:
        g = K * (~L) * R * L * (~K)
        assert abs(g(z)-w) < EPSILON, "g(z)=%s"%(g(z),)
        return g

SL = lambda a, b, c : Mobius(a, b, c, (b*c+1)/a)

I = Mobius()


# does not need hashable operators
def mulclose(gen, g0=I, verbose=False, maxsize=None):
    ops = list(gen)
    bdy = gen
    while bdy:
        _bdy = []
        for g in bdy:
            for h in gen:
                k = g*h
                if k not in ops:
                    ops.append(k)
                    _bdy.append(k)
            #if maxsize and len(ops) >= maxsize:
            #    break
            #_bdy.sort(key = lambda g : g.dist(g0))
        bdy = _bdy
        if verbose:
            print("mulclose:", len(ops))
        if maxsize and len(ops) >= maxsize:
            break

    return ops




def mktriangle(a, b, c):
    # ref [2] Proposition 2.5
    assert a>1 and b>1 and c>1 
    ca, sa = cos(pi/a), sin(pi/a)
    cb, sb = cos(pi/b), sin(pi/b)
    cc, sc = cos(pi/c), sin(pi/c)
    ga = Mobius(ca, sa, -sa, ca)
    assert ga.is_sl()
    lmbda = (ca * cb + cc) / (sa*sb) # factor of 2 in denominator ?
    #print("lmbda =", lmbda)
    assert lmbda >= 1., lmbda
    mu = lmbda + (lmbda**2-1)**0.5
    #print("1/mu =", 1/mu)
    gb = Mobius(cb, mu*sb, -sb/mu, cb)
    #print(gb)
    assert gb.is_sl()
    return [ga, gb]


class Generator(object):
    def __init__(self, gen, verbose=False, maxsize=None):
        tree = kdtree.create(dimensions=Mobius.REAL_DIMS)
        self.tree = tree
        self.size = 0
        for op in gen:
            tree.add(op)
            self.size += 1
        if I not in gen:
            tree.add(I)
            self.size += 1
        self.gen = list(gen)
        if maxsize is not None:
            self.generate(verbose, maxsize)

    def __len__(self):
        return self.size

    def add(self, op):
        if op not in self:
            self.tree.add(op)
            self.size += 1

    def generate(self, verbose=False, maxsize=None):
        tree, size = self.tree, self.size
        bdy = gen = self.gen
        while bdy:
            _bdy = []
            for g in bdy:
                for h in gen:
                    k = g*h
                    near, dist = tree.search_nn(k)
                    #print(near.data)
                    if dist > EPSILON:
                        _bdy.append(k)
                        tree.add(k)
                        size += 1
                if maxsize and size >= maxsize:
                    break
            if maxsize and size >= maxsize:
                break
            bdy = _bdy
            if verbose:
                print("[%s]" % size, end="", flush=True)
        if verbose:
            print("[%s]" % size)
        tree = tree.rebalance()
        self.tree = tree
        self.size = size
    
    def getops(self):
        ops = [node.data for node in self.tree.inorder()]
        return ops

    def get(self, op):
        near, dist = self.tree.search_nn(op)
        if dist < EPSILON:
            assert near.data == op
            return near.data
        return None

    def __contains__(self, op):
        near, dist = self.tree.search_nn(op)
        if dist < EPSILON:
            assert near.data == op
            return True
        return False

    def __iter__(self):
        for node in self.tree.inorder():
            yield node.data




def test():
    i = 1j
    I = Mobius()
    g = Mobius(1., 2., 3., 7)
    assert g*~g == I
    K = Mobius.cayley()
    assert K*~K == I
    h = g.todisc()
    assert g.is_sl()
    #print(h)
    for z in (h(1), h(-1), h(i), h(-i)):
        assert abs(z*z.conjugate() - 1.) < EPSILON
    assert h.is_su()

    assert abs(h.det - 1.) < EPSILON, h.det

    theta = 2*pi / 5
    g = Mobius.rotate(theta)
    assert g.order() == 5

    p0, q0, r0 = 1., 2+i, 3+1.1*i
    p1, q1, r1 = 5., 1.1*i, 7*i
    g = Mobius.send(p0, q0, r0, p1, q1, r1)
    assert abs(g(p0) - p1) < EPSILON
    assert abs(g(q0) - q1) < EPSILON
    assert abs(g(q0) - q1) < EPSILON

    rnd = lambda : 2*random() - 1.
    z = rnd() + rnd()*1j
    C = Mobius(1, 0, 0, 1, True) # conjugation
    assert abs(C(z) - z.conjugate()) < EPSILON

    # test group action property
    gs = [Mobius(rnd(), rnd(), rnd(), rnd(), bool(random()>0.5)) for i in range(10)]
    for g in gs:
      for h in gs:
        z = rnd() + rnd()*1j
        lhs = (g*h)(z)
        rhs = g(h(z))
        assert abs(lhs - rhs) < EPSILON

    for g in gs:
        h = ~g
        assert g*h == I
        z = rnd() + rnd()*1j
        lhs = (g*h)(z)
        rhs = g(h(z))
        assert abs(lhs - rhs) < EPSILON
        assert abs(lhs - z) < EPSILON

    N = argv.get("N", 100)
    gens = [Mobius(1, 1, 0, 1), Mobius(0, 1, -1, 0)]
    G = Generator(gens, verbose=False, maxsize=N)

    #gens = [Mobius(1, 2, 0, 1), Mobius(1, 0, -2, 1)]
    #H = Generator(gens, verbose=False, maxsize=N)

    ops = []
    for g in G:
        if g.a.real%3 == 1 and g.c.real%3 == 0 and g.d.real%3 == 1:
            ops.append(g)
    print(len(ops))
    H = Generator(ops, maxsize=N)

    if 0:
        # ------------ {5,5} ------------------------
    
        N = argv.get("N", 100)
        gens = mktriangle(5, 2, 5)
        G = Generator(gens, verbose=True, maxsize=N)
    
        a, b = gens
        assert a.order() == 10
        assert b.order() == 4
        assert (a*b).order() == 5
    
        h = a*a*b
        h = h*h
        assert h.order() == None
    
        H = Generator([h])
        while len(H) < N:
            for g in G:
                h1 = g * h * (~g)
                H.add(h1)
            ops = list(H)
            for h0 in ops:
              for h1 in ops:
                H.add(h0*h1)
            print("|H| =", len(H))


    # ------------ todisc ------------------------

    gens = [g.todisc() for g in gens]
    a, b = gens
    print(a.inner_fixed())
    print(b.inner_fixed())
    print((a*b).inner_fixed())
    #z1 = (a*b).inner_fixed() # pentagon center



if __name__ == "__main__":
    test()

    print("OK\n")





