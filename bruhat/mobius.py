#!/usr/bin/env python3

"""
Mobius transforms and triangle groups.

"""

import sys
from random import shuffle
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


"""

The automorphisms of the disc is SU(1,1), which
lives inside SL(2, C) as:
    [[u v]
     [v* u*]]
    such that u.u* - v.v* == 1.

Refs:
[1] "Indra's Pearls"
[2] https://arxiv.org/pdf/1311.2081v2.pdf
"""

class Mobius(object):
    """
    [[a b]
     [c d]]
    """
    def __init__(self, a=1., b=0., c=0., d=1.):
        #assert abs(det - 1.) < EPSILON, det
        self.a = complex(a)
        self.b = complex(b)
        self.c = complex(c)
        self.d = complex(d)

    @property
    def det(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        return a*d - b*c

    @property
    def trace(self):
        return self.a + self.d

    def __eq__(g, h):
        return abs(g.a-h.a)+abs(g.b-h.b)+ abs(g.c-h.c)+abs(g.d-h.d) < EPSILON
    #def __ne__(g, h):
    #    return not g.__eq__(h)

    def __str__(self):
        flds = (self.a, self.b, self.c, self.d)
        return "[[%s %s]\n [%s %s]]"%flds
    __repr__ = __str__

    def __mul__(g, h):
        """
        [[a b]  * [[a b]
         [c d]]    [c d]]
        """
        a = g.a*h.a + g.b*h.c 
        b = g.a*h.b + g.b*h.d 
        c = g.c*h.a + g.d*h.c 
        d = g.c*h.b + g.d*h.d 
        return Mobius(a, b, c, d)

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
        a, b, c, d = (self.a, self.b, self.c, self.d)
        if abs(c) < EPSILON and abs(b) < EPSILON:
            return (0.j, 0.j)
        elif abs(c) < EPSILON:
            return None # fixes point at infinity
        disc = (a+d)**2 - 4*(a*d - b*c)
        #assert disc >= 0.
        return (a-d+disc**0.5)/(2*c), (a-d-disc**0.5)/(2*c)

    def inner_fixed(self):
        z, w = self.fixed()
        if abs(z) < 1+EPSILON:
            return z
        return w

    def is_translate(self):
        if self == I:
            return True
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
        a, b, c, d = (self.a, self.b, self.c, self.d)
        det = a*d - b*c
        return Mobius(d/det, -b/det, -c/det, a/det)
    __invert__ = inv

    @classmethod
    def rotate(cls, theta):
        a = cmath.exp(1j*theta)
        d = cmath.exp(-1j*theta)
        return Mobius(a, 0., 0., d)

    def __call__(self, z):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        top, bot = (a*z + b), (c*z + d)
        if abs(bot) < EPSILON:
            return None
        w = top / bot
        return w

    def trafo(self, x, y):
        z = x + 1.j*y
        w = self(z)
        assert w is not None
        return w.real, w.imag

    def is_sl(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        return abs(self.det - 1.) < EPSILON

    def is_su(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        # u, v, v*, u*
        return abs(complex(d).conjugate() - a) < EPSILON and \
            abs(complex(c).conjugate() - b) < EPSILON

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
        a, b, c, d = (g.a-h.a, g.b-h.b, g.c-h.c, g.d-h.d)
        return (a*a.conjugate() + b*b.conjugate() 
            + c*c.conjugate() + d*d.conjugate()).real**0.5

    @classmethod
    def get_translate(cls, z, w):
        "find Mobius translation that sends z to w"
        # First send to the upper half-plane:
        K = Mobius.cayley()
        u, v = (~K)(z), (~K)(w)
        # Now find the geodesic that goes through u and v:
        top = v.real**2 + v.imag**2 - u.real**2 - u.imag**2
        bot = 2*(v.real - u.real)
        ax = top / bot
        radius = ((u.real - ax)**2 + u.imag**2)**0.5
        #print("ax:", ax, "radius:", radius)
        # The left and right endpoints of the geodesic:
        left, right = ax-radius, ax+radius
        #disc.show_point(K(left))
        #disc.show_point(K(right))
        # Now send this geodesic to the +ve imaginary axis
        L = Mobius(0., -1., 1., -right)
        L = Mobius(1., -L(left), 0., 1.) * L
        # L will send left -> 0, right -> inf
        assert abs(L(left)) < EPSILON
        # R will send L(u) -> L(v)
        R = Mobius(L(v).imag / L(u).imag, 0., 0., 1.)
        # The final translation:
        g = K * (~L) * R * L * (~K)
        assert abs(g(z)-w) < EPSILON
        return g

SL = lambda a, b, c : Mobius(a, b, c, (b*c+1)/a)

I = Mobius()

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


class Cayley(object):
    def __init__(self, table):
        self.table = dict(table)

    def __len__(self):
        return len(self.table)

    def dump(self):
        table = self.table
        vals = list(set(table.values()))
        vals.sort()
        print('     '+''.join("%3d"%i for i in vals))
        print('     '+'-'*3*len(vals))
        for i in vals:
          print("%3d |"%i, end="")
          for j in vals:
            k = table.get((i, j))
            if k is None:
                s = "   "
            else:
                s = "%3d"%k
            print(s, end="")
          print()

    def _deduce(self):
        table = self.table
        vals = list(set(table.values()))
        vals.sort()
        for i in vals:
            found = {}
            for j in vals:
                k = table.get((i, j))
                if k is None:
                    continue
                if k not in found:
                    found[k] = j
                    continue
                j1 = found[k]
                # i*j == i*j1 == k => j==j1
                self.rewrite(j1, j)
                return True
        return False

    def deduce(self):
        while self._deduce():
            pass

    def rewrite(self, i, j):
        # send j to i
        table = self.table
        keys = list(table.keys())
        for ii, jj in keys:
            kk = table[ii, jj]
            changed = False
            if ii == j or jj == j:
                del table[ii, jj]
                if ii == j:
                    ii = i
                if jj == j:
                    jj = i
                changed = True
            if kk == j:
                kk = i
                changed = True
            if changed:
                table[ii, jj] = kk
    



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


    ops = list(G)
    for idx, op in enumerate(ops):
        op.idx = idx
    table = {}
    for g in ops:
      for h in ops:
        k = g*h
        if k not in G:
            continue
        k = G.get(k)
        table[g.idx, h.idx] = k.idx

    table = Cayley(table)
    i = 0
    while i < len(ops):
        g = ops[i]
        j = i+1
        while j < len(ops):
            h = ops[j]
            if (~g)*h in H:
                table.rewrite(g.idx, h.idx)
                ops.pop(j)
            else:
                j += 1
        i += 1
        print(len(table), len(ops))

    table.dump()
    table.deduce()
    print()
    table.dump()

    return

    done = False
    while not done:
        done = True
        i = 0
        while i < len(ops):
            g = ops[i]
            j = i+1
            while j < len(ops):
                h = ops[j]
                if (~g)*h in H:
                    ops.pop(j)
                    done = False
                else:
                    j += 1
            i += 1
        print(done)
    print(len(ops))
    return

    #for g in ops:
    #  for h in ops:
    #    if g==h:
    #        continue
    #    print(int(relation(g, h)), end="")
    #  print()

    equs = [equ.Equ(op) for op in ops]
    n = len(ops)
    for i in range(n):
        print(len(set(e.top for e in equs)), end=" ", flush=True)
        ei = equs[i]
        for j in range(i+1, n):
            ej = equs[j]
            if ei.eq(ej):
                continue
            if relation(ops[i], ops[j]):
                ei.merge(ej)
    #print(len([e for e in equs if e.top is e]))
    print(len(set(e.top for e in equs)))
    
    return

    K = Generator([I])
    for g in G:
      for h in H:
        gh = g*h
        if gh in K:
            break
      else:
        K.add(g)
    print(len(K))

    return

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





