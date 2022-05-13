#!/usr/bin/env python

"""
some quadratic number rings...

Ref [1]
https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/quadratic-integers-and-coxeter-groups/CF262D475903A0104145D1294DA80EF9
Quadratic Integers and Coxeter Groups
1999 Norman W. Johnson and Asia Ivić Weiss

"""


from functools import total_ordering


#from bruhat.action import Perm, Group, burnside
from bruhat.gset import Perm, Group
from bruhat.element import GenericElement, Ring
from bruhat.lin import Space, Lin
from argv import argv




class NumberRing(Ring):
    def __init__(self, reduct=(-1, 0)):
        self.zero = Number((0, 0), self)
        self.one = Number((1, 0), self)
        self.i = Number((0, 1), self)
        self.reduct = reduct

    def promote(self, a):
        if isinstance(a, Number):
            assert a.tp == self
            return a
        if isinstance(a, int):
            return Number((a, 0), self)
        try:
            a = int(a)
            return Number((a, 0), self)
        except ValueError:
            pass
        except TypeError:
            pass
        return None

    def __truediv__(self, other):
        assert isinstance(other, Ideal)
        return QuotientRing(self, other)

    def norm(self, a):
        a = a.value
        return a[0]**2 + a[1]**2

    def add(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]+b[0], a[1]+b[1]), self)

    def sub(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]-b[0], a[1]-b[1]), self)

    def mul(self, a, b):
        a, b = a.value, b.value
        reduct = self.reduct
        c = a[1]*b[1]
        re = a[0]*b[0] + reduct[0]*c
        im = a[0]*b[1] + a[1]*b[0] + reduct[1]*c
        return Number((re, im), self)

    def neg(self, a):
        a = a.value
        return Number((-a[0], -a[1]), self)

    def mod(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]%b[0], a[1]%b[1]), self)

    def conjugate(self, element):
        assert element.tp == self
        a, b = element.value
        return Number((a, -b), self)

    def i_find(self):
        yield self.zero
        one, i = self.one, self.i
        n = 1
        while 1:
            x = n*one
            while x.value[0] > 0:
                yield x
                x += i-1
            while x.value[1] > 0:
                yield x
                x += -i-1
            while x.value[0] < 0:
                yield x
                x += -i+1
            while x.value[1] < 0:
                yield x
                x += i+1
            n += 1

    def find(self, count):
        items = self.i_find()
        for i in range(count):
            yield items.__next__()


@total_ordering
class Number(GenericElement):

    def __str__(self):
        a, b = self.value
        if b==0:
            s = str(a)
        elif a==0 and b==1:
            s = "i"
        elif a==0 and b==-1:
            s = "-i"
        elif a==0:
            s = "%d*i"%b
        else:
            s = "(%d+%d*i)"%(a, b)
            s = s.replace("+-", "-")
        return s

    # choose some total_ordering ...
    def __lt__(a, b):
        assert isinstance(b, Number)
        assert a.tp == b.tp
        tp = a.tp
        ra = tp.norm(a)
        rb = tp.norm(b)
        if ra < rb:
            return True
        if ra > rb:
            return False
        # now ra==rb
        a, b = a.value, b.value
        if a[0] > b[0]:
            return True
        if a[0] < b[0]:
            return False
        # now a[0] == b[0]
        if a[1] > b[1]:
            return True
        return False


class Ideal(object):
    "_principle ideal"
    def __init__(self, ring, a):
        assert isinstance(ring, NumberRing)
        a = ring.promote(a)
        assert a is not None
        assert a.tp == ring
        self.ring = ring
        #assert a.value[0]>0 and a.value[1]>0
        self.a = a
        self.cache = {}
        i = self.i_find()
        items = []
        while len(items) < 10:
            x = i.__next__()
            if x != 0:
                items.append(x)
        self.items = items

    def __contains__(self, b):
        pass

    def i_find(self):
        ring = self.ring
        a = self.a
        for x in ring.i_find():
            yield x*a

    def reduce(self, src):
        ring = self.ring
        cache = self.cache
        if src in cache:
            return cache[src]
        items = self.items
        done = False
        tgt = src
        while not done:
            done = True
            for x in items:
                if tgt + x < tgt:
                    tgt = tgt+x
                    done = False
        cache[src] = tgt
        return tgt


class QuotientRing(Ring):
    def __init__(self, ring, ideal):
        self.ring = ring
        self.ideal = ideal
        self.one = GenericElement(ring.one, self) # wrapped
        self.zero = GenericElement(ring.zero, self) # wrapped
        self.i = GenericElement(ring.i, self) # wrapped

        # let's hope it's finite and not too big.
        self.items = self._get_items()

    def _get_items(self, N=5):
        one, i = self.one, self.i
        items = set()
        for i in range(-N, N):
          for j in range(-N, N):
            items.add(i*one + j*i)
        items = list(items)
        items.sort(key=str)
        return items

    def __getitem__(self, i):
        return self.items[i]

    def promote(self, a):
        ring = self.ring
        if isinstance(a, GenericElement):
            assert a.tp == self
            return a
        if isinstance(a, Number):
            assert a.tp == ring
            return GenericElement(a, self)
        if isinstance(a, int):
            a = ring.promote(a)
            return GenericElement(a, self)
        try:
            a = int(a)
            a = ring.promote(a)
            return GenericElement(a, self)
        except ValueError:
            pass
        except TypeError:
            pass
        return None

    def add(self, a, b):
        ring, ideal = self.ring, self.ideal
        value = ring.add(a.value, b.value)
        value = ideal.reduce(value)
        a = GenericElement(value, self)
        return a

    def sub(self, a, b):
        ring, ideal = self.ring, self.ideal
        value = ring.sub(a.value, b.value)
        value = ideal.reduce(value)
        a = GenericElement(value, self)
        return a

    def mul(self, a, b):
        ring, ideal = self.ring, self.ideal
        value = ring.mul(a.value, b.value)
        value = ideal.reduce(value)
        a = GenericElement(value, self)
        return a

    def conjugate(self, element):
        assert 0, "don't do this...!"
        ring, ideal = self.ring, self.ideal
        value = ring.conjugate(element.value)
        value = ideal.reduce(value)
        element = GenericElement(value, self)
        return element

    def get_gl(self, V):
        assert V.n == 2, "todo.."
        lins = []
        for a in self:
         for b in self:
          for c in self:
           for d in self:
               if a*d-b*c != 0:
                   f = Lin(V, V, [[a, b], [c, d]])
                   lins.append(f)
        return lins

    def get_sl(self, V):
        assert V.n == 2, "todo.."
        lins = []
        for a in self:
         for b in self:
          for c in self:
           for d in self:
               if a*d-b*c == 1:
                   f = Lin(V, V, [[a, b], [c, d]])
                   lins.append(f)
        return lins




class Mobius(object):
    """ 
    Mobius transform (& conjugate's) over some ring...
    """
    def __init__(self, ring, a=1, b=0, c=0, d=1, conjugate=False):
        assert isinstance(ring, Ring)
        self.ring = ring
        self.a = ring.promote(a)
        self.b = ring.promote(b)
        self.c = ring.promote(c)
        self.d = ring.promote(d)
        assert self.a is not None, repr(a)
        assert self.b is not None, repr(b)
        assert self.c is not None, repr(c)
        assert self.d is not None, repr(d)
        self.conjugate = conjugate
        self._key = None

    @property
    def det(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        return a*d - b*c

    @property
    def trace(self):
        return self.a + self.d

    def key(self):
        if self._key:
            return self._key
        a, b, c, d = (self.a, self.b, self.c, self.d)
        key = []
        i = self.ring.i
        for r in [1, i, -1, -i]:
            item = str((r*a, r*b, r*c, r*d))
            key.append(item)
        key.sort()
        key.append(self.conjugate)
        key = tuple(key)
        self._key = key
        return key

    def __hash__(self):
        return hash(self.key())

    def __eq__(g, h):
        return g.key() == h.key()
        if g.conjugate != h.conjugate:
            return False
        if g.a==h.a and g.b==h.b and g.c==h.c and g.d==h.d:
            return True
        i = g.ring.i
        for r in [-1, i, -i]:
            if g.a==r*h.a and g.b==r*h.b and g.c==r*h.c and g.d==r*h.d:
                return True
        # XXX do this properly...
        return False
        assert 0, "todo"

    def __str__(self):
        flds = (self.a, self.b, self.c, self.d)
        s = "<[[%s %s]\n  [%s %s]]>"%flds
        if self.conjugate:
            s = "*"+s
        return s
    __repr__ = __str__

    def __mul__(g, h):
        """
        [[a b]  * [[a b]
         [c d]]    [c d]]
        """
        # See ref. [1], eq (28)
        assert g.ring == h.ring
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
        return Mobius(g.ring, a, b, c, d, conjugate)

    def __getitem__(self, idx):
        return (self.a, self.b, self.c, self.d)[idx]

    def inv(self):
        a, b, c, d = (self.a, self.b, self.c, self.d)
        det = a*d - b*c
        return Mobius(self.ring, d/det, -b/det, -c/det, a/det, self.conjugate)
    __invert__ = inv

    @classmethod
    def conjugate(cls):
        return Mobius(self.ring, 1, 0, 0, 1, True)

#    def __call__(self, top, bot):
#        a, b, c, d = (self.a, self.b, self.c, self.d)
#        if z is None: # infinity
#            top, bot = a, c
#            if bot == 0:
#                return None # infinity
#            return top / bot
#        z = self.ring.promote(z)
#        if self.conjugate:
#            z = z.conjugate()
#        top, bot = (a*z + b), (c*z + d)
#        if abs(bot) < EPSILON:
#            return None # infinity
#        w = top / bot
#        return w

#    @classmethod
#    def _send(cls, p, q, r):
#        # [1] pg 93
#        top, bot = q-r, q-p
#        g = Mobius(top, -p*top, bot, -r*bot)
#        return g
#
#    @classmethod
#    def send(cls, p0, q0, r0, p1, q1, r1):
#        "Mobius transforms are triply transitive"
#        g = cls._send(p0, q0, r0)
#        h = cls._send(p1, q1, r1)
#        return (~h)*g

    def __pow__(g, e):
        if e < 0:
            return (~g).__pow__(-e) # recurse
        op = Mobius()
        for i in range(e):
            op = g*op
        return op

    def order(g, maxorder=999):
        I = Mobius(g.ring)
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

# does not need hashable operators
def mulclose(gen, verbose=False, maxsize=None):
    ops = list(gen)
    bdy = gen
    while bdy:
        _bdy = []
        for g in bdy:
            for h in gen:
                k = g*h
                for j in ops:
                    if j==k:
                        break
                else:
                    ops.append(k)
                    _bdy.append(k)
        bdy = _bdy
        if verbose:
            print("mulclose:", len(ops))
        if maxsize and len(ops) >= maxsize:
            break
    return ops


def cayley(elements):
    "build the regular permutation representation"
    elements = list(elements)
    lookup = dict((g, idx) for (idx, g) in enumerate(elements))
    assert len(lookup)==len(elements)
    perms = []
    for e in elements:
        perm = []
        for f in elements:
            g = e*f
            k = lookup[g]
            perm.append(k)
        assert len(set(perm)) == len(perm)
        perm = Perm(perm)
        perms.append(perm)
    G = Group(perms)
    hom = dict(zip(elements, G))
    #for g in elements:
    #  for h in elements:
    #    gh = g*h
    #    if hom[gh] == hom[g]*hom[h]:
    #        continue
    #    print(g)
    #    print(h)
    #    print(gh)
    #    assert 0
    return G, hom




def main():

    # ----------------------------
    # Gaussian integers

    R = NumberRing((-1, 0))
    one = R.one
    i = R.i

    assert i*i == -one
    assert i**3 == -i

    assert (1+i)**2 == 2*i
    assert (2+i)*(2-i) == 5
    assert (3+2*i)*(3-2*i) == 13

    items = []
    i = R.i_find()
    while len(items) < 20:
        items.append(i.__next__())
    items.sort()
    #print(' '.join(str(x) for x in items))

    # reflection group
    i = R.i
    I = Mobius(R)
    a = Mobius(R, 1, 0, 0, 1, True)
    b = Mobius(R, i, 0, 0, 1, True)
    c = Mobius(R, -1, 0, 1, 1, True)
    d = Mobius(R, 0, 1, 1, 0, True)
    gens = [a, b, c, d]

    if 0:
        # yes, works...
        G = mulclose(gens, maxsize=100)
        for f in G:
         for g in G:
          for h in G:
            assert f*(g*h) == (f*g)*h
        return

    # Coxeter reflection group: a--4--b--4--c--3--d
    for g in gens:
        assert g.order() == 2
    assert (a*b).order() == 4
    assert (a*c).order() == 2
    assert (a*d).order() == 2
    assert (b*c).order() == 4
    assert (b*d).order() == 2
    assert (c*d).order() == 3

    # rotation group
    a = Mobius(R, -i, 0, 0, 1)
    b = Mobius(R, -i, 0, i, 1)
    c = Mobius(R, 1, 1, -1, 0)

    assert a.order() == 4
    assert b.order() == 4
    assert c.order() == 3

    # ----------------------------
    # Eisenstein integers
    R = NumberRing((-1, -1))
    one = R.one
    i = R.i

    assert i**3 == 1
    assert 1 + i + i**2 == 0

    # ----------------------------
    # Kleinian integers
    R = NumberRing((-2, -1))
    one = R.one
    i = R.i

    assert i * (-1-i) == 2
    assert i * (1+i) == -2

    # ----------------------------

    R = NumberRing((-1, 0))
    I = Ideal(R, 2-R.i)
    S = R/I

    zero, one, i = S.zero, S.one, S.i
    
    assert 1+i != 0
    assert 2-i == 0
    assert 2+i != 0
    assert 1-i == 2+i

    V = Space(S, 2)
    SL = S.get_sl(V)
    GL = S.get_gl(V)

    def lin(a, b, c, d):
        a = S.promote(a)
        b = S.promote(b)
        c = S.promote(c)
        d = S.promote(d)
        return Lin(V, V, [[a, b], [c, d]])
        
    # generators for SL
    gen = lin(0, 1, -1, 0), lin(1, 0, 1, 1), lin(1, 0, i, 1)
    G = mulclose(gen)
    assert len(G) == 120

    i = S.i
    I = Mobius(S)

    if 0:
        # This just does not work over S
        a = Mobius(S, 1, 0, 0, 1, True)
        b = Mobius(S, i, 0, 0, 1, True)
        c = Mobius(S, -1, 0, 1, 1, True)
        d = Mobius(S, 0, 1, 1, 0, True)
        gens = [a, b, c, d]
        G = mulclose(gens)
        inv = {}
        for g in G:
          for h in G:
            if g*h == I:
                inv[g] = h
                inv[h] = g
            assert (g==h) == (hash(g)==hash(h))
            for f in G:
                print(g)
                print(h)
                print(f)
                assert (g*h)*f == g*(h*f)
                print()
        assert len(G) == 240
        return
    
        G, hom = cayley(G)
        assert len(G) == 240
        G.do_check()
    
        return

        # rotation group
        gens = [a*b, b*c, c*d]


def test_gaussian():
    # ----------------------------

    from bruhat.action import mulclose

    R = NumberRing((-1, 0))
    z = 2+3*R.i # |G| = 2184
    z = 3*R.one # |G| = 720
    z = 2+1*R.i # |G| = 120
    #z = -1+R.i # fail
    I = Ideal(R, z)
    S = R/I

    zero, one, i = S.zero, S.one, S.i
    
    a = Mobius(S, -i, 0, 0, 1)
    b = Mobius(S, -i, 0, i, 1)
    c = Mobius(S, 1, 1, -1, 0)

    gens = [a, b, c]
    assert a.order() == 4
    assert b.order() == 4
    assert c.order() == 3

    #print((a*b).order()) # == 2
    #print((b*a).order()) # == 2
    
    G = mulclose(gens)
    #assert len(G) == 120, len(G)
    print(len(G))
    #a, b, c = G[:3]
    #assert a.order() == 4
    #assert b.order() == 4
    #assert c.order() == 3

    G, hom = cayley(G)
    a, b, c = hom[a], hom[b], hom[c]
    assert a.order() == 4
    assert b.order() == 4
    assert c.order() == 3

    def make_cayley(a, b):
        H = Group.generate([a, b])
        print(len(H), len(G)//len(H))

        import string
        names = dict(zip(H, string.ascii_lowercase[:len(H)]))
        f = open("cayley_graph.dot", "w")
        print("digraph\n{\n", file=f)
        for g in H:
            print("    %s -> %s [color=red];"% (names[g], names[a*g]), file=f)
            print("    %s -> %s [color=green];"% (names[g], names[b*g]), file=f)
        print("}", file=f)
        f.close()

    if argv.cayley_graph:
        #make_cayley(a, b)
        make_cayley(c, b)

    if 0:
        H = Group.generate([b, c])
        gset = G.action_subgroup(H)
        print(gset.get_orbits())

    #H = Group.generate([c])
    #o_edges = G.left_cosets(H)
    #print(len(o_edges))

    if argv.burnside:
        burnside(G)

    for H in G.subgroups():
        if len(H)==60:
            assert H.is_simple()
            break



if __name__ == "__main__":

    main()
    test_gaussian()



