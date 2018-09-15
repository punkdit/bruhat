#!/usr/bin/env python3

"""
"""

import sys, os

from action import mulclose, Perm, Group

from util import cross
from argv import argv


class Type(object):
    pass



class Element(object):
    def __init__(self, tp):
        assert isinstance(tp, Type)
        self.tp = tp
        #self.promote = tp.promote # ??

    def __add__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.add(self, other)
        return a

    def __sub__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.sub(self, other)
        return a

    def __rsub__(self, other):
        tp = self.tp
        other = tp.promote(other)
        a = tp.sub(other, self)
        return a

    def __mul__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.mul(self, other)
        return a

    def __rmul__(self, value):
        tp = self.tp
        other = tp.promote(value)
        a = tp.mul(other, self)
        return a

    def __radd__(self, value):
        tp = self.tp
        other = tp.promote(value)
        a = tp.add(other, self)
        return a

    def __neg__(self):
        tp = self.tp
        a = tp.neg(self)
        return a

    def __truediv__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.div(self, other)
        return a

    def __rtruediv__(self, other):
        tp = self.tp
        other = tp.promote(other)
        a = tp.div(other, self)
        return a

    def __pow__(self, n):
        assert int(n)==n
        assert n>=0
        p = self.tp.one
        for i in range(n):
            p = self*p
        return p


class Keyed(object):

    def __init__(self, key):
        self.key = key

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        if self.__class__ is not other.__class__:
            return False
        return self.key == other.key

    def __ne__(self, other):
        if self.__class__ is not other.__class__:
            return True
        return self.key != other.key

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.key)
    __repr__ = __str__



class GenericElement(Element):
    """
        Element with immutable, hashable cargo 
    """

    def __init__(self, value, tp):
        Element.__init__(self, tp)
        self.value = value

    def __eq__(self, other):
        tp = self.tp
        other = tp.promote(other)
        return self.value == other.value

    def __ne__(self, other):
        tp = self.tp
        other = tp.promote(other)
        return self.value != other.value

    def __hash__(self):
        return hash((self.tp, self.value))

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return "%s(%r, %r)"%(self.__class__.__name__, self.value, self.tp)



# ----------------------------------------------------------------------------


class Ring(Type):

    def neg(self, a):
        zero = self.zero
        a = self.sub(zero, a)
        return a

    def __repr__(self):
        return "%s()"%(self.__class__.__name__,)

    def __hash__(self):
        return hash(self.__class__)

    def __eq__(self, other):
        assert isinstance(other, Type)
        return self.__class__ == other.__class__ # ?

    def __ne__(self, other):
        assert isinstance(other, Type)
        return self.__class__ != other.__class__ # ?


class IntegerRing(Ring):

    """
        The ring of integers.
    """

    # XXX these should be singletons (so we can test equality with "is")

    def __init__(self):
        Ring.__init__(self)
        self.one = Integer(1, self)
        self.zero = Integer(0, self)

    def add(self, a, b):
        assert a.tp is self # do we need this here?
        assert b.tp is self # do we need this here?
        return Integer(a.i + b.i, self)
    
    def sub(self, a, b):
        return Integer(a.i - b.i, self)
    
    def mul(self, a, b):
        return Integer(a.i * b.i, self)
    
    def neg(self, a):
        return Integer(-a.i, self)

    def div(self, a, b):
        if a.i%b.i:
            raise Exception
        i = a.i // b.i
        return Integer(i, self)
    
    def promote(self, value):
        if isinstance(value, Integer):
            assert value.tp is self
            return value
        if isinstance(value, int):
            return Integer(value, self)
        return None


class FiniteField(Ring):

    """
        The field of integers modulo a prime.
    """

    # XXX these should be singletons (so we can test equality with "is")

    def __init__(self, p):
        Ring.__init__(self)
        self.p = p
        self.one = FieldElement(1, self)
        self.zero = FieldElement(0, self)

    def nonzero_els(self):
        for i in range(1, self.p):
            yield FieldElement(i, self)

    def get_elements(self):
        for i in range(self.p):
            yield FieldElement(i, self)

    def __hash__(self):
        return hash((self.__class__, self.p))

    def __eq__(self, other):
        if self.__class__ is not other.__class__:
            return False
        return self.p == other.p

    def __ne__(self, other):
        if self.__class__ is not other.__class__:
            return True
        return self.p != other.p

    def add(self, a, b):
        assert a.tp is self
        assert b.tp is self
        p = self.p
        return FieldElement((a.i + b.i)%p, self)
    
    def sub(self, a, b):
        p = self.p
        return FieldElement((a.i - b.i)%p, self)
    
    def mul(self, a, b):
        p = self.p
        return FieldElement((a.i * b.i)%p, self)
    
    def neg(self, a):
        p = self.p
        return FieldElement((p-a.i)%p, self)
    
    def inverse(self, a):
        p = self.p
        assert 0<a.i<p
        for j in range(1, p):
            if (j*a.i)%p == 1:
                break
        else:
            assert 0
        return FieldElement(j, self)

    def div(self, a, b):
        b = self.inverse(b)
        return self.mul(a, b)
    
    def promote(self, value):
        if isinstance(value, FieldElement):
            assert value.tp is self
            return value
        if not isinstance(value, int):
            return None
        value = value % self.p
        return FieldElement(value, self)


class Integer(Element):
    def __init__(self, i, tp): # XXX make all Element's use generic "value" attr ?? XXX
        Element.__init__(self, tp)
        assert int(i)==i
        self.i = i

    def __eq__(self, other):
        tp = self.tp
        other = tp.promote(other)
        return self.i == other.i

    def __ne__(self, other):
        tp = self.tp
        other = tp.promote(other)
        return self.i != other.i

    def __hash__(self):
        return hash((self.tp, self.i))

    def __str__(self):
        return str(self.i)

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.i)


class FieldElement(Integer):
    def __init__(self, i, tp):
        Element.__init__(self, tp)
        assert int(i)==i
        assert 0<=i<tp.p
        self.i = i



# ----------------------------------------------------------------------------


class PolynomialRing(Ring):
    """
        Ring of polynomials, over some other base ring.
    """

    def __init__(self, base):
        Ring.__init__(self)
        self.base = base
        self.zero = Polynomial({}, self)
        one = self.base.one
        self.one = Polynomial({0:one}, self)
        self.x = Polynomial({1:one}, self)

    def __hash__(self):
        return hash((self.__class__, self.base))

    def __eq__(self, other):
        assert self.__class__ is other.__class__
        return self.base == other.base

    def __ne__(self, other):
        assert self.__class__ is other.__class__
        return self.base != other.base

    def __div__(self, mod):
        assert mod.tp == self
        return ModuloRing(self, mod)

    def promote(self, value):
        if isinstance(value, Polynomial):
            assert value.tp == self
            return value
        if isinstance(value, Element) and value.tp == self:
            return value
        _value = self.base.promote(value)
        if _value is None:
            #print(self, "cannot promote", value)
            return None
        return Polynomial({0:_value}, self)

    def add(self, a, b):
        cs = dict(a.cs)
        for deg, coeff in b.cs.items():
            cs[deg] = cs.get(deg, 0) + coeff
        return Polynomial(cs, self)

    def sub(self, a, b):
        cs = dict(a.cs)
        for deg, coeff in b.cs.items():
            cs[deg] = cs.get(deg, 0) - coeff
        return Polynomial(cs, self)

    def mul(self, a, b):
        cs = {}
        for d0, c0 in a.cs.items():
          for d1, c1 in b.cs.items():
            deg = d0+d1
            cs[deg] = cs.get(deg, 0) + c0*c1
        return Polynomial(cs, self)



class Polynomial(Element):
    def __init__(self, cs, tp):
        Element.__init__(self, tp)
        self.cs = {} # use a tuple here and switch to GenericElement ?
        self.deg = -1 # zero Polynomial has degree -1... right?
        for deg, coeff in cs.items():
            if coeff==0: # strip these out
                continue
            _coeff = tp.base.promote(coeff)
            assert _coeff is not None, repr(cs)
            self.cs[deg] = _coeff
            self.deg = max(self.deg, deg)
        self.base = tp.base

    def __str__(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort(reverse=True)
        terms = []
        for deg in keys:
            coeff = cs[deg]
            assert coeff!=0
            if coeff==1 and deg==1:
                term = "x"
            elif deg==0:
                term = str(coeff)
            elif deg==1:
                term = "%s*x" % (coeff,)
            elif coeff==1:
                term = "x**%d" % deg
            else:
                term = "%s*x**%d" % (coeff, deg)
            terms.append(term)
        s = "+".join(terms) or "0"
        s = s.replace("+-", "-")
        s = s.replace("-1*", "-")
        return s

    def __repr__(self):
        return "Polynomial(%s, %s)"%(self.cs, self.tp)

    def __eq__(self, other): # XXX put this in the tp ?
        other = self.tp.promote(other)
        assert self.tp == other.tp
        return self.cs == other.cs

    def __ne__(self, other): # XXX put this in the tp ?
        other = self.tp.promote(other)
        assert self.tp == other.tp
        return self.cs != other.cs

    def __hash__(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()
        value = tuple((deg, cs[deg]) for deg in keys)
        return hash(value)

    def __getitem__(self, idx):
        return self.cs.get(idx, self.base.zero)

    def __call__(self, v):
        v = self.base.promote(v)
        r = self.base.zero
        for deg, coeff in self.cs.items():
            r = r + coeff * v**deg
        return r

    def reduce(a, b):
        """
            a.reduce(b): use a to reduce b.
            return (div, rem) such that div*a+rem==b.
        """
        tp = a.tp
        base = a.base
        if b.deg < a.deg:
            # already reduced
            return base.zero, b
        x = tp.x
        r = base.zero
        #print("reduce: %s, %s" % (a, b))
        assert a != 0
        while b.deg >= a.deg:
            b0 = b[b.deg]
            assert b0 is not None, b
            a0 = a[a.deg]
            assert a0 is not None, repr(a)
            coeff = b0/a0
            assert coeff != 0
            m = coeff * x**(b.deg - a.deg)
            #print("m = %s"%m)
            ma = m*a
            assert ma.deg == b.deg
            _b = b - ma
            #print("b = %s"%_b)
            assert _b.deg < b.deg
            b = _b
            r = r + m
        return r, b


class ModuloRing(Keyed, Ring):
    """
        Ring modulo a principle ideal.
    """

    def __init__(self, ring, mod):
        Ring.__init__(self)
        self.ring = ring
        assert mod.tp == ring
        self.mod = mod
        key = (self.__class__, self.ring, self.mod)
        Keyed.__init__(self, key)
        self.reduce = mod.reduce
        self.zero = ModuloElement(ring.zero, self)
        self.one = ModuloElement(ring.one, self)

#    def __hash__(self):
#        return hash(self.key)
#
#    def __eq__(self, other):
#        if self.__class__ is not other.__class__:
#            return False
#        return self.key == other.key
#
#    def __ne__(self, other):
#        if self.__class__ is not other.__class__:
#            return True
#        return self.key != other.key

    def promote(self, value):
        if isinstance(value, Element) and value.tp==self:
            return value
        value = self.ring.promote(value)
        value = self.reduce(value)[1]
        a = ModuloElement(value, self)
        return a

    def add(self, a, b):
        ring = self.ring
        value = ring.add(a.value, b.value)
        value = self.reduce(value)[1]
        a = ModuloElement(value, self)
        return a

    def sub(self, a, b):
        ring = self.ring
        value = ring.sub(a.value, b.value)
        value = self.reduce(value)[1]
        a = ModuloElement(value, self)
        return a

    def mul(self, a, b):
        ring = self.ring
        value = ring.mul(a.value, b.value)
        value = self.reduce(value)[1]
        a = ModuloElement(value, self)
        return a


class GaloisField(ModuloRing):
    def __init__(self, ring, mod):
        ModuloRing.__init__(self, ring, mod)
        assert isinstance(ring, PolynomialRing)
        assert mod[mod.deg] == 1 # monic
        dim = mod.deg
        assert dim>1
        one = ring.one # unwrapped
        x = ring.x # unwrapped
        basis = [one, x]
        for i in range(dim-2):
            basis.append(basis[-1]*x)
        self.basis = basis # unwrapped
        elements = []
        els = list(ring.base.get_elements())
        for idx in cross((els,)*dim):
            a = ring.zero # unwrapped
            for i, coeff in enumerate(idx):
                a = a + coeff * basis[i]
            a = ModuloElement(a, self) # wrap
            elements.append(a)
        self.elements = elements
        assert len(elements) == len(els) ** dim, len(elements)

        self.zero = ModuloElement(ring.zero, self)
        self.one = ModuloElement(one, self)
        self.x = ModuloElement(x, self)

    def div(self, a, b):
        ring = self.ring
        mod = self.mod
        if a==self.zero:
            return self.zero
        if b==self.zero:
            return None # fail
        for div in self.elements:
            if a == b*div:
                break
        else:
            assert 0
        return div


class ModuloElement(GenericElement):
    @property
    def deg(self):
        return self.value.deg


# ----------------------------------------------------------------------------


class Linear(Keyed, Type):
    
    def __init__(self, n, base):
        "Algebraic group: (n, n) matrices over a base ring"
        self.n = n
        self.base = base
        key = (self.n, self.base)
        Keyed.__init__(self, key)
        Type.__init__(self)
        one = self.base.one
        zero = self.base.zero
        self.zero = LinearElement(
            tuple(tuple(zero for j in range(n)) for i in range(n)), self)
        self.one = LinearElement(
            tuple(tuple(
            (one if i==j else zero)
            for j in range(n)) for i in range(n)), self)

    def promote(self, value):
        if isinstance(value, Element) and value.tp==self:
            return value
        if isinstance(value, Element) and isinstance(value.tp, Linear):
            assert 0, "%s cannot promote from %s" % (self, value.tp)
            return None
        value = self.base.promote(value)
        n = self.n
        zero = self.base.zero
        a = tuple(tuple(
            (value if i==j else zero)
            for j in range(n)) for i in range(n))
        a = LinearElement(a, self)
        return a
    
    def add(self, a, b):
        a = a.value # unwrap
        b = b.value # unwrap
        n = self.n
        value = tuple(tuple(
                a[i][j] + b[i][j]
            for j in range(n)) for i in range(n))
        return LinearElement(value, self)
    
    def sub(self, a, b):
        a = a.value # unwrap
        b = b.value # unwrap
        n = self.n
        value = tuple(tuple(
                a[i][j] - b[i][j]
            for j in range(n)) for i in range(n))
        return LinearElement(value, self)

    def neg(self, a):
        a = a.value # unwrap
        n = self.n
        value = tuple(tuple(-a[i][j]
            for j in range(n)) for i in range(n))
        return LinearElement(value, self)

    def mul(self, a, b):
        a = a.value # unwrap
        b = b.value # unwrap
        n = self.n
        value = tuple(tuple(
                sum(a[i][k] * b[k][j] for k in range(n))
            for j in range(n)) for i in range(n))
        return LinearElement(value, self)

    def get(self, value):
        n = self.n
        assert len(value)==n
        for row in value:
            assert len(row)==n
        base = self.base
        value = tuple(tuple(
            base.promote(value[i][j]) 
            for j in range(n)) for i in range(n))
        return LinearElement(value, self)


class LinearElement(GenericElement):
    def __init__(self, value, tp):
        n = tp.n
        assert len(value)==n
        assert len(value[0])==n
        base = tp.base
        value = tuple(tuple(
            base.promote(value[i][j]) 
            for j in range(n)) for i in range(n))
        GenericElement.__init__(self, value, tp)

    def __getitem__(self, key):
        i, j = key
        return self.value[i][j]


# ----------------------------------------------------------------------------


def test():

    f = FiniteField(5)
    one = f.one
    zero = f.zero

    assert one==1
    assert str(one) == "1"

    a = zero
    for i in range(1000):
        a = one + a
        if a==zero:
            break
    assert i==4, i

    # -------------------------

    Z = IntegerRing()
    one = Z.one
    two = one + one
    assert one/one == one
    assert two/two == one
    assert two/(-one) == -two

    ring = PolynomialRing(Z)

    one = ring.one
    assert ring.x == Polynomial({1:1}, ring)
    x = ring.x

    assert (x+one)**5 == x**5+5*x**4+10*x**3+10*x**2+5*x+1
    assert str((x+one)**3) == "x**3+3*x**2+3*x+1"

    a = x**3 + x + 1
    b = x**6
    div, rem = a.reduce(b)
    assert a*div + rem == b
    #print("(%s) * (%s) + %s = %s" % (a, div, rem, b))

    # -------------------------

    field = FiniteField(5)
    ring = PolynomialRing(field)

    one = ring.one
    x = ring.x

    assert (x+one)**5 == x**5+1

    assert x**0 == one

    # -------------------------

    # check we can build galois field GF(8)

    field = FiniteField(2)
    ring = PolynomialRing(field)

    one = ring.one
    x = ring.x

    f = x**3 - x - 1
    assert f == x**3+x+1
    assert f(0) != 0
    assert f(1) != 0
    
    b = x**5
    div, rem = f.reduce(b)
    assert f*div + rem == b

    group = []
    for i in range(2):
     for j in range(2):
      for k in range(2):
        a = i + j*x + k*x**2
        if a != 0:
            group.append(a)

    div = {}
    for a in group:
      for b in group:
        c = f.reduce(a*b)[1]
        div[c, a] = b
        div[c, b] = a

    # all non-zero pairs elements of GF(8) should be divisable:
    assert len(div) == 49

    # --------------------------------------------------

    # GF(4)
    x = ring.x
    GF = GaloisField(ring, (x**2 + x + 1))

    x = GF.x
    one = GF.one
    a, b = x, x+one
    assert a*a == b
    assert a*b == b*a
    assert a*b == one
    assert b*b == a

    assert one/x == one+x

    # --------------------------------------------------

    # GF(8)
    x = ring.x
    GF = GaloisField(ring, (x**3 + x + 1))

    one = GF.one
    x = GF.x

    assert one/x == x**2+1

    # --------------------------------------------------

    # GF(64)
    x = ring.x
    GF = GaloisField(ring, (x**6 + x + 1))

    one = GF.one
    x = GF.x

    assert len(GF.elements) == 64

    # -------------------------

    p = 7
    field = FiniteField(p)
    ring = PolynomialRing(field)
    x = ring.x

    if p % 4 == 3:
        r = p-1
    else:
        assert 0

    GF = GaloisField(ring, x**2 - r)

    one = GF.one
    x = GF.x

    for a in GF.elements:
        b = a**(p+1)
        assert b.deg <= 0, repr(b)

    # -------------------------

    Z2 = FiniteField(2)
    for GL in [Linear(1, Z), Linear(2, Z), Linear(2, Z2)]:

        one = GL.one
        zero = GL.zero
        
        assert one-one == zero
        assert zero+one == one
        assert one+one != one
        assert one*one == one
        assert zero*one == zero
        assert one+one == 2*one

        if GL.base == Z2:
            assert one+one == zero
    
    GL2 = Linear(2, Z)
    A = GL2.get([[1, 1], [0, 1]])

    # -------------------------

    field = FiniteField(3)
    GL = Linear(2, field)
    A = GL.get([[1, 1], [0, 1]])
    B = GL.get([[1, 0], [1, 1]])

    SL2_3 = mulclose([GL.one, A, B])
    assert len(SL2_3) == 24
    
    # -------------------------

    field = FiniteField(7)
    GL = Linear(2, field)
    A = GL.get([[1, 1], [0, 1]])
    B = GL.get([[1, 0], [1, 1]])

    SL2_7 = mulclose([GL.one, A, B])
    assert len(SL2_7) == 336 == 168*2
    



if __name__ == "__main__":

    test()

    print("OK")





