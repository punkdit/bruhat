#!/usr/bin/env python3

"""
"""

import sys, os

from bruhat.action import mulclose, Perm, Group, burnside

from bruhat.util import cross
from bruhat.argv import argv


def divisors(n):
    divs = [1] 
    for i in range(2, n): 
        if n%i == 0:
            divs.append(i)
    if n>1:
        divs.append(n)
    return divs


class Type(object):

    def promote(self, a):
        if isinstance(a, Element) and a.tp == self:
            return a
        return None

    def eq(self, a, b):
        return a.value == b.value

    def ne(self, a, b):
        return a.value != b.value

    def bool(self, a):
        return True

    def gcd(self, a, b):
        raise TypeError("Not implemented for Type %s"%self)

    def truediv(self, a, b):
        raise TypeError("Not implemented for Type %s"%self)

    def floordiv(self, a, b):
        raise TypeError("Not implemented for Type %s"%self)

    # etc...


class Element(Type):
    def __init__(self, tp):
        assert isinstance(tp, Type)
        self.tp = tp
        #self.promote = tp.promote # ??

    def __eq__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        return tp.eq(self, other)

    def __ne__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        return tp.ne(self, other)

    def __bool__(self):
        return self.tp.bool(self)

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
        if other is None:
            return NotImplemented
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
        if other is None:
            return NotImplemented
        a = tp.mul(other, self)
        return a

    def __matmul__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.matmul(self, other)
        return a

    def __rmatmul__(self, value):
        tp = self.tp
        other = tp.promote(value)
        if other is None:
            return NotImplemented
        a = tp.matmul(other, self)
        return a

    def __radd__(self, value):
        tp = self.tp
        other = tp.promote(value)
        if other is None:
            return NotImplemented
        a = tp.add(other, self)
        return a

    def __neg__(self):
        tp = self.tp
        a = tp.neg(self)
        return a

    def __floordiv__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.floordiv(self, other)
        return a

    def __rfloordiv__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.floordiv(other, self)
        return a

    def __truediv__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.truediv(self, other)
        return a

    def __rtruediv__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.truediv(other, self)
        return a

    def __mod__(self, other):
        tp = self.tp
        other = tp.promote(other)
        if other is None:
            return NotImplemented
        a = tp.mod(self, other)
        return a

    def __pow__(self, n):
        assert int(n)==n
        assert n>=0
        p = self.tp.one
        for i in range(n):
            p = self*p
        return p


class Keyed(object): # use for Type's

    def __init__(self, key):
        if not type(key) is tuple:
            key = (key,)
        self.key = (self.__class__,) + key

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.key == other.key

    def __ne__(self, other):
        return self.__class__ != other.__class__ or self.key != other.key

    def __str__(self):
        return str(self.key)
    __repr__ = __str__



class GenericElement(Element):
    """
        Element with immutable, hashable cargo 
    """

    def __init__(self, value, tp):
        Element.__init__(self, tp)
        self.value = value

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

    def bool(self, a):
        return a != self.zero


class IntegerRing(Ring):

    """
        The ring of integers.
    """

    # XXX these should be singletons (so we can test equality with "is")

    def __init__(self):
        Ring.__init__(self)
        self.one = Integer(1, self)
        self.zero = Integer(0, self)
    
    def promote(self, value):
        if isinstance(value, Integer):
            assert value.tp == self
            return value
        if isinstance(value, int):
            return Integer(value, self)
        try:
            value = int(value)
            return Integer(value, self)
        except ValueError:
            pass
        except TypeError:
            pass
        return None

    def add(self, a, b):
        assert a.tp is self # do we need this here?
        assert b.tp is self # do we need this here?
        return Integer(a.value + b.value, self)
    
    def sub(self, a, b):
        return Integer(a.value - b.value, self)
    
    def mul(self, a, b):
        return Integer(a.value * b.value, self)
    
    def neg(self, a):
        return Integer(-a.value, self)

    def truediv(self, a, b):
        if a.value%b.value:
            raise Exception("cannot divide %r by %r" % (a, b))
        i = a.value // b.value
        return Integer(i, self)

    def floordiv(self, a, b):
        i = a.value // b.value
        return Integer(i, self)

    def mod(self, a, b):
        a = a.value
        b = b.value
        return Integer(a%b, self)
    
    def gcd(self, a, b):
        assert a.tp == self, repr(a)
        assert b.tp == self, repr(b)
        a = a.value
        b = b.value
        while b != 0:
            a, b = b, a%b
        factor = Integer(a, self)
        return factor


class Integer(GenericElement):

    def reduce(a, b):
        """
            a.reduce(b): use a to reduce b.
            return (div, rem) such that div*a+rem==b.
        """

        a = a.value # unwrap
        b = b.value # unwrap

        div = a//b
        rem = a%b
        return div, rem



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
        return FieldElement((a.value + b.value)%p, self)
    
    def sub(self, a, b):
        p = self.p
        return FieldElement((a.value - b.value)%p, self)
    
    def mul(self, a, b):
        p = self.p
        return FieldElement((a.value * b.value)%p, self)
    
    def neg(self, a):
        p = self.p
        return FieldElement((p-a.value)%p, self)
    
    def inverse(self, a):
        p = self.p
        assert 0<a.value<p
        for j in range(1, p):
            if (j*a.value)%p == 1:
                break
        else:
            assert 0
        return FieldElement(j, self)

    def truediv(self, a, b):
        b = self.inverse(b)
        return self.mul(a, b)
    floordiv = truediv

    def mod(self, a, b):
        return self.zero
    
    def promote(self, value):
        if isinstance(value, FieldElement):
            assert value.tp is self
            return value
        #if not isinstance(value, int):
        #    return None
        try:
            value = int(value)
        except:
            return None
        value = value % self.p
        return FieldElement(value, self)



class FieldElement(Integer):
    def __init__(self, value, tp):
        GenericElement.__init__(self, value, tp)
        assert int(value)==value
        assert 0<=value<tp.p


# ----------------------------------------------------------------------------


class PolynomialRing(Keyed, Ring):
    """
        Ring of polynomials, over some other base ring.
    """

    def __init__(self, base):
        Keyed.__init__(self, base)
        Ring.__init__(self)
        self.base = base
        self.zero = Polynomial({}, self)
        one = self.base.one
        self.one = Polynomial({0:one}, self)
        self.x = Polynomial({1:one}, self)

    def __truediv__(self, mod):
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

    def eq(self, a, b):
        return a.cs == b.cs

    def ne(self, a, b):
        return a.cs != b.cs

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

    def truediv(self, a, b):
        p, rem = b.reduce(a)
        if rem != self.zero:
            return None
        return p

    def floordiv(self, a, b):
        p, rem = b.reduce(a)
        return p

    def mod(self, a, b):
        div, rem = b.reduce(a)
        return rem

    def content(self, a):
        # Definition 3.2.7, p116 (Cohen, 1993)
        cs = list(a.cs.values())
        ring = self.base
        if not cs:
            return ring.zero
        d = cs[0]
        for e in cs[1:]:
            d = ring.gcd(d, e)
        return d

    def primitive(self, A):
        if A==self.zero:
            return self.zero
        a = self.content(A)
        return A/a

#    def gcd(self, a, b):
#        assert a.tp == self
#        assert b.tp == self
#        args = (a, b)
#        print("gcd(%s, %s):"%(a, b))
#        seen = set([(a, b)])           # remove me
#        while b != self.zero:
#            a, b = b, a%b
#            print("\t%s, %s"%(a, b))
#            assert (a, b) not in seen, "gcd(%s, %s): %s"%args  # remove me
#            seen.add((a, b))           # remove me
#            #print("gcd", a, b)
#        return a

    def pseudo_div(self, A, B):
        # algorithm 3.1.2, p112 (Cohen, 1993)
        assert A.tp == self
        assert B.tp == self
        assert B!=0
        m, n = A.deg, B.deg
        d = B[n]
        zero = self.zero
        x = self.x
        Q = zero
        R = A
        while R.deg >= B.deg:
            S = R[R.deg] * x**(R.deg - B.deg)
            Q = d*Q + S
            R = d*R - S*B

        return Q, R

    def gcd(self, A, B):
        # algorithm 3.2.10, p117 (Cohen, 1993)
        assert A.tp == self
        assert B.tp == self
        args = (A, B)
#        print("gcd(%s, %s):"%(A, B))
        if B==0:
            return A
        ring = self.base
        a = self.content(A)
        b = self.content(B)
        d = ring.gcd(a, b)
        A = self.primitive(A)
        B = self.primitive(B)

        while 1:
#            print("\tA=%s, B=%s" % (A, B))
            #Q, R = B.reduce(A)
            Q, R = self.pseudo_div(A, B)
#            print("\tQ=%s, R=%s" % (Q, R,))
            if Q != 0:
                lhs = B[B.deg] ** (A.deg - B.deg + 1) * A
                rhs = B*Q + R
                assert lhs == rhs
            else:
                assert R==A
            if R==0:
                break
            if R.deg == 0:
                B = 1
                break
            A = B
            B = R/self.content(R)

        #print("\td=%s, B=%s"%(d, B), ring)
        #print("\t", type(d*B))
        B = self.promote(d*B)
#        print("\tB=%s"%B)
        assert B is not None
        return B



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
#        print("reduce: %s, %s" % (a, b))
        assert a != 0
        while b.deg >= a.deg:
            b0 = b[b.deg]
            assert b0 is not None, b
            a0 = a[a.deg]
            assert a0 is not None, repr(a)
            div = b0//a0
            if div==0:
                break
            rem = b0%a0
#            print("div = %s, rem = %s" % (div, rem))
            assert div != 0
            m = div * x**(b.deg - a.deg)
#            print("m = %s"%m)
            ma = m*a
            assert ma.deg == b.deg
            _b = b - ma
#            print("b = %s"%_b)
            assert _b.deg <= b.deg
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
        key = (self.ring, self.mod)
        Keyed.__init__(self, key)
        self.reduce = mod.reduce
        self.zero = ModuloElement(ring.zero, self)
        self.one = ModuloElement(ring.one, self)

    def promote(self, value):
        if isinstance(value, Element) and value.tp==self:
            return value
        value = self.ring.promote(value)
        if value is None:
            return None
        value = self.reduce(value)[1]
        a = ModuloElement(value, self)
        return a

    def __getattr__(self, attr):
        value = getattr(self.ring, attr)
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

    def gcd(self, a, b):
        assert a.tp == self
        assert b.tp == self
        ring = self.ring
        a = ring.gcd(a.value, b.value)
        a = self.promote(a)
        return a

    def truediv(self, a, b):
        assert a.tp == self
        assert b.tp == self
        ring = self.ring
        a = a.value / b.value
        #a = ring.truediv(a.value, b.value)
        a = self.promote(a)
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

    def truediv(self, a, b):
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


class FieldOfFractions(Ring):
    
    def __init__(self, base):
        Ring.__init__(self)
        self.base = base
        one = base.one
        zero = base.zero
        self.one = Fraction((one, one), self)
        self.zero = Fraction((zero, one), self)

    def promote(self, value):
        if isinstance(value, Fraction):
            assert value.tp == self
            return value
        value = self.base.promote(value)
        if value is None:
            return None
        value = Fraction((value, self.base.one), self)
        return value

    def eq(self, a, b):
        atop, abot = a.value
        btop, bbot = b.value
        return atop * bbot == btop * abot
    
    def ne(self, a, b):
        atop, abot = a.value
        btop, bbot = b.value
        return atop * bbot != btop * abot
    
    def add(self, a, b):
        atop, abot = a.value
        btop, bbot = b.value
        top = atop * bbot + btop * abot
        bot = abot * bbot
        return Fraction((top, bot), self)
    
    def sub(self, a, b):
        atop, abot = a.value
        btop, bbot = b.value
        top = atop * bbot - btop * abot
        bot = abot * bbot
        return Fraction((top, bot), self)
    
    def mul(self, a, b):
        atop, abot = a.value
        btop, bbot = b.value
        top = atop * btop
        bot = abot * bbot
        return Fraction((top, bot), self)
    
    def neg(self, a):
        atop, abot = a.value
        return Fraction((-atop, abot), self)

    def truediv(self, a, b):
        atop, abot = a.value
        btop, bbot = b.value
        top = atop * bbot
        bot = abot * btop
        if bot == self.zero:
            raise Exception("cannot divide %s by %s"%(a, b))
        return Fraction((top, bot), self)
    

class Fraction(GenericElement):
    def __init__(self, value, tp):
        top, bot = value
        base = tp.base
        assert top.tp == base
        assert bot.tp == base
        factor = base.gcd(top, bot)
        top = top/factor
        bot = bot/factor
        value = (top, bot)
        GenericElement.__init__(self, value, tp)

    def __str__(self):
        top, bot = self.value
        one = self.tp.base.one
        if bot==one:
            return str(top)
        return "(%s/%s)"%(top, bot)




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
            # XXX TODO
            assert 0, "TODO: %s cannot promote from %s" % (self, value.tp)
            return None
        value = self.base.promote(value)
        if 0:
            n = self.n
            zero = self.base.zero
            a = tuple(tuple(
                (value if i==j else zero)
                for j in range(n)) for i in range(n))
            a = LinearElement(a, self)
            return a
        return value
    
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
        n = self.n
        if a.tp == self.base:
            b = b.value # unwrap
            value = tuple(tuple(a * b[i][j] for j in range(n)) for i in range(n))
        elif b.tp == self.base:
            a = a.value # unwrap
            value = tuple(tuple(a[i][j] * b for j in range(n)) for i in range(n))
        else:
            a = a.value # unwrap
            b = b.value # unwrap
            value = tuple(tuple(
                    sum(a[i][k] * b[k][j] for k in range(n))
                for j in range(n)) for i in range(n))
        return LinearElement(value, self)

    def truediv(self, a, b):
        assert b.tp == self.base
        a = a.value # unwrap
        n = self.n
        value = tuple(tuple(a[i][j] / b for j in range(n)) for i in range(n))
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

    def __getitem__(self, key):
        i, j = key
        return self.value[i][j]

    def __str__(self):
        n = self.tp.n
        value = self.value
        rows = [[str(value[i][j]) for j in range(n)] for i in range(n)]
        w = 1
        for row in rows:
            for col in row:
                w = max(w, len(col))
        rows = ['[%s]'%' '.join(s.rjust(w) for s in row) for row in rows]
        lines = []
        for i, row in enumerate(rows):
            if i==0:
                row = "["+row
            else:
                row = " "+row
            if i==len(rows)-1:
                row = row+"]"
            else:
                row = row+","
            lines.append(row)
        return '\n'.join(lines)

    def trace(self):
        x = self.tp.base.zero
        value = self.value
        for i in range(self.tp.n):
            x = x + value[i][i]
        return x 


# ----------------------------------------------------------------------------


def cayley(elements):
    "build the regular permutation representation"
    elements = list(elements)
    lookup = {}
    for idx, e in enumerate(elements):
        lookup[e] = idx
    items = list(range(len(elements)))
    perms = []
    for i, e in enumerate(elements):
        perm = {}
        for j, f in enumerate(elements):
            g = e*f
            k = lookup[g]
            perm[j] = k
        perm = Perm(perm, items)
        perms.append(perm)
    G = Group(perms, items)
    return G


_cyclotomic_cache = {}
def cyclotomic(ring, n): # make this a method ?

    key = (ring, n)
    if key in _cyclotomic_cache:
        return _cyclotomic_cache[key]

    divs = divisors(n)

    x = ring.x
    one = ring.one

    if n==1:
        p = x-1

    elif len(divs)==2:
        # prime
        p = sum(x**i for i in range(n))

    else:

        p = x**n - 1

        assert divs[-1] == n
        for i in divs[:-1]:
            p = p / cyclotomic(ring, i)

    _cyclotomic_cache[key] = p
    return p

# ----------------------------------------------------------------------------

Z = IntegerRing()
Q = FieldOfFractions(Z)


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

    assert ring.content(10*x + 2) == 2
    assert ring.content(5*x**2 + 10*x + 2) == 1

    for (a, b) in [
        (4, 10),
        (4, 10*x),
        (2*x+1, 3*x),
        (x**3 + x + 1, x**6),
    ]:
        a = ring.promote(a)
        b = ring.promote(b)
        div, rem = a.reduce(b)
        assert a*div + rem == b
        #print("(%s) * (%s) + %s = %s" % (a, div, rem, b))

        c = ring.gcd(a, b)
        #print("gcd(%s, %s) == %s"%(a, b, c))

    #return

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

    # AKA. binary tetrahedral group
    SL2_3 = mulclose([A, B])
    assert len(SL2_3) == 24

    if argv.SL2_3:
        burnside(cayley(SL2_3))

    # -------------------------
    # https://people.maths.bris.ac.uk/~matyd/GroupNames/1/GL(2,3).html

    C = GL.get([[2, 0], [0, 1]])
    D = GL.get([[1, 0], [0, 2]])

    GL2_3 = mulclose([A, B, C, D])
    assert len(GL2_3) == 48

    if argv.GL2_3:
        burnside(cayley(GL2_3))

    # -------------------------
    # https://people.maths.bris.ac.uk/~matyd/GroupNames/1/Q8.html
    # quaternion group, Q_8 

    A = GL.get([[2, 2], [2, 1]])
    B = GL.get([[0, 2], [1, 0]])

    Q8 = mulclose([A, B])
    assert(len(Q8)) == 8

    if argv.Q8 or argv.Q_8:
        burnside(cayley(Q8))

    # -------------------------
    # https://people.maths.bris.ac.uk/~matyd/GroupNames/97/SL(2,5).html
    # aka binary icosahedral group

    field = FiniteField(5)
    GL = Linear(2, field)

    A = GL.get([[4, 2], [4, 1]])
    B = GL.get([[3, 3], [4, 1]])

    SL2_5 = mulclose([A, B])
    assert len(SL2_5) == 120

    if argv.SL2_5:
        burnside(cayley(SL2_5))

    # -------------------------

    field = FiniteField(7)
    GL = Linear(2, field)
    A = GL.get([[1, 1], [0, 1]])
    B = GL.get([[1, 0], [1, 1]])

    SL2_7 = mulclose([GL.one, A, B])
    assert len(SL2_7) == 336 == 168*2


    # https://people.maths.bris.ac.uk/~matyd/GroupNames/1/CSU(2,3).html
    # aka binary octahedral group
    A = GL.get([[2, 1], [2, 5]])
    B = GL.get([[1, 2], [6, 6]])
    C = GL.get([[4, 0], [2, 2]])
    D = GL.get([[2, 5], [6, 5]])

    BO = mulclose([A, B, C, D])
    assert len(BO) == 48
    if argv.BO:
        burnside(cayley(BO))

    # -------------------------

    ring = PolynomialRing(Z)
    x = ring.x
    assert cyclotomic(ring, 1) == x-1
    assert cyclotomic(ring, 2) == x+1
    assert cyclotomic(ring, 3) == x**2+x+1
    assert cyclotomic(ring, 4) == x**2+1
    assert cyclotomic(ring, 5) == x**4+x**3+x**2+x+1
    assert cyclotomic(ring, 6) == x**2-x+1
    assert cyclotomic(ring, 7) == x**6+x**5+x**4+x**3+x**2+x+1
    assert cyclotomic(ring, 8) == x**4+1
    assert cyclotomic(ring, 9) == x**6+x**3+1
    assert cyclotomic(ring, 10) == x**4-x**3+x**2-x+1
    assert cyclotomic(ring, 24) == x**8-x**4+1

    # -------------------------

    for n in range(2, 9):

        p = cyclotomic(PolynomialRing(Z), 2*n)
    
        ring = PolynomialRing(Z) / p
        x = ring.x
    
        x1 = x**(2*n-1) # inverse
        
        assert x**(2*n) == 1
        assert x1 != 1
    
        # ----------------------------
        # Binary dihedral group
        # https://people.maths.bris.ac.uk/~matyd/GroupNames/dicyclic.html
    
        GL = Linear(2, ring)
    
        A = GL.get([[x, 0], [0, x1]])
        B = GL.get([[0, -1], [1, 0]])
    
        Di = mulclose([A, B])
        assert len(Di) == 4*n

        if argv.Dic and n==argv.get("n"):
            burnside(cayley(Di))

    # -------------------------

    one = Q.one
    two = one+one
    assert two*one == two
    half = one/two
    assert 2*half == 1

    assert 4*one/4 == one
    assert str(18*one/4) == "(9/2)"

    # -------------------------

    ring = PolynomialRing(Z)
    field = FieldOfFractions(ring)
    one = field.one
    two = one+one
    assert two*one == two
    half = one/two
    assert 2*half == 1
    assert 4*one/4 == one

    # -------------------------
    # Gaussian integers
    # http://math.ucr.edu/home/baez/week216.html

    ring = PolynomialRing(Z)
    x = ring.x
    gints = ring / (x**2 + 1)

    one = gints.one
    i = gints.x

    assert i*i == -one
    assert i**3 == -i

    assert (1+i)**2 == 2*i
    assert (2+i)*(2-i) == 5
    assert (3+2*i)*(3-2*i) == 13

    # -------------------------
    # Eisenstein integers
    
    ring = PolynomialRing(Z)
    x = ring.x
    ints = ring / cyclotomic(ring, 3)

    one = ints.one
    i = ints.x

    assert i**3 == 1

    # -------------------------
    # Kleinian integers
    
    ring = PolynomialRing(Z)
    x = ring.x
    ints = ring / (x**2 + x + 2)

    one = ints.one
    i = ints.x

    assert i * (-1-i) == 2
    assert i * (1+i) == -2

    # -------------------------
    #

    ring = PolynomialRing(Z)
    p = cyclotomic(PolynomialRing(Z), 8)
    ring = ring / p
    field = FieldOfFractions(ring)

    x = field.promote(ring.x)
    one = field.one

    r2 = x + x**7
    assert r2**2 == 2  # sqrt(2)

    i = x**2
    assert i**2 == -one # sqrt(-1)

    GL = Linear(2, field)
    I = GL.get([[1, 0], [0, 1]])
    J = GL.get([[1, 0], [0, i]])
    
    M = (r2/2) * GL.get([[1, 1], [1, -1]]) # Hadamard

    assert J**4 == I
    assert M**2 == I
    
    G = mulclose([J, M])
    assert len(G) == 192

    assert M in G

    def get_order(A):
        count = 1
        B = A
        while B!=I:
            B = A*B
            count += 1
        return count
            
    #for g in G:
    #    if get_order(g)==4:
    #        print(g)
    assert get_order(M)==2

    # Metaplectic Hadamard
    M1 = GL.get([[1, -i], [-i, 1]])
    r = x / r2
    M1 = r*M1

    assert get_order(M1)==4

#    #G = cayley(G)
#    #orders = [g.order() for g in G]
#    orders = [get_order(A) for A in G]
#    orders.sort()
#    print(orders)



if __name__ == "__main__":

    test()

    print("OK")





