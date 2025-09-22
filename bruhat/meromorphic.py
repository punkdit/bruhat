#!/usr/bin/env python

"""

"""


from bruhat.argv import argv
from bruhat.poly import Poly
from bruhat.action import mulclose

from bruhat.element import Z, PolynomialRing


if 0:
    ring = Z

else:
    base = PolynomialRing(Z, "I")
    x = base.x
    ring = base / (x**2 + 1)

    one = ring.one
    I = ring.x

    # XXX FAIL
    # XXX need to make a GaussianIntegerRing, etc.





class Fraction:
    def __init__(self, top, bot=1):
        top = Poly.promote(top, ring)
        bot = Poly.promote(bot, ring)
        r, s = bot.reduce(top)
        if s==0:
            top = r
            bot = Poly.promote(1, ring)
        self.top = top
        self.bot = bot

    def __str__(self):
        if self.bot == 1:
            s = str(self.top)
        else:
            s = "(%s)/(%s)"%(self.top, self.bot)
        return s
    __repr__ = __str__

    @classmethod
    def promote(cls, item):
        if isinstance(item, Fraction):
            return item
        return Fraction(item)

    def __hash__(self):
        return hash((self.top, self.bot))

    def __eq__(self, other):
        other = Fraction.promote(other)
        return self.top*other.bot == other.top*self.bot

    def __mul__(self, other):
        other = Fraction.promote(other)
        top = self.top * other.top
        bot = self.bot * other.bot
        return Fraction(top, bot)
    __rmul__ = __mul__

    def __add__(self, other):
        other = Fraction.promote(other)
        atop, abot = self.top, self.bot
        btop, bbot = other.top, other.bot
        top = atop * bbot + btop * abot
        bot = abot * bbot
        return Fraction(top, bot)
    __radd__ = __add__

    def __sub__(self, other):
        other = Fraction.promote(other)
        atop, abot = self.top, self.bot
        btop, bbot = other.top, other.bot
        top = atop * bbot - btop * abot
        bot = abot * bbot
        return Fraction(top, bot)

    def __rsub__(self, other):
        other = Fraction.promote(other)
        return other - self

    def __neg__(self):
        return Fraction(-self.top, self.bot)

    def __floordiv__(self, other):
        other = Fraction.promote(other)
        top = self.top * other.bot
        bot = self.bot * other.top
        return Fraction(top, bot)
    __truediv__ = __floordiv__

    def __pow__(self, n):
        top = self.top**n
        bot = self.bot**n
        return Fraction(top, bot)

    def get_vars(self):
        vs = self.top.get_vars() + self.bot.get_vars()
        vs = list(set(vs))
        vs.sort()
        return vs

    def substitute(self, items):
        ns = dict(items)
        ns["I"] = I
        etop = self.top.python_str()
        top = eval(etop, ns)
        ebot = self.bot.python_str()
        bot = eval(ebot, ns)
        top = Fraction.promote(top)
        bot = Fraction.promote(bot)
        return top/bot


class Meromorphic:
    def __init__(self, p):
        p = Fraction.promote(p)
        vs = p.get_vars()
        self.p = p
        self.vs = vs

    def __str__(self):
        return "Meromorphic(%s:%s)"%(self.p, self.vs)

    def __eq__(self, other):
        return self.p == other.p

    def __hash__(self):
        return hash(self.p)

    def __mul__(self, other):
        assert isinstance(other, Meromorphic)
        vs = self.vs
        assert len(vs) == 1
        v = vs[0]
        p = self.p.substitute({v:other.p})
        return Meromorphic(p)



def poly(arg):
    p = Poly(arg, ring)
    p = Fraction(p)
    return p


def test():

    zero = poly({})
    one = poly({():1})
    x = poly("x")
    y = poly("y")
    z = poly("z")

    assert x*y == y*x
    assert x*(y+z) == x*y + x*z

    assert str(x**4 + y) == "y + x^4"

    R = (x+y)**4

    S = (x+y)**2

    assert S**2 == R

    #assert str( 10*x / (2*x) ) == "5" 

    assert ( R / S ) == S

    assert R == x**4 + 4*x**3*y + 6*x**2*y**2 + 4*x*y**3 + y**4

    items = {"x":2, "y":x}
    assert R.substitute(items) == 16 + 32*x + 24*x**2 + 8*x**3 + x**4

    assert ( (x+I)**2 ) == -1 + 2*I*x + x**2

    # Clifford group

    i = Meromorphic( z )
    H = Meromorphic( (z+1)/(z-1) )
    S = Meromorphic( -I*z )
    Z = Meromorphic( -z )

    assert S*S == Z

    SH = S*H
    print(S*H)
    print(S*H*S)
    print(S*H*S*H)

    return

    G = mulclose([S, H], verbose=True, maxsize=100)
    print(len(G))
    for g in G:
        print(g)





if __name__ == "__main__":

    from time import time

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


    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")


