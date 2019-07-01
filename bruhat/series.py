#!/usr/bin/env python3

from bruhat.element import Q
from bruhat.util import factorial, cross


ring = Q


class Series(object):

    def __init__(self):
        self._cache = {}

    def __str__(self):
        return "Series(%s, ...)" % (
            ', '.join(str(value) for value in [self[i] for i in range(6)]))

    def getitem(self, idx):
        return ring.zero

    def __getitem__(self, idx):
        assert isinstance(idx, int)
        assert idx >= 0
        value = self._cache.get(idx)
        if value is None:
            value = self.getitem(idx)
            self._cache[idx] = value
        return value

    @classmethod
    def promote(cls, item):
        if isinstance(item, Series):
            return item
        return RMul(item, One())

    def __add__(a, b):
        b = Series.promote(b)
        return Add(a, b)

    def __neg__(a):
        return RMul(-ring.one, a)

    def __sub__(a, b):
        b = Series.promote(b)
        return Add(a, RMul(-ring.one, b))

    def __rmul__(a, r):
        return RMul(r, a)

    def __mul__(a, b):
        b = Series.promote(b)
        return Mul(a, b)

    def __pow__(a, n):
        assert isinstance(n, int)
        if n==0:
            return One() # what about 0**0 ?
        if n==1:
            return a
        if n==2:
            return a*a
        return Pow(a, n)

    def __call__(a, b):
        b = Series.promote(b)
        return Compose(a, b)

    def eq(a, b, idx=10):
        b = Series.promote(b)
        for i in range(idx):
            if a[i] != b[i]:
                return False
        return True


class RMul(Series):
    def __init__(self, r, child):
        assert isinstance(child, Series)
        Series.__init__(self)
        self.r = r
        self.child = child

    def getitem(self, idx):
        return self.r * self.child[idx]


class Binop(Series):
    def __init__(self, a, b):
        assert isinstance(a, Series)
        assert isinstance(b, Series)
        Series.__init__(self)
        self.a = a
        self.b = b


class Add(Binop):
    def getitem(self, idx):
        return self.a[idx] + self.b[idx]


class Sub(Binop):
    def getitem(self, idx):
        return self.a[idx] - self.b[idx]


class Mul(Binop):
    def getitem(self, idx):
        a, b = self.a, self.b
        value = ring.zero
        for j in range(idx+1):
            value += a[j] * b[idx-j]
        return value


class Pow(Series):
    def __init__(self, child, n):
        assert isinstance(child, Series)
        Series.__init__(self)
        self.n = n
        assert n>2
        self.child = child

    def getitem(self, idx):
        idxs = list(range(idx+1))
        n = self.n
        child = self.child
        value = ring.zero
        for jdxs in cross([idxs]*(n-1)):
            total = 0
            walue = ring.one
            for j in jdxs:
                walue *= child[j]
                total += j
            if total <= idx:
                walue *= child[idx-total]
                value += walue
        return value


class Compose(Binop):
    def __init__(self, a, b):
        assert b[0] == ring.zero
        Binop.__init__(self, a, b)

    def getitem(self, idx):
        a, b = self.a, self.b
        value = ring.zero
        for i in range(idx+1):
            value += a[i] * ((b**i)[idx])
        return value


# ----------------------------------------

class Zero(Series):
    pass

class One(Series):
    def getitem(self, idx):
        if idx == 0:
            return ring.one
        return ring.zero

class X(Series):
    def getitem(self, idx):
        if idx == 1:
            return ring.one
        return ring.zero


class Exp(Series):
    def getitem(self, idx):
        return ring.one / factorial(idx)


# ----------------------------------------

#class F(Series):
#    def getitem(self, idx):
        

def main():

    zero = Zero()
    one = One()
    x = X()

    assert ((x+1)**2).eq(x**2 + 2*x + 1)

    a = 2*x+1
    assert (a**3).eq( a*a*a )
    assert (a**4).eq( a*a*a*a )

    exp = Exp()

    #print(exp(zero))
    #print(exp(x+x))
    #print(exp(x**2))
    #print(exp(-x))

    f = exp - 1
    #f = exp - 1 + x**3 # works too...
    print("f =", f)

    # -----------------------------
    # From: http://math.ucr.edu/home/baez/trimble/trimble_lie_operad.pdf
    # Page 8.
    # find compositional inverse

    assert f[0] == 0
    assert f[1] == 1
    
    def delta(h):
        return h(f) - h

    g = zero
    d = x
    for n in range(5):
        s = (-1)**n * d
        d = delta(d)
        g += s

    print("g =", g)

    print("f(g) = ", f(g))
    print("g(f) = ", g(f))

    # ------------------------------------------------------------------------------------------
    # Verifying https://golem.ph.utexas.edu/category/2018/01/more_secrets_of_the_associahed.html

    C = f
    D = g

    c1, c2, c3, c4 = C[2], C[3], C[4], C[5]
    d1, d2, d3, d4 = D[2], D[3], D[4], D[5]
    
    assert(d1 == -c1)
    assert(d2 == -c2 + 2*c1**2)
    assert(d3 == -c3 + 5*c2*c1 - 5*c1**3)
    assert(d4 == -c4 + 6*c3*c1 + 3*c2**2 - 21*c2*c1**2 + 14*c1**4)
    
    




if __name__ == "__main__":

    main()


        


