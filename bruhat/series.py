#!/usr/bin/env python3

from bruhat.element import Q
from bruhat.util import factorial, cross


ring = Q


class Series(object):

    def __init__(self):
        self._cache = {}

    def __str__(self):
        return "%s(%s, ...)" % (
            self.__class__.__name__,
            ', '.join(str(value) for value in [self[i] for i in range(5)]))

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
    def getitem(self, idx):
        assert 0, "TODO"


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

def main():

    zero = Zero()
    one = One()
    x = X()

    a = 2*x+1
    assert (a**3).eq( a*a*a )
    assert (a**4).eq( a*a*a*a )

    f = Exp()
    print(f)

    




if __name__ == "__main__":

    main()


        


