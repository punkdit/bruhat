#!/usr/bin/env python3

"""
Build Jacobi-Mumford theta series.

"""

from bruhat.smap import SMap


class Ring(object):
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return name

Z = Ring("Z")
Z.one = 1
Z.zero = 0


class PolyRing(Ring):
    def __init__(self, base, variables):
        assert isinstance(base, Ring)
        self.base = base
        self.variables = list(variables)
        rank = len(variables)
        self.rank = rank # number of variables
        self.zero = FiniteSeries(self, {})
        one = base.one
        self.one = FiniteSeries(self, {(0,)*rank : one})
        for i in range(rank):
            v = variables[i]
            cs = [0]*rank
            cs[i] = 1
            s = FiniteSeries(self, {tuple(cs) : one})
            setattr(self, v, s)

    def __str__(self):
        return "%s[%s]"%(self.base.name, ','.join(self.variables))


class Series(object):

    def __init__(self, ring):
        self.ring = ring
        self.rank = ring.rank # number of variables


class FiniteSeries(Series):
    def __init__(self, ring, _coeffs):
        Series.__init__(self, ring)
        coeffs = {}
        for k,v in _coeffs.items():
            assert len(k) == self.rank
            if v == 0:
                continue
            coeffs[k] = v
        self.coeffs = coeffs

    def __len__(self):
        return len(self.coeffs)

    def __str__(self):
        coeffs = self.coeffs
        keys = list(coeffs.keys())
        def fn(key):
            k = tuple(abs(k) for k in reversed(key))
            return (sum(k), k)
        keys.sort(key=fn)
        vs = self.ring.variables
        base = self.ring.base
        items = []
        for key in keys:
            assert len(key) == len(vs)
            term = ["%s**%s"%(b,e)+" " for (b,e) in zip(vs, key) if e!=0]
            term = ''.join(term)
            term = term.replace("**1 ", " ")
            term = term.replace(" ", "")
            v = coeffs[key]
            assert v != base.zero
            if v != base.one or not term:
                term = "%s%s"%(v, term)
            items.append(term)
        expr = " + ".join(items) or "0"
        return expr

    def smap(self):
        smap = SMap()
        assert self.rank == 2
        kmin, jmin = 0, 0
        kmax, jmax = 0, 0
        coeffs = self.coeffs
        for (k,j) in coeffs.keys():
            kmin = min(k, kmin)
            jmin = min(j, jmin)
            kmax = max(k, kmax)
            jmax = max(j, jmax)
        for (key,v) in coeffs.items():
            k,j = key
            smap[j-jmin+1,k-kmin] = str(v)
        for k in range(kmin, kmax):
            smap[0,k-kmin] = '-' if k!=0 else "*"
        return smap

    def __eq__(self, other):
        ring = self.ring
        if type(other) is int:
            other = FiniteSeries(ring, {(0,)*ring.rank : other})
        assert isinstance(other, FiniteSeries)
        assert ring == other.ring
        return self.coeffs == other.coeffs

    def __add__(self, other):
        ring = self.ring
        if type(other) is int:
            other = FiniteSeries(ring, {(0,)*ring.rank : other})
        assert isinstance(other, FiniteSeries)
        assert ring == other.ring
        zero = ring.base.zero
        lhs, rhs = self.coeffs, other.coeffs
        coeffs = dict(lhs)
        for k,v in rhs.items():
            coeffs[k] = coeffs.get(k, zero) + v
        return FiniteSeries(ring, coeffs)

    __radd__ = __add__

    def __rmul__(self, other):
        assert type(other) is int
        coeffs = dict((k, other*v) for (k,v) in self.coeffs.items())
        return FiniteSeries(self.ring, coeffs)

    def __mul__(self, other):
        assert self.ring == other.ring
        ring = self.ring
        if type(other) is int:
            other = FiniteSeries(ring, {(0,)*ring.rank : other})
        assert isinstance(other, FiniteSeries)
        assert ring == other.ring
        zero = ring.base.zero
        lhs, rhs = self.coeffs, other.coeffs
        coeffs = {}
        for (lk,lv) in lhs.items():
          for (rk,rv) in rhs.items():
            k = tuple(l+r for (l,r) in zip(lk,rk))
            coeffs[k] = coeffs.get(k, zero) + lv*rv
        return FiniteSeries(ring, coeffs)

    def __pow__(self, n):
        if n<0:
            if len(self) != 1:
                assert 0, "can only invert monomials"
            k,v = list(self.coeffs.items())[0]
            k = tuple(n*e for e in k)
            a = FiniteSeries(self.ring, {k:v})
        else:
            a = self.ring.one
            for i in range(n):
                a = self*a
        return a


def test():

    ring = PolyRing(Z, 'pq')
    one, zero = ring.one, ring.zero
    p, q = ring.p, ring.q

    M = Series(ring)

    assert str(ring.zero) == "0", str(ring.zero)
    assert str(ring.one) == "1", str(ring.one)
    assert str(ring.p) == "p"
    assert str(ring.q) == "q"

    #assert str(p+q) == "q + p"
    #assert str(p+q+7) == "7 + q + p"

    assert (p+q)*(p+q) == q**2+2*p*q+p**2
    assert (p+q)*(p+q) == (p+q)**2
    #assert str((p+q)**2) == "q**2 + 2pq + p**2", str((p+q)**2) 

    assert (q**-2) * q * q == 1


    JM10 = 1
    for j in range(1, 10):
        k = j * (j-1) // 2
        JM10 += p**j * q**k
        k = j * (j+1) // 2
        JM10 += p**-j * q**k

    print(JM10.smap())

    JM10_2 = JM10*JM10
    print(JM10_2.smap())



if __name__ == "__main__":

    test()

    print("OK")



