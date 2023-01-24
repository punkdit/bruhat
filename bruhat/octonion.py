#!/usr/bin/env python3

"""
Implement the Cayley-Dickson _construction to
get complex numbers, quaternions, octonions, sedenions, etc.

copied from sedenions.py

"""

import math, os
from functools import reduce
from operator import add, mul
from time import time
start_time = time()

from bruhat import element
from bruhat.util import choose, cross
#from bruhat.action import mulclose, Perm, Group
from bruhat.gset import mulclose, Perm, Group
from bruhat.argv import argv
from bruhat.isomorph import Point, Graph, search



class Number(object):
    def __init__(self, a):
        self.a = a
        self.shape = ()
        self.flat = (a,)

    def __str__(self):
        return str(self.a)

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.a)

    def __hash__(self):
        assert 0
        print("__hash__", self.__class__.__name__)
        return hash(str(self))

    def promote(self, a):
        if not isinstance(a, Number):
            a = Number(a) # wrap it up
        return a

    def __add__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a + other.a)

    def __sub__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a - other.a)

    def __mul__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a * other.a)

    def __rmul__(self, other):
        return Number(self.a * other)

    def __truediv__(self, other):
        return self.__class__(self.a / other)

    def __eq__(self, other):
        assert self.__class__ is other.__class__
        return self.a == other.a

    def __ne__(self, other):
        assert self.__class__ is other.__class__
        return self.a != other.a

    def __neg__(self):
        return Number(-self.a)

    def conj(self):
        return self

    def is_real(self):
        return True

    def is_zero(self):
        return self.a == 0

    def get_zero(self):
        return Number(0)


class Double(Number):
    def __init__(self, a, b):
        if not isinstance(a, Number):
            a = Number(a)
        if not isinstance(b, Number):
            b = Number(b)
        assert a.shape == b.shape
        self.shape = (a.shape, b.shape)
        assert isinstance(a, Number)
        self.a = a
        self.b = b
        self.flat = a.flat + b.flat

    def get_zero(self):
        return Double(self.a.get_zero(), self.b.get_zero())

    def promote(self, a):
        if not isinstance(a, Number):
            a = Number(a)
        if a.shape == self.shape:
            return a
        assert str(a.shape) in str(self.shape)
        a = self.a.promote(a)
        return Double(a, self.b.get_zero())

    def __repr__(self):
        #a, b = self.pair
        return "%s(%s, %s)"%(self.__class__.__name__, self.a, self.b)

    def __str__(self):
        return "(%s, %s)"%(self.a, self.b)

    def __add__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a = self.a + other.a
        b = self.b + other.b
        return self.__class__(a, b)

    def __sub__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a = self.a - other.a
        b = self.b - other.b
        return self.__class__(a, b)

    def __mul__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a, b = self.a, self.b
        c, d = other.a, other.b
        x = self.__class__(a*c - d.conj()*b, d*a + b*c.conj())
        return x

    def __rmul__(self, other):
        return self.__class__(self.a * other, self.b * other)

    def __truediv__(self, other):
        return self.__class__(self.a / other, self.b / other)

    def __eq__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        return self.a == other.a and self.b == other.b

    def __ne__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        return self.a != other.a or self.b != other.b

    def __hash__(self):
        return hash(str(self))

    def __neg__(self):
        return self.__class__(-self.a, -self.b)

    def conj(self):
        return self.__class__(self.a.conj(), -self.b)

    def norm2(self):
        return self.conj() * self

    def is_real(self):
        return self.a.is_real() and self.b.is_zero()

    def is_zero(self):
        return self.a.is_zero() and self.b.is_zero()



def is_commutative(items):
    for a in items:
      for b in items:
        if a*b != b*a:
            return False
    return True


def is_anticommutative(items):
    for a in items:
      for b in items:
        if a!=b and a*b != -b*a:
            return False
    return True


def is_associative(items):
    for a in items:
      for b in items:
        for c in items:
            if a*(b*c) != (a*b)*c:
                return False
    return True


def is_alternative(items):
    for a in items:
      for b in items:
        if a*(b*b) != (a*b)*b:
            return False
    return True

def dot(a, b): 
    ab = reduce(add, [ai*bi for (ai,bi) in zip(a.flat,b.flat)])
    return ab



def get_geometry(imag):
    graph = Graph()
    N = len(imag)
    triples = set()
    cycles = []
    for idx in range(N):
        graph.add('p')
    for idx in range(N):
      for jdx in range(N):
        if idx==jdx:
            continue
        k = imag[idx]*imag[jdx]
        if k not in imag:
            continue
        kdx = imag.index(k)
        key = [idx, jdx, kdx]
        cycle = list(key)
        key.sort()
        key = tuple(key)
        if key in triples:
            continue
        triples.add(key)
        cycles.append(cycle)
        p = graph.add('l')
        graph.join(idx, p)
        graph.join(jdx, p)
        graph.join(kdx, p)

    return graph, cycles



def find_perms(imag):
    bag0, cycles = get_geometry(imag)
    bag1, cycles = get_geometry(imag)
    #print(cycles)

    N = 3
    struct = []
    for cycle in cycles:
        items = []
        for i in range(N):
            items.append(tuple(cycle[(i+j)%N] for j in range(N)))
        items.sort()
        #print items
        struct.append(items)
    struct.sort()

    for cycle in cycles:
        cycle = [bag0[i] for i in cycle]
        nbd = set(cycle[0].nbd)
        for point in cycle:
            nbd = nbd.intersection(point.nbd)
        assert len(nbd)==1

    #print(struct)

    perms = []
    count = 0
    total = 0
    for f in search(bag0, bag1):
        _struct = [[tuple(f[i] for i in cycle) for cycle in items] for items in struct ]
        for items in _struct:
            items.sort()
        _struct.sort()
        #print(_struct)
        if struct==_struct:
            count += 1
            #print("*")
        total += 1
        perms.append(f)

    return perms


def main():

    ring = element.Q
    one, zero = ring.one, ring.zero

    x = Number(2*one)
    y = Double(2*one, zero)
    assert x==y

    # ----------- double: complex --------------------

    one = Double(ring.one, ring.zero)
    i = Double(ring.zero, ring.one)
    zero = Double(ring.zero, ring.zero)
    assert i*i == -1

    cplex = [one, i]
    assert is_commutative(cplex)
    assert is_associative(cplex)

    # ----------- double: quaternions --------------------

    zero, one, i, j, k = [
        Double(zero, zero),
        Double(one, zero),
        Double(i, zero),
        Double(zero, one),
        Double(zero, i),
    ]

    for x in [i, j, k]:
        assert x*x == -1
        for y in [i, j, k]:
            if x==y:
                continue
            assert x*y == -y*x
    assert i*j == -j*i
    assert i*j*k == -1

    quaternions = [one, i, j, k]
    assert not is_commutative(quaternions)
    assert is_anticommutative(quaternions[1:])
    assert is_associative(quaternions)

    # ----------- double: octonions --------------------

    basis = [
        Double(one, zero), 
        Double(zero, one),
        Double(i, zero), 
        Double(j, zero), 
        Double(zero, i), 
        Double(k, zero),
        Double(zero, k),
        Double(zero, j), 
    ]
    zero = Double(zero, zero)

    one = basis[0]
    imag = basis[1:]
    assert one + one == 2*one
    for i in imag:
        assert i*i == -one

    assert not is_commutative(basis)
    assert not is_associative(basis)
    assert is_anticommutative(basis[1:])
    assert is_alternative(basis)

    i0, i1, i2, i3, i4, i5, i6 = imag
    inf = one

    assert i1*i2 == i4
    assert i2*i3 == i5
    assert i3*i4 == i6
    assert i2*i1 == -i4
    assert i3*i2 == -i5
    assert i4*i3 == -i6

    assert i2*i4 == i1
    assert i3*i5 == i2
    assert i4*i6 == i3

    assert i4*i1 == i2

    def get(desc):
        x = zero
        sign = 1
        for a in desc:
            if a=='-':
                sign = -1
            elif a=='i':
                x = x+sign*inf
                sign = +1
            else:
                a = int(a)
                x = x+sign*imag[a]
                sign = +1
        return x / 2
    assert i6 == imag[6]

    lhs = (inf + i0 + i2 + i6) / 2
    assert lhs == get("i026")
    rhs = (inf + i0 + i1 + i3) / 2
    assert rhs == get("i013")
    print(lhs*rhs)
    z = lhs*rhs
    print(get("0235"))
    for idxs in choose(list("0123456"), 4):
        #print(idxs)
        x = get(idxs)
        if x==z:
            print(idxs)
            print(x)
    assert lhs*rhs == get("0156")

    assert get("0235") == (i0 + i2 + i3 + i5)/2
    assert get("02-35") == (i0 + i2 - i3 + i5)/2

    lhs = get('i356')
    rhs = get('-i356')
    assert lhs*lhs == rhs
    #assert get('i356') * get('0463') == get('-i461') # ???

    gen = "0124 0235 0346 i450 0561 i602 i013 i356 146i 25i1 3612 4i23 5134 6245"
    gen = gen.split()
    gen = [get(a) for a in gen]
    units = mulclose(gen)
    assert len(units) == 240
    units = set(units)
    assert one in units

    for i in basis:
        assert i in units
        assert -i in units

    print()
    print(len(units))

    #for a in units:
    #  for b in units:
    #    assert a*b in units

#    # j-frame, page 138
#    frame = [
#        one,
#        get("-0124"), get("0-124"), get("01-24"), -i6,
#        get("012-4"), -i3, -i6
#    ]

    fmt = "(" + 8*"%4s," + ")"

    def send(x):
        xs = x.flat
        items = [xi*j for (xi,j) in zip(xs, frame)]
        return reduce(add, items)

    def act(perm, x):
        xs = x.flat
        xs = [xs[perm[i]+1] for i in range(7)]
        result = reduce(add, [xi*yi for (xi,yi) in zip(xs, imag)])
        result += x.flat[0]*one
        return result

    perms = find_perms(imag)
    for perm in perms:
        auto = True
        for a in basis:
            a1 = act(perm, a)
            assert a1 in basis
            for b in basis:
                b1 = act(perm, b)
                if a1*b1 != act(perm, a*b):
                    print(a1*b1)
                    print(act(perm, a*b))
                    print()
                    auto = False
                    break
        if auto:
            print("found")

    return

    found = []
    for a in units:
        if dot(a, one) == 0:
            found.append(a)
    assert len(found) == 126

    found = []
    for a in units:
        if a==one:
            continue
        if a*a==one:
            continue
        if a*a*a==one:
            found.append(a)
    print(len(found))

    perms = []
    lookup = {unit:idx for (idx, unit) in enumerate(units)}
    n = len(units)
    for x in found:
        xc = x.conj()
        perm = [None]*n
        for a in units:
            b = xc*a*x
            perm[lookup[a]] = lookup[b]
        assert None not in perm
        perm = Perm(perm)
        perms.append(perm)
    #G = Group(None, perms)
    G = mulclose(perms, verbose=True)
    print(len(G))

if __name__ == "__main__":

    from argv import argv
    name = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%name)
    else:
        fn = eval(name)
        fn()

    print("OK: %.3f seconds\n"%(time() - start_time))






