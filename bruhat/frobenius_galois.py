#!/usr/bin/env python3
"""
investigate frobenius automorphism of finite field
"""

from functools import reduce
from operator import matmul

import numpy

from bruhat.element import FiniteField, PolynomialRing, GaloisField, CyclotomicRing
from bruhat.action import Perm, Group, mulclose
from bruhat.util import factorize, cross, all_primes, is_prime
from bruhat.argv import argv


def atpow(A, n):
    return reduce(matmul, [A]*n)


def GF(q):
    ps = factorize(q)
    assert len(set(ps)) == 1
    p = ps[0]
    r = len(ps)
    assert q==p**r
    #print(p, r)
    
    field = FiniteField(p)
    if r==1:
        return field

    ring = PolynomialRing(field)

    zero = ring.zero
    one = ring.one
    x = ring.x

    itemss = [tuple(range(p)) for i in range(r)]
    for idxs in cross(itemss):
        poly = x**r
        for i, idx in enumerate(idxs):
            poly += idx*(x**i)
        #print(poly)
        for i in range(p):
            if poly(i) == zero:
                #print("i=", i)
                break
        else:
            break
    #print("poly:", poly)
    #print([str(poly(i)) for i in range(p)])

    F = GaloisField(ring, poly)
    def frobenius(a):
        return a**p
    def hermitian(a, b):
        return a*(b**p)
    def trace_hermitian(a, b):
        return (a**p)*b + a*(b**p)
    F.frobenius = frobenius
    F.hermitian = hermitian

    one = F.one
    zero = F.zero
    def order(u):
        count = 1
        v = u
        assert v != zero
        while v != one:
            v = u*v
            count += 1
            assert count < 2*q
        return count
    for u in F.elements:
        if u==zero:
            continue
        #print(u, order(u))
        if order(u) == q-1:
            F.omega = u # generator of multiplicative group
            break
    else:
        assert 0
    return F


class Code(object):
    def __init__(self, field, gen):
        self.field = field


def test():

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
    F = GaloisField(ring, (x**2 + x + 1))

    omega = F.x
    one = F.one
    a, b = omega, omega+one
    assert a*a == b
    assert a*b == b*a
    assert a*b == one
    assert b*b == a

    assert one/omega == one+omega

    frobenius = lambda a : a**2
    assert frobenius(one) == one
    assert frobenius(one+one) == one+one
    assert frobenius(omega) == omega+1
    assert frobenius(omega+1) == omega

#    # --------------------------------------------------
#
#    for p in all_primes(20):
#        assert p>1
#        field = FiniteField(p)
#        ring = PolynomialRing(field)
#        zero = ring.zero
#        one = ring.one
#        x = ring.x
#        poly = 1+x+x**2
#        for i in range(1, p):
#            assert poly(i) != zero, "p=%d, poly=%s, i=%d"%(p, poly, i)

    # --------------------------------------------------

    p = argv.get("p", 2)
    r = argv.get("r", 2)

    print("q =", p**r)
    F = GF(p**r)
    print(F.mod)
    zero = 0

    els = F.elements
    assert len(els) == len(set(els)) == p**r

    for a in els:
      for b in els:
        if b!=0:
            c = a/b # XXX fails for p=3, r=4
            assert c*b==a

    # build the hexacode
    w = F.x
    w2 = w**2
    words = []
    for a in els:
     for b in els:
      for c in els:
        f = lambda x : a*x**2 + b*x + c
        v = [a, b, c, f(1), f(w), f(w2)]
        words.append(v)
    code = numpy.array(words)
    for word in code:
      print(' '.join("%4s"%(c) for c in word))

    print(code.shape)
    assert len(code)==4**3

    def inner(w0, w1):
        assert len(w0)==len(w1)
        n = len(w0)
        #r = sum([F.hermitian(a, b) for (a, b) in zip(w0, w1)])
        r = sum([a*F.frobenius(b) for (a, b) in zip(w0, w1)])
        return r

    for w0 in code:
      for w1 in code:
        print(inner(w0, w1), end=" ")
      print()


def test_galois():
    p = argv.get("p", 2)
    r = argv.get("r", 2)

    print("q =", p**r)
    F = GF(p**r)
    print("GF(%d) = GF(%d)[x]/(%s)" % (p**r, p, F.mod))
    print("omega =", F.omega)
    omega = F.omega
    #assert omega**2 == omega + 1

    els = F.elements
    assert len(els) == len(set(els)) == p**r

    for a in els:
      for b in els:
        if b!=0:
            c = a/b # XXX fails for p=3, r=4
            assert c*b==a

    def orbit(a):
        items = set([a])
        while 1:
            a = F.frobenius(a)
            if a in items:
                break
            items.add(a)
        return items

    lookup = {e:i for (i,e) in enumerate(els)}
    for i,e in enumerate(els):
        print("%s: %s" % (i, e))
    g = {lookup[e] : lookup[F.frobenius(e)] for e in els}
    #print(["%s:%s"%(k,v) for (k,v) in g.items()])
    values = set(g.values())
    assert len(values)==len(g)
    items = list(range(len(els)))
    g = Perm(g, items)
    I = Perm.identity(items)
    G = [I]
    h = g
    print(g)
    while h != I:
        G.append(h)
        h = g*h
        #print(len(G))
        assert len(G) <= p**r
    print("|G| =", len(G))
    G = Group.generate([g])
    print("|G| =", len(G))

    for i in els:
        n = len(orbit(i))
        if n==1:
            print("*", end="")
        print(n, end=" ")
    print()

    def inner(u, v):
        r = u*F.frobenius(v) - F.frobenius(u)*v
        return r

    print("els:", end=" ")
    for u in els:
        print(u, end=" ")
    print()
    print("inner(u,u):", end=" ")
    for u in els:
        print(inner(u,u), end=" ")
    print()
    print("inner(u,v):")
    for u in els:
      for v in els:
        print(inner(u,v), end=" ")
      print()

    print()
    for u in els:
      #print([int(inner(u,v)==0) for v in els])
      for v in els:
        i = int(inner(u,v)==0)
        print(i or '.', end=" ")
      print()
        

    print("test_galois: OK")


def test_quadratic():
    p = argv.get("p", 2)
    r = argv.get("r", 2)

    q = p**(2*r)
    F = GF(q)
    print("GF(%d) = GF(%d)[x]/(%s)" % (q, p, F.mod))
    print("omega =", F.omega)
    omega = F.omega
    #assert omega**2 == omega + 1

    frobenius = lambda a : a**(p**r)
    print("frobenius = a**%d"%(p**r,))

    print("|F| =", len(F.elements))
    print("fixed:", len([a for a in F.elements if frobenius(a)==a]))



class Matrix(object):
    def __init__(self, ring, d, A=None):
        self.ring = ring
        self.d = d
        if A is None:
            A = numpy.zeros((d, d), dtype=object)
        else:
            assert A.shape == (d, d)
            A = A.copy()
        self.A = A
        self.shape = (d, d)

    @classmethod
    def identity(cls, ring, d):
        m = Matrix(ring, d)
        for i in range(d):
            m[i, i] = ring.one
        return m

    def promote(self, x):
        if isinstance(x, Matrix):
            return x
        x = self.ring.promote(x)
        ring, d = self.ring, self.d
        m = Matrix(ring, d)
        for i in range(d):
            m[i, i] = x
        return m

    def __pow__(self, n):
        if n==0:
            return self.identity(self.ring, self.d)
        if n==1:
            return self
        A = reduce(matmul, [self.A]*n)
        return Matrix(self.ring, self.d, A)

    def __str__(self):
        vs = numpy.empty(self.shape, dtype=object)
        for idx in numpy.ndindex(self.shape):
            vs[idx] = str(self[idx])
        s = str(vs)
        s = s.replace("'", "")
        return s

    def __hash__(self):
        return hash(str(self))

    def __getitem__(self, idx):
        return self.A[idx]

    def __setitem__(self, idx, value):
        value = self.ring.promote(value)
        self.A[idx] = value

    def __eq__(self, other):
        #other = self.promote(other)
        assert isinstance(other, Matrix), other
        A = self.A == other.A
        return numpy.alltrue(A)

    def __ne__(self, other):
        return not self==other

    def __add__(self, other):
        assert isinstance(other, Matrix)
        B = self.A + other.A
        return Matrix(self.ring, self.d, B)

    def __sub__(self, other):
        assert isinstance(other, Matrix)
        B = self.A - other.A
        return Matrix(self.ring, self.d, B)

    def __mul__(self, other):
        #assert isinstance(other, Matrix)
        other = self.promote(other)
        B = numpy.dot(self.A, other.A)
        return Matrix(self.ring, self.d, B)

    def __matmul__(self, other):
        d = self.d * other.d
        A = numpy.kron(self.A, other.A)
        return Matrix(self.ring, d, A)

    def __rmul__(self, other):
        #assert isinstance(other, Matrix)
        other = self.promote(other)
        B = numpy.dot(other.A, self.A)
        return Matrix(self.ring, self.d, B)

    def transpose(self):
        A = self.A.transpose().copy()
        return Matrix(self.ring, self.d, A)



def test_hw():
    # Heisenberg-Weyl group

    # copied from dev/qudit.py
    
    d = argv.get("d", 3) # qudit dimension
    assert d>=2
    assert is_prime(d)
    
    if d==2: 
        ring = CyclotomicRing(4)
        #ring = element.Z
        gamma = ring.x
        w = gamma**2
        assert w == -1
        argv.double = True
    elif 1:
        ring = CyclotomicRing(d**2)
        gamma = ring.x
        w = gamma**d
    else:
        ring = CyclotomicRing(d)
        w = ring.x
    
    I = Matrix(ring, d)
    wI = Matrix(ring, d)
    X = Matrix(ring, d)
    Z = Matrix(ring, d)
    Zdag = Matrix(ring, d)
    S = Matrix(ring, d)
    Sdag = Matrix(ring, d)
    T = Matrix(ring, d)
    Tdag = Matrix(ring, d)

    for j in range(d):
        I[j, j] = 1
        wI[j, j] = w
        X[j, (j+1)%d] = 1
        Z[j, j] = w**j
        Zdag[j, j] = w**(d-j)

        if d==2:
            val = gamma**(j**2)
            ival = gamma**(d**2 - (j**2)%(d**2))
        else:
            val = w**(j**2)
            ival = w**(d - (j**2)%d)
        assert val*ival == 1

        S[j, j] = val
        Sdag[j, j] = ival

        if d in [2, 3, 6]:
            val = gamma**(j**3)
            ival = gamma**(d**2 - (j**3)%(d**2))
        else:
            val = w**(j**3)
            ival = w**(d - (j**3)%d)
        assert val*ival == 1

        T[j, j] = val
        Tdag[j, j] = ival

    Y = w*X*Z # ?
    Xdag = X.transpose()

    assert S*Sdag == I

    assert Z*Zdag == I
    assert X*Xdag == I

    if d>2:
        for j in range(1, d+1):
            assert (j==d) == (X**j==I)
            assert (j==d) == (Z**j==I)
    else:
        assert X != I
        assert Z != I
        assert X*X == I
        assert Z*Z == I

    assert Z*X == (w**(d-1))*X*Z

    if argv.double:
        pauli = mulclose([X, Z, gamma*I])
    else:
        pauli = mulclose([X, Z])
    print("pauli:", len(pauli))

    if d==2:
        phases = mulclose([gamma*I])
    else:
        phases = mulclose([w*I])
    print("phases:", len(phases))

    assert Zdag in pauli
    assert Xdag in pauli

    if d<=3:
        # slow..
        lhs = atpow(X, d)
        rhs = atpow(Z, d)
        assert lhs * rhs == rhs * lhs

    equ = {e : {g*e for g in phases} for e in pauli}

    G = []
    for i in range(d):
      for j in range(d):
        m = (X**i) * (Z**j)
        G.append(m)
    print(len(G))

    for g in G:
      for h in G:
        i = int(g*h == h*g)
        print(i or '.', end=" ")
      print()



if __name__ == "__main__":

    name = argv.next() or "test"
    fn = eval(name)
    fn()

