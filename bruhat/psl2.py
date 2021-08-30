#!/usr/bin/env python

"""
PSL(2,Q)

"""

from random import randint

import numpy

from bruhat.argv import argv
if argv.fast:
    from bruhat import _element as element
else:
    from bruhat import element
from bruhat import elim

#from bruhat.chain import Space, Lin

"""
Aut(H) = {
z |->  az + b
       ------
       cz + d
}
= PSL(2, R)

"""

class Mat(object):
    def __init__(self, ring, a, b, c, d):
        self.ring = ring
        pos = (a, b, c, d)
        neg = (-a, -b, -c, -d)

        # we need to make this canonical for hash( )
        zero = ring.zero
        if a == zero and b > zero:
            pass # OK
        elif a == zero and b < zero:
            pos, neg = neg, pos # swap
        elif a == zero:
            assert False, (a, b)
        elif a < zero:
            pos, neg = neg, pos # swap
        self.pos = pos
        self.neg = neg

    @classmethod
    def construct(cls, ring, a, b, c, d):
        a = ring.promote(a)
        b = ring.promote(b)
        c = ring.promote(c)
        d = ring.promote(d)
        x = a*d - b*c
        assert x == ring.one
        return cls(ring, a, b, c, d)

    @property
    def a(self):
        return self.pos[0]

    @property
    def b(self):
        return self.pos[1]

    @property
    def c(self):
        return self.pos[2]

    @property
    def d(self):
        return self.pos[3]


    def check(self):
        a, b, c, d = self.pos
        x = a*d - b*c
        assert x==self.ring.one

    def inv(self):
        ring = self.ring
        a, b, c, d = self.pos
        m = Mat(ring, d, -b, -c, a)
        m.check()
        return m

# um... need complex numbers here...
#    def send(self, z):
#        "Apply mobius transform"
#        z = self.ring.promote(z)
#        a, b, c, d = self.pos
#        z = (a*z + b) / (c*z + d)
#        return z
#    __call__ = send

    @classmethod
    def rand(cls, ring, n=5):
        assert n>0
        while 1:
            a = ring.promote(randint(-n, n))
            b = ring.promote(randint(-n, n))
            if a != ring.zero or b != ring.zero: # easy
                break
        if a != ring.zero:
            c = ring.promote(randint(-n, n))
            d = ( ring.one + b*c ) / a
        else:
            d = ring.promote(randint(-n, n))
            c = (a*d - 1)/b
        m = Mat(ring, a, b, c, d)
        m.check()
        return m

    def __mul__(self, other):
        assert isinstance(other, Mat)
        assert self.ring is other.ring
        a, b, c, d = self.pos
        e, f, g, h = other.pos
        m = Mat(self.ring, a*e + b*g, a*f+b*h, c*e+d*g, c*f+d*h)
        return m

    def __pow__(self, i):
        assert i>=0
        assert int(i)==i
        if i==0:
            return self.construct(1, 0, 0, 1)
        if i==1:
            return self
        A = self
        while i>1:
            A = A*self
            i -= 1
        return A

    def __eq__(self, other):
        assert isinstance(other, Mat)
        assert self.ring is other.ring
        return self.pos == other.pos or self.pos == other.neg

    def __hash__(self):
        return hash(self.pos)

    def __str__(self):
        return "Mat(%s, %s, %s, %s)"%(self.pos)
    __repr__ = __str__

    def is_modular(self, n):
        a, b, c, d = self.pos
        result = (a.top%n==a.bot%n)
        result = result and (d.top%n==d.bot%n)
        result = result and (d.top%n==d.bot%n)

    is_modular = lambda m : m.a%n==1 and m.d%n==1 and m.b%n==0 and m.c%n==0

def mulclose_fast(gen, verbose=False, maxsize=None):
    els = set(gen)
    bdy = list(els)
    changed = True
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B
                if C not in els:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return els


def mulclose_subgroup(ring, gen, test, verbose=False, maxsize=None):
    "test is a callback: is the element in the subgroup ?"
    els = set(g for g in gen if not test(g))
    bdy = list(els)
    changed = True
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        # check loop invariant
        for a in els:
          for b in els:
            if a is b:
                continue
            assert not test(a.inv()*b)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B
                if C in els or test(C):
                    continue
                Ci = C.inv()
                for D in els:
                    if test(Ci*D): # represents same coset
                        break
                else:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    els.add(Mat.construct(ring, 1, 0, 0, 1))
    return els




def test_2q():

    # --------------------------------------------------------------
    # PSL(2, Q)

    ring = element.Q
    construct = lambda *args : Mat.construct(ring, *args)
    I = construct(1, 0, 0, 1)
    nI = construct(-1, 0, 0, -1)
    A = construct(1, 1, 0, 1)
    Ai = construct(1, -1, 0, 1)
    #print(construct(2, 0, 0, 2))
    #assert I == construct(2, 0, 0, 2)

    assert I == nI

    assert A == A
    assert A != I
    assert A*A == construct(1, 2, 0, 1)
    assert A*Ai == I

    for i in range(1000):
        A = Mat.rand(ring, 3)
        B = Mat.rand(ring, 3)
        C = A*B

        lhs = (A==B) 
        rhs = (hash(A) == hash(B))
        assert not lhs or rhs

        assert A.inv()*A == I


def test_psl2z():

    # --------------------------------------------------------------
    # PSL(2, Z)

    # https://en.wikipedia.org/wiki/Modular_group
    # G = <S, T | S*S==I, (S*T)**3 == I>

    ring = element.Z
    construct = lambda *args : Mat.construct(ring, *args)

    I = construct(1, 0, 0, 1)
    S = construct(0, -1, 1, 0)
    T = construct(1, 1, 0, 1)
    assert S*S == I
    assert (S*T)**3 == I

    N = 100
    G = mulclose_fast([S, T], maxsize=N)
    print(len(G))
    assert len(G) >= N

    n = argv.get("n", 3)
    is_modular = lambda m : (
        (m.a%n==1 and m.d%n==1 or (-m.a)%n==1 and (-m.d)%n==1)
        and m.b%n==0 and m.c%n==0
    )

    for g in G:
      assert is_modular(g) == is_modular(g.inv())
      for h in G:
        if is_modular(g) and is_modular(h):
            #print()
            #print(g)
            #print(h)
            #print(g*h)
            assert is_modular(g*h)
            assert is_modular(h*g)

    J = mulclose_subgroup(ring, [S, T], is_modular, verbose=True)
    print(len(J))

    #for m in J:
    #    print(m)

def astr(A):
    s = str(A)
    s = s.replace(" 0", " .")
    return s

def test_bring():
    "Bring's curve"

    faces = list(range(1,13))
    pairs = (
        "1,2 1,3  1,4  1,5  1,6 "
        "2,4 2,11 2,1  2,10 2,5 "
        "3,7 3,6  3,5  3,8  3,1 "
        "4,9 4,2  4,6  4,10 4,1 "
        "5,7 5,1  5,11 5,3  5,2 "
        "6,8 6,4  6,3  6,9  6,1 "
        "7,8 7,5  7,12 7,3  7,11 "
        "8,9 8,3  8,12 8,6  8,7 "
        "9,10 9,6  9,12 9,4 9,8 "
        "10,11 10,4 10,12 10,2 10,9 "
        "11,12 11,5 11,10 11,7 11,2 "
        "12,11 12,9 12,7 12,10 12,8 "
    ).split()
    pairs = [eval(item) for item in pairs]
    assert len(pairs) == len(set(pairs))

    edges = set()
    for (a,b) in pairs:
        if a > b:
            a, b = b, a
        if (a,b) in edges:
            continue
        edges.add((a,b))
    #print("edges:", len(edges))
    edges = list(edges)
    edges.sort()
    #print(edges)
    elookup = dict((e, i) for (i, e) in enumerate(edges))

    n = len(edges)
    mz = len(faces)
    Hz = numpy.zeros((n, mz), dtype=int)
    for (a,b) in pairs:
        j = a-1
        sign = 1
        if a > b:
            sign = -1
            a, b = b, a
        i = elookup[(a, b)]
        assert Hz[i, j] == 0
        Hz[i, j] = sign

    #print(astr(Hz))

    for i in range(n):
        assert Hz[i].sum() == 0, ("lost edge: %s" % edges[i])

    #print("bdy:")
    bdy = dict((i,[]) for i in faces)
    for (a,b) in pairs:
        bdy[a].append(b)
    #print(bdy)

    def clockwise(face, other):
        "move clockwise around face starting at other face"
        assert face in faces
        assert other in faces
        (a, b) = face, other
        if a>b:
            a, b = b, a
        assert (a,b) in edges
        others = bdy[face]
        i = others.index(other)
        nxt = others[(i-1)%len(others)]
        return nxt

    assert clockwise(1, 2) == 6
    assert clockwise(12, 8) == 10
    
    # each vert will be a tuple of face's
    VERT = 5
    verts = []
    for face in faces:
        for other in bdy[face]:
            vert = [face]
            while other != face:
                vert.append(other)
                other = clockwise(vert[-1], vert[-2])
            assert len(vert) == VERT
            i = vert.index(min(vert))
            vert = tuple(vert[(i+j)%VERT] for j in range(VERT))
            verts.append(vert)
    verts = set(verts)

    mx = len(verts)
    Hx = numpy.zeros((n, mx), dtype=int)
    for row, vert in enumerate(verts):
        for i in range(VERT):
            a, b = vert[i], vert[(i+1)%VERT]
            sign = 1
            if a > b:
                sign = -1
                a, b = b, a
            col = elookup[a, b]
            Hx[col, row] = sign

    #print(astr(Hx))

    A = numpy.dot(Hx.transpose(), Hz)
    assert numpy.alltrue(A==0)


#
#
#    F2 = element.FiniteField(2)

    
    
    


if __name__ == "__main__":
    test_bring()







