#!/usr/bin/env python

"""

Here we implement chain complexes using Kapranov-Voevodsky 2-vector spaces.

A chain:
    C0 <--- C1 <--- C2 ...

is encoded as a "circulant" 2-matrix:

    C0  C1  C2               .   C0  C1 
    .   C0  C1   --------->  .   .   C0 
    .   .   C0               .   .   .

with components the boundary map.

"""


from time import time
start_time = time()



from random import choice
from functools import reduce

import numpy

from bruhat import elim
from bruhat.lin import Lin, Space, AddSpace, element
from bruhat.smap import SMap
from bruhat.argv import argv
from bruhat.vect2 import Rig, Cell0, Cell1, Cell2


def get_shift(rig, cell0):
    N, K = rig.zero, rig.one
    X = [[(K if i-1==j else N) for i in cell0] for j in cell0]
    X = Cell1(cell0, cell0, X)
    return X


class KVChain(object):

    def __init__(self, bdy):
        assert isinstance(bdy, Cell2)
        n = bdy.src.src
        assert bdy.src.tgt == n
        self.bdy = bdy
        self.rig = bdy.rig
        self.n = n
        self.X = get_shift(self.rig, self.n)
        self.check()

    def check(self):
        bdy = self.bdy
        m, n = bdy.shape
        assert m==n
        tgt, src = bdy.hom
        for i in range(n-1):
            for j in range(n):
                assert bdy[j,i+1].tgt == bdy[j,i].src
                assert (bdy[j,i] * bdy[j,i+1]).is_zero()
        assert bdy.tgt == (self.X<<bdy.src).normalized

    @classmethod
    def from_lins(cls, rig, lins):
        prev = None
        for lin in lins:
            assert isinstance(lin, Lin)
            tgt, src = lin.hom
            if prev is not None:
                assert prev.tgt == src
                assert (lin*prev).is_zero()
            prev = lin

        grades = [] if not lins else [lin.src for lin in lins] + [lins[-1].tgt]
        grades.sort(key = lambda space:space.grade)
        #print("grades:", grades)
        bdys = list(reversed(lins))
        n = len(grades)+1
        # f : B <--- A
        A = numpy.empty((n, n), dtype=object)
        B = numpy.empty((n, n), dtype=object)
        f = numpy.empty((n, n), dtype=object)
        for i in range(n):
          for j in range(n):
            A[i, j] = rig.zero
            B[i, j] = rig.zero
        for i in range(n):
          for j in range(i, n):
            idx = j-i
            if idx < len(grades):
                A[i, j] = grades[idx]
            if idx>0:
                B[i, j] = grades[idx-1]
        for i in range(n):
          for j in range(n):
            if 0<=j-i-1<len(bdys):
                f[i,j] = bdys[j-i-1]
                #print("f[i,j]", i, j)
            else:
                f[i,j] = Lin.zero(B[i,j], A[i,j])
        cell = Cell0(rig, n, "n")
        src = Cell1(cell, cell, A)
        tgt = Cell1(cell, cell, B)
        bdy = Cell2(tgt, src, f)
        return KVChain(bdy)

    def __eq__(self, other):
        assert isinstance(other, KVChain)
        return self.bdy.hom == other.bdy.hom and self.bdy == other.bdy

    def identity(self):
        bdy = self.bdy
        f = Cell2.identity(bdy.src)
        return KVChainMap(self, self, f)

    mulcache = {} 
    def __matmul__(self, other): # KVMulChain ??
        assert isinstance(other, KVChain)
        assert self.rig == other.rig
        key = id(self), id(other)
        if key in KVChain.mulcache:
            return KVChain.mulcache[key]

        rig = self.rig
        bdyC = self.bdy
        XC, C = bdyC.hom
        n = C.src
        assert C.tgt == n
        assert XC.hom == (n, n)
    
        N, K = rig.zero, rig.one
        X = [[(K if i-1==j else N) for i in n] for j in n]
        X = Cell1(n, n, X)
        assert XC == (X<<C).normalized
        assert XC == (C<<X).normalized # commutes !
    
        bdyD = other.bdy
        XD, D = bdyD.hom
        assert D.hom == (n, n)
        assert XD.hom == (n, n)
    
        CD = C<<D
        CD = CD.normalized
    
        front = bdyC << Cell2.identity(D)
        front = front.normalized
        assert front.src == CD
    
        back = Cell2.identity(C) << bdyD
        back = back.normalized
        assert back.src == CD
    
        XCD = (XC << D).normalized
        assert front.tgt == XCD
        assert back.tgt == XCD
    
        mid = front + back
        left = Cell1.codiagonal(n)
        right = Cell1.diagonal(n)
        mid = left << mid << right
        mid = mid.normalized
    
        bot = Cell2.diagonal(CD)
    
        cell = mid * bot
    
        top = Cell2.codiagonal(XCD)
        assert top.tgt == XCD
    
        bdy = top * cell

        C = KVChain(bdy)
        KVChain.mulcache[key] = C
        return C


class KVChainMap(object):
    def __init__(self, tgt, src, f):
        assert isinstance(tgt, KVChain)
        assert isinstance(src, KVChain)
        assert isinstance(f, Cell2)
        assert f.tgt == tgt.bdy.src
        assert f.src == src.bdy.src
        self.hom = (tgt, src)
        self.tgt = tgt
        self.src = src
        self.f = f
        self.check()

    def check(self):
        tgt, src = self.hom
        f = self.f
        X = src.X
        lhs = (X<<f) * src.bdy
        rhs = tgt.bdy * f
        assert lhs.normalized == rhs.normalized

    @classmethod
    def zero(cls, tgt, src):
        assert isinstance(src, KVChain)
        assert isinstance(tgt, KVChain)
        f = Cell2.zero(tgt.bdy.src, src.bdy.src)
        return KVChainMap(tgt, src, f)

    @classmethod
    def get_symmetry(cls, C, D):
        tgt = D@C
        src = C@D
        #print(tgt.bdy.src)
        #print(src.bdy.src)
        f = Cell2.zero(tgt.bdy.src, src.bdy.src) # XXX todo
        return KVChainMap(tgt, src, f)

    def __eq__(self, other):
        assert isinstance(other, KVChainMap)
        return self.f == other.f

    def __mul__(self, other):
        assert isinstance(other, KVChainMap)
        assert self.src == other.tgt
        f = self.f * other.f
        return KVChainMap(self.tgt, other.src, f)



def test():

    p = 2
    ring = element.FiniteField(p)
    rig = Rig(ring)

    U = Space(ring, 4, 1, "U")
    V = Space(ring, 4, 0, "V")

    A = elim.parse(ring, "11.. .11. ..11 1..1")
    f = Lin(V, U, A)

    C = KVChain.from_lins(rig, [f])

    i = C.identity()
    assert i*i == i

    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 3, 0, "C_0")
    C = KVChain.from_lins(rig, [Lin(C0, C1, [[1,1,0],[0,1,1],[1,0,1]])])

    #CC = C@C


    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 2, 0, "C_0")
    C = KVChain.from_lins(rig, [Lin(C0, C1, [[1,1,0],[0,1,1]])])

    D1 = Space(ring, 4, 1, "D_1")
    D0 = Space(ring, 3, 0, "D_0")
    D = KVChain.from_lins(rig, [Lin(D0, D1, [[1,1,0,0],[0,1,1,0],[0,0,1,1]])])

#    CD = C@D
#    DC = D@C
#    assert CD != DC

    C2 = Space(ring, 1, 2, "C_2")
    C1 = Space(ring, 1, 1, "C_1")
    C0 = Space(ring, 1, 0, "C_0")
    C = KVChain.from_lins(rig, [Lin(C1, C2, [[0]]), Lin(C0, C1, [[0]])])

    D2 = Space(ring, 1, 2, "D_2")
    D1 = Space(ring, 1, 1, "D_1")
    D0 = Space(ring, 1, 0, "D_0")
    D = KVChain.from_lins(rig, [Lin(D1, D2, [[0]]), Lin(D0, D1, [[0]])])

    E2 = Space(ring, 1, 2, "E_2")
    E1 = Space(ring, 1, 1, "E_1")
    E0 = Space(ring, 1, 0, "E_0")
    E = KVChain.from_lins(rig, [Lin(E1, E2, [[0]]), Lin(E0, E1, [[0]])])

    CDE = (C@D)@E
    CDE = C@(D@E)
    bdy = CDE.bdy
    src = bdy.src
    print(src[0,2].name)

    #s = KVChainMap.get_symmetry(D, E)

    return

    assert s.src == CD
    assert s.tgt == DC

    si = KVChainMap.get_symmetry(D, C)
#    assert si*s == CD.identity()
#    assert s*si == DC.identity()
    assert si*s == KVChainMap.zero(CD, CD)
    assert s*si == KVChainMap.zero(DC, DC)


if __name__ == "__main__":


    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))



