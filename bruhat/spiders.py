#!/usr/bin/env python3

from operator import add, mul, matmul
from functools import reduce
from time import time
start_time = time()

import numpy

from bruhat.lin import Space, Lin, element
from bruhat.argv import argv

def tpow(v, n, one=None):
    return reduce(matmul, [v]*n) if n>0 else one


class SMC(object):
    def __init__(self, ring, n):
        self.K = Space(ring, 1, 0, "K")
        self.V = Space(ring, n, 0, "V")
        self.one = Lin(self.K, self.K, [[1]])
        self.ring = ring
        self.n = n
    
    def ket(self, i):
        ring = self.ring
        K, V = self.K, self.V
        v = [ring.zero]*self.n
        v[i] = ring.one
        v = numpy.array(v)
        v.shape = (V.n, K.n)
        return Lin(V, K, v)
    
    def bra(self, i):
        v = self.ket(i)
        v = v.transpose()
        return v

    def green(self, m, n):
        # m legs <--- n legs
        K, V = self.K, self.V
        mid = Lin(tpow(K, m, K), tpow(K, n, K), [[1]])
        G = reduce(add, [
            tpow(self.ket(i), m, self.one) * mid * tpow(self.bra(i), n, self.one)
            for i in range(self.n)])
        return G
    
    def red(self, m, n):
        # m legs <--- n legs
        # hmmmmm...
        return R


def main():

    dim = 2
    ring = element.Z
    c = SMC(ring, dim)

    ket, bra = c.ket, c.bra
    green, red = c.green, c.red
    K, V = c.K, c.V
    one = c.one

    v0 = ket(0) # |0>
    v1 = ket(1) #    <1|
    u0 = bra(0) # |0>
    u1 = bra(1) #    <1|

    X = Lin(V, V, [[0, 1], [1, 0]])
    Z = Lin(V, V, [[1, 0], [0, -1]])

    I = v0*u0 + v1*u1
    assert green(1, 1) == I

    assert tpow(ket(0), 1) == v0
    assert tpow(ket(0), 2) == v0@v0
    assert green(2, 2) == v0 @ v0 @ u0 @ u0 + v1 @ v1 @ u1 @ u1 
    a = ket(0) * bra(0)
    assert green(1,2).shape == (v0@u0@u0).shape

    mul = green(1, 2)
    comul = green(2, 1)
    unit = green(1, 0)
    counit = green(0, 1)

    assert mul * (I @ unit) == I
    assert mul * (I @ mul) == mul * (mul @ I)
    assert mul * comul == I

    cap = counit * mul
    cup = comul * unit

    assert (I @ cap)*(cup @ I) == I

    assert ( mul * (v0 @ v0) ) == v0
    assert ( mul * (v1 @ v1) ) == v1
    assert ( mul * (v0 @ v1) ).is_zero()
    assert ( mul * (v1 @ v0) ).is_zero()

    A = counit * mul * (X @ I) * comul
    assert A.is_zero()

    assert mul * (X@X) == X*mul
    assert (X@X) * comul == comul*X
    assert counit * X == counit
    assert X * unit == unit

    print( green(2, 1) * green(1, 2) )
    print( (X@I) * green(2, 1) * green(1, 2) * (X@I) )

    return

    _mul = red(1, 2)
    _comul = red(2, 1)
    _unit = red(1, 0)
    _counit = red(0, 1)



if __name__ == "__main__":
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed:", _seed)
        numpy.random.seed(_seed)

    fn = argv.next() or "main"

    print("%s()"%fn)

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))





