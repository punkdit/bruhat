#!/usr/bin/env python3

from operator import add, mul, matmul
from functools import reduce
from time import time
start_time = time()

from bruhat.argv import argv

from qupy.dense import Qu, Gate, bitvec
I, X, Z, H = Gate.I, Gate.X, Gate.Z, Gate.H

DIM = 2

one = Qu((), '', 1.)

def ket(*bits):
    return bitvec(*bits)

def bra(*bits):
    return ~bitvec(*bits)

def tpow(v, n):
    return reduce(matmul, [v]*n) if n>0 else one

def green(m, n):
    # m legs <--- n legs
    G = reduce(add, [tpow(ket(i), m) * tpow(bra(i), n) for i in range(DIM)])
    return G

def red(m, n):
    # m legs <--- n legs
    G = green(m, n)
    R = tpow(H, m) * G * tpow(H, n)
    return R

def main():

    assert Z*X == -X*Z

    v0 = ket(0) # |0>
    v1 = ket(1)
    assert v0 * ~v0 + v1*~v1 == I

    assert tpow(ket(0), 1) == v0
    assert tpow(ket(0), 2) == v0@v0
    assert green(2, 2) == v0 @ v0 @ ~v0 @ ~v0 + v1 @ v1 @ ~v1 @ ~v1 
    assert green(1,2).shape == (v0@v0@ ~v0).shape

    #print(red(2, 2).shortstr() )

    inv = H 

    mul = green(1, 2)
    comul = green(2, 1)
    unit = green(1, 0)
    counit = green(0, 1)

    _mul = red(1, 2)
    _comul = red(2, 1)
    _unit = red(1, 0)
    _counit = red(0, 1)

    assert mul * (I @ unit) == I
    assert mul * (I @ mul) == mul * (mul @ I)

    cap = counit * mul
    cup = comul * unit

    assert (I @ cap)*(cup @ I) == I

    assert counit * H == _counit
    assert H * unit == _unit

#    A = counit * mul * (H @ I) * comul
#    print(A)
#    A = _counit * mul * (H @ I) * comul
#    print(A)

    assert ( mul * (v0 @ v0) ) == v0
    assert ( mul * (v1 @ v1) ) == v1
    assert ( mul * (v0 @ v1) ).is_zero()
    assert ( mul * (v1 @ v0) ).is_zero()

    A = counit * mul * (X @ I) * comul
    print(A.is_zero())

    assert mul * (X@X) == X*mul
    assert (X@X) * comul == comul*X
    assert counit * X == counit
    assert X * unit == unit



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





