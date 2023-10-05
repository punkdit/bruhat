#!/usr/bin/env python3
"""
Use sage's CyclotomicField's and Matrix's to
build spiders with phase a multiples of pi/4.
"""

from operator import add, mul, matmul
from functools import reduce

from bruhat.clifford_sage import Matrix, CyclotomicField, Clifford
from bruhat.argv import argv

dim = 2 # qubits
K = CyclotomicField(8)
w8 = K.gen()
w4 = w8**2
r2 = w8 + w8**7
assert r2**2 == 2

matrix = lambda rows : Matrix(K, rows)
I = matrix([[1, 0], [0, 1]])
H = (r2/2) * matrix([[1, 1], [1, -1]])
X = matrix([[0, 1], [1, 0]])
Z = matrix([[1, 0], [0, -1]])
Y = matrix([[0, -w4], [w4, 0]])


def parse(s):
    assert len(s)
    ops = [{"I":I, "X":X, "Z":Z, }[c] for c in s]
    return reduce(matmul, ops)


def tpow(v, n, one=matrix([[1]])):
    return reduce(matmul, [v]*n) if n>0 else one

def ket(i):
    v = [0]*dim
    v[i] = 1
    v = matrix([[w] for w in v])
    return v

def bra(i):
    v = ket(i)
    v = v.transpose()
    return v

def green(m, n, phase=0):
    # m legs <--- n legs
    assert phase in [0, 1, 2, 3]
    r = w4**phase
    G = reduce(add, [
        (r**i) * tpow(ket(i), m) * tpow(bra(i), n)
        for i in range(dim)])
    return G

def red(m, n, phase=0):
    # m legs <--- n legs
    assert phase in [0, 1, 2, 3]
    G = green(m, n, phase)
    R = tpow(H, m) * G * tpow(H, n)
    return R


def test_spider():

    v0 = ket(0) # |0>
    v1 = ket(1) #    |1>
    u0 = bra(0) # <0|
    u1 = bra(1) #    <1|

    I = v0*u0 + v1*u1
    assert green(1, 1) == I
    assert red(1, 1) == I

    assert tpow(ket(0), 1) == v0
    assert tpow(ket(0), 2) == v0@v0
    assert green(2, 2) == v0 @ v0 @ u0 @ u0 + v1 @ v1 @ u1 @ u1 
    a = ket(0) * bra(0)
    assert green(1,2).shape == (v0@u0@u0).shape

    assert v0 == (r2/2)*red(1,0)
    assert v1 == (r2/2)*red(1,0,2)

    assert green(1, 1, 2) == Z
    assert red(1, 1, 2) == X

    # green frobenius algebra
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

    assert green(2, 1) * green(1, 2) + (X@I) * green(2, 1) * green(1, 2) * (X@I) == I@I

    # red frobenius algebra
    _mul = red(1, 2)
    _comul = red(2, 1)
    _unit = red(1, 0)
    _counit = red(0, 1)

    u0 = H*ket(0) # |->
    u1 = H*ket(1) #    |+>

    assert u0 == (r2/2)*green(1,0)
    assert u1 == (r2/2)*green(1,0,2)

    assert _mul * (I @ _unit) == I
    assert _mul * (I @ _mul) == _mul * (_mul @ I)
    assert _mul * _comul == I

    cap = _counit * _mul
    cup = _comul * _unit

    assert (I @ cap)*(cup @ I) == I

    assert ( _mul * (u0 @ u0) ) == u0
    assert ( _mul * (u1 @ u1) ) == u1
    assert ( _mul * (u0 @ u1) ).is_zero()
    assert ( _mul * (u1 @ u0) ).is_zero()

    A = _counit * _mul * (Z @ I) * _comul
    assert A.is_zero()

    assert _mul * (Z@Z) == Z*_mul
    assert (Z@Z) * _comul == _comul*Z
    assert _counit * Z == _counit
    assert Z * _unit == _unit

    assert red(2, 1) * red(1, 2) + (Z@I) * red(2, 1) * red(1, 2) * (Z@I) == I@I


def test_clifford():
    n = 2
    space = Clifford(n)
    lhs = space.CX()
    rhs = (green(1, 2) @ I) * (I @ red(2, 1))
    assert lhs == (2/r2)*rhs

    lhs = space.CZ()
    rhs = (green(1, 2) @ I) * (I @ H @ I) * (I @ green(2, 1))
    assert lhs == (2/r2)*rhs

    n = 5
    space = Clifford(n)

    E = tpow(H, n) * green(n, 1)
    for i in range(n):
        E = space.CZ(i, (i+1)%n) * E

    D = E.transpose()

    assert D*E == I

    stabs = "XZZXI IXZZX XIXZZ ZXIXZ" # five qubit code
    for op in stabs.split():
        stab = parse(op)
        assert stab * E == E, op


def test_gen():

    swap = Clifford(2).SWAP()

    gen = []
    for left in [0, 1, 2]:
      for right in [0, 1, 2]:
        for phase in [0, 1, 2, 3]:
            gen.append(green(left, right, phase))
            gen.append(red(left, right, phase))

    phase = 1
    gen.append( green(1, 1, phase) @ I )
    gen.append( I @ green(1, 1, phase) )
    gen.append( red(1, 1, phase) @ I )
    gen.append( I @ red(1, 1, phase) )
    gen.append( swap )

    gen = set(gen)
    print("gen:", len(gen))

    tgt = {w4*H, H}

    found = set(gen)
    bdy = {i:[] for i in [1, 2, 4]}
    for g in gen:
        bdy[g.shape[0]].append(g)
    for _ in range(6):
        _bdy = {i:[] for i in [1, 2, 4]}
        for A in gen:
          for B in bdy[A.shape[1]]:
            C = A*B
            if C in found:
                continue
            #if C.shape == (1,1) and C[0,0]==-1:
            #    print("-1 !")
            if C in tgt:
                print("hit tgt !!!!!!!!!!")
                return
            found.add(C)
            _bdy[C.shape[0]].append(C)

        bdy = _bdy
        print("found:", len(found), [len(b) for b in bdy.values()])



def test():
    test_spider()
    test_clifford()



if __name__ == "__main__":
    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed:", _seed)
        numpy.random.seed(_seed)

    fn = argv.next() or "test"

    print("%s()"%fn)

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))





