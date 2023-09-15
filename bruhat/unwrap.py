#!/usr/bin/env python

from time import time
start_time = time()
from random import shuffle, choice, randint
from functools import reduce
from operator import add

import numpy

from bruhat.action import Perm

from qumba.solve import (parse, shortstr, linear_independent, eq2, dot2, identity2,
    rank, rand2, pseudo_inverse, kernel, direct_sum, zeros2, solve2, normal_form)
from qumba.qcode import QCode, SymplecticSpace
from qumba import construct
from qumba.autos import get_isos
from qumba.csscode import find_zx_duality, find_autos
from qumba.argv import argv


def unwrap(code, check=True):
    H0 = code.deepH
    m, n, _ = H0.shape
    Sx = H0[:, :, 0]
    Sz = H0[:, :, 1]
    Sxz = numpy.concatenate((Sx, Sz), axis=1)
    Szx = numpy.concatenate((Sz, Sx), axis=1)
    H = zeros2(2*m, 2*n, 2)
    H[:m, :, 0] = Sxz
    H[m:, :, 1] = Szx
    code = QCode(H, check=check)
    return code


def zxcat(code, duality=None):
    #print(duality)
    pairs = []
    perm = []
    for (i, j) in enumerate(duality):
        if i==j:
            return None
        assert i!=j
        if i < j:
            pairs.append((i, j))
            perm.append(i)
            perm.append(j)
    assert len(pairs)*2 == len(duality)
    #print(pairs)

    right = code.apply_perm(perm)
    right = right.to_qcode()

    inner = construct.get_422()
    left = len(pairs) * inner

    #print(left)
    #print(right)
    right = QCode.trivial(left.n - right.n) + right

    code = left << right
    return code


def test_randcat():
    #css = construct.toric(2, 2)
    css = QCode.fromstr("""
.XX.XX..
X..XXX..
..ZZZ.Z.
..ZZ.Z.Z
X.X...XX
ZZ...ZZ.
    """)
    css = construct.reed_muller()
    css = construct.get_10_2_3()

    code = css.to_qcode()
    assert code.n%2 == 0
    inner = construct.get_422()

    idxs = list(range(code.n))

    for trial in range(1000):
        remain = list(idxs)
        pairs = []
        while remain:
            idx = remain.pop(randint(0, len(remain)-1))
            jdx = remain.pop(randint(0, len(remain)-1))
            pairs.append((idx, jdx))
        duality = [None]*code.n
        for (i,j) in pairs:
            duality[i] = j
            duality[j] = i
    
        #lhs = reduce(add, [inner]*(code.n//2))
        #print(lhs)
    
        dode = zxcat(code, duality)
        lhs = dode.is_selfdual()

        dode = code.apply_perm(duality)
        dode = dode.apply_H()
        rhs = dode.equiv(code)

        assert lhs == rhs
        if lhs:
            print("/", flush=True, end="")
        else:
            print(".", flush=True, end="")

    print()
    



def test_biplanar():
    code = construct.toric(4, 4)
    #code = construct.toric(2, 2, 1)
    #code = construct.biplanar(12, 6)

    print(code)
    #print(code.distance())

    items = list(range(code.n))
    perm = lambda f:Perm(dict(enumerate(f)), items)

    isos = find_autos(code.Ax, code.Az)
    isos = [perm(iso) for iso in isos]
    print(len(isos))

    duality = find_zx_duality(code.Ax, code.Az)
    duality = perm(duality)

    for iso in isos:
        d = duality * iso
        fixed = d.fixed()
        if len(fixed):
            continue
        if (d*d).is_identity():
            print("found", d)
            break
    else:
        print("not found")
        return

    pairs = [tuple(o) for o in d.orbits()]
    pairs.sort()
    print(pairs)
    n = len(pairs)
    m = code.mx
    H = zeros2(m, 2*n)
    lookup = {}

    cols = reduce(add, pairs)
    Hx = code.Hx
    print("Hx:", Hx.shape)
    H = Hx[:, cols]
    qcode = QCode(H)
    print(qcode)
    print(qcode.longstr())
    #print(qcode.get_params())

    code = code.to_qcode()

    dode = zxcat(code, [d[i] for i in range(code.n)])
    print()
    print(dode, dode.get_params())
    assert dode.is_selfdual()
    dode = dode.to_css()
    Hx, Hz = dode.Hx, dode.Hz
    print("row weights:", Hx.sum(1))


def test_zx():
    for code in QCode.load_codetables():
        if code.n > 14:
            break
        if code.k == 0:
            continue
        n = code.n
        print()
        code2 = unwrap(code)
        iso = list(range(n, 2*n)) + list(range(n))
        assert code2.get_dual().equiv(code2.apply_perm(iso))

        code2 = code2.to_css()
        #code2.get_params()
        print(code, "CSS" if code.is_css() else " ", end=" ", flush=True)
        print("--->", code2, end=" ", flush=True)
        d = code2.x_distance()
        print("d =", d)
        #iso = code2.find_zx_duality()
        #print("iso", iso)
        
        code = zxcat(code2, iso)
        if code is None:
            print("FAIL")
            continue
        assert code.is_selfdual()
        print("\t", code, "self-dual")

        if code.n <= 10:
            code = code.to_css()
            d = code.x_distance()
            print("\td =", d)


#test = test_zx

def get_jaunty():
    # https://arxiv.org/abs/2010.06628
    code = QCode.fromstr("""
    XXX..X........
    ..XX..XX......
    .......XX..XX.
    X........X..XX
    .YYYY.........
    ..Y..YY...Y...
    ....Y...YY..Y.
    ..........YYYY
    ZZ..Z....Z....
    ...ZZ..ZZ.....
    ......ZZ..ZZ..
    """, check=True)
    return code


def get_prism():
    code = QCode.fromstr("""
    XXXX....
    ZZ..ZZ..
    ....XXXX
    ..ZZ..ZZ
    Y.Z.X...
    .Y.Z.X..
    ..X.Z.Y.
    """, check=True)
    return code


def test():

    #code = get_jaunty()
    #code = get_prism()

    code = QCode.fromstr("ZYX. .XZY YZ.X")

    print(code)
    print(code.get_params())
    print(code.longstr())

    code2 = unwrap(code)
    print(code2)
    print(code2.longstr())
    print(code2.is_selfdual())
    print(code2.to_css().distance())

    toric = construct.toric(2, 2, 1)
    toric = toric.to_qcode()
    assert code2.is_isomorphic(toric)



if __name__ == "__main__":

    start_time = time()


    profile = argv.profile
    name = argv.next() or "test"
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()


    t = time() - start_time
    print("OK! finished in %.3f seconds\n"%t)





