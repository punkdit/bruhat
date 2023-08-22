#!/usr/bin/env python

from time import time
start_time = time()
from random import shuffle, choice

import numpy

from qumba.solve import (parse, shortstr, linear_independent, eq2, dot2, identity2,
    rank, rand2, pseudo_inverse, kernel, direct_sum, zeros2, solve2, normal_form)
from qumba.qcode import QCode, SymplecticSpace
from qumba import construct
from qumba.autos import get_isos
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


def zxcat(code, duality):
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
        print(code, "CSS" if code.is_css() else " ")
        print("--->", code2, end=" ", flush=True)
        d = code2.x_distance()
        print("distance =", d)
        #iso = code2.find_zx_duality()
        #print("iso", iso)
        
        code = zxcat(code2, iso)
        if code is None:
            print("FAIL")
            continue
        assert code.is_selfdual()
        print(code, "self-dual")


test = test_zx


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





