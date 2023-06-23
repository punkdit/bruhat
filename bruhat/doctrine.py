#!/usr/bin/env python
"""
doctrines for various code families.
"""

from functools import cache
import numpy

from bruhat.dev.geometry import all_codes
from bruhat.sp_pascal import i_grassmannian
from bruhat.qcode import QCode
from bruhat.solve import zeros2, enum2, dot2, shortstr
#from bruhat.algebraic import normal_form_p as normal_form
from bruhat.orthogonal import normal_form # numba version
from bruhat.util import choose
from bruhat.argv import argv

from time import time
start_time = time()

CHECK = False

def equal(a_code, b_code):
    lhs = a_code.flatH
    rhs = b_code.flatH
    assert lhs.shape == rhs.shape
    m, n = lhs.shape
    #print()
    #print("equal")
    #print(shortstr(lhs))
    #print('-'*n)
    #print(shortstr(rhs))
    rr_lhs = normal_form(lhs)
    rr_rhs = normal_form(rhs)
    result = numpy.alltrue(rr_lhs==rr_rhs)
    if CHECK:
        lhs = [dot2(v, lhs).tobytes() for v in enum2(m)]
        rhs = [dot2(v, rhs).tobytes() for v in enum2(m)]
        lhs.sort()
        rhs.sort()
        if (lhs==rhs) != result:
            print("fail", lhs==rhs)
            print(shortstr(rr_lhs))
            print('-'*n)
            print(shortstr(rr_rhs))
            assert 0
    return result


def test_equal():
    lhs = QCode.fromstr("XX YY")
    print(lhs)

    rhs = lhs.apply_S(0).apply_S(1)
    print(rhs)
    assert equal(lhs, rhs)
    return


def search(n, m):

    if m==0:
        return 1

    nn = 2*n
    refl = zeros2(nn, nn)
    cols = []
    for i in range(n):
        cols.append(i)
        cols.append(2*n-i-1)
    #print(cols)

    count = 0
    #for H in all_codes(m, nn):
    for piv,H in i_grassmannian(n, m):
        #print(H)
        H1 = H[:, cols]
        #print(H1)
        H1.shape = m, n, 2
        code = QCode(H1)
        c = code
        for idx in range(n):
            c = c.apply_S(idx)
        #if str(code) != str(c):
        #    continue
        if 0:
            print(code)
            print("-"*n)
            print(c)
            print()
        if not equal(code, c):
            continue
        #if str(code) != str(c):
        #    print(code)
        #    print("-"*n)
        #    print(c)
        #    print()
        count += 1
    return count


def main():
    n = argv.get("n")
    m = argv.get("m", 1)
    if n is not None:
        count = search(n, m)
        print(count)
        return

    #for n in range(1, 6):
    for n in range(6, 7):
      for m in range(n+1):
        count = search(n, m)
        print("%3d"%count, end=" ", flush=True)
      print()


if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()



