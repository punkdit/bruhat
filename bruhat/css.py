#!/usr/bin/env python

"""
build CSS codes and self-dual CSS codes
"""


from time import time
start_time = time()

from functools import lru_cache
cache = lru_cache(maxsize=None)

import numpy
from numpy import alltrue

from bruhat.solve import find_kernel, dot2, shortstr
from bruhat.algebraic import (
    qchoose_2, all_codes, normal_form_p, Algebraic, Matrix)
from bruhat.sp_pascal import grassmannian
from bruhat.argv import argv


@cache
def find_css(n, mx, mz):
    if mx < mz:
        return find_css(n, mz, mx) # recurse
    count = 0
    #found = set()
    for Hx in qchoose_2(n, mx):
        Jz = find_kernel(Hx)
        assert mx + len(Jz) == n
        assert dot2(Hx, Jz.transpose()).sum() == 0
        for Kz in qchoose_2(len(Jz), mz):
            Hz = dot2(Kz, Jz)
            Hz = normal_form_p(Hz)
            assert dot2(Hx, Hz.transpose()).sum() == 0
            key = str((Hx, Hz))
            #assert key not in found
            #found.add(key)
            count += 1
    return count


def main_css():

    n = argv.get("n", 5)
    m = argv.get("m", n//2)
    mx = argv.get("mx", m)
    mz = argv.get("mz", mx)
    assert mx + mz <= n

    for n in range(10):
        print("n:", n)
        for mx in range(n+1):
          print("   ", end="")
          for mz in range(n+1):
            if mx+mz>n:
                continue
            count = find_css(n, mx, mz)
            print("(%s,%s):%8d"%(mx, mz, count,), end=" ", flush=True)
          print()



def dump(Hx):
    from qupy.ldpc import solve
    solve.int_scalar = numpy.int64
    from qupy.ldpc.css import CSSCode
    code = CSSCode(Hx=Hx, Hz=Hx)
    dx = code.x_distance()
    assert dx>1
    if dx==2:
        print(".", flush=True, end="")
    else:
        print(dx)
        print(Hx)


def main_selfdual():

    n = argv.get("n", 5)
    m = argv.get("m", n//2)
    show = argv.show
    assert 2*m <= n

    count = 0
    found = 0
    for Hx in qchoose_2(n, m):
        if dot2(Hx, Hx.transpose()).sum() == 0:
            if Hx.sum(0).min() > 0:
                if show:
                    #print(Hx)
                    dump(Hx)
                found += 1
            count += 1

    print()
    print("found: %d, total: %d"%(found, count))


def search(n, m, G, M):
    F = G.invariant_form
    found = set()
    for g in G:
# no solution...
#        A = g.A
#        if A[n:2*n, :2*n+1].max():
#            continue
#        if A[:n, 2*n+1:].max():
#            continue
#        if A[3*n+1:4*n+1, :2*n+1].max():
#            continue
#        if A[4*n+1:, 2*n+1:].max():
#            continue
        M1 = M*g
        assert M.shape == M1.shape
        A = M1.A
        if A[:m, 2*n+1:].max():
            continue
        if A[m:, :2*n+1].max():
            continue
        #print('.', end='', flush=True)
        H1 = A[:m, :]
        H2 = A[m:, :]
        H11 = dot2(H1, F)
        if numpy.alltrue(H11==H2):
            #print('M', end='', flush=True)
            #print(shortstr(g))
            #print()
            #return g
            found.add(g)
    print("found:", len(found))
    return found


def find_selfdual():

    n = argv.get("n", 2)
    m = argv.get("m", 1)

    G = Algebraic.Sp(4*n+2)
    print(len(G))
    F = G.invariant_form


    sol = None
    for _, H in grassmannian(n, m):
        A = H[:, :n]
        B = H[:, n:]
        print(A, B)

        M = numpy.zeros((2*m, 2*(2*n+1)), dtype=int)
        M[:m, :n] = A
        M[m:2*m, n:2*n] = A

        M[m:2*m, 2*n+2:2*n+2+n] = B
        M[:m, 2*n+2+n:] = B
        M[:, 2*n] = 1
        M[:, 2*n+1] = 1

        M = Matrix(M)
        print(M)

        # isotropic
        lhs = M * F * M.transpose()
        assert lhs.is_zero()

        found = search(n, m, G, M)
        sol = found if sol is None else sol.intersection(found)
        print("sol:", len(sol))




if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()


    print("finished in %.3f seconds"%(time() - start_time))


