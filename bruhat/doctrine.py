#!/usr/bin/env python
"""
doctrines for various code families.
"""

from functools import cache, reduce
from operator import matmul, mul
import numpy

from bruhat.dev.geometry import all_codes
from bruhat.sp_pascal import i_grassmannian
from bruhat.qcode import QCode
from bruhat.solve import zeros2, enum2, dot2, shortstr, array2
from bruhat.util import choose
from bruhat.argv import argv

if argv.numba:
    from bruhat.orthogonal import normal_form # numba version
else:
    from bruhat.algebraic import normal_form_p as normal_form


from time import time
start_time = time()

CHECK = False

def equal(lhs, rhs):
    assert lhs.shape == rhs.shape
    m, n, k = lhs.shape
    lhs = lhs.view()
    rhs = rhs.view()
    lhs.shape = (m, 2*n)
    rhs.shape = (m, 2*n)
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


def search(n, m, accept, verbose=False):

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
        if accept(H1):
            if verbose:
                print(".", end='', flush=True)
            count += 1
    if verbose:
        print()
    return count


S = array2([[1,1],[0,1]])
def has_transversal_S_upto_sign(H1):
    H2 = dot2(H1, S)
    return equal(H1, H2)


# slightly faster with cache... not the bottleneck 
#@cache
def get_op(desc):
    from qupy.dense import Gate
    lookup = {(0, 0) : Gate.I, (1, 0) : Gate.X, (0, 1) : Gate.Z, (1, 1) : Gate.Y}
    op = reduce(matmul, [lookup[d] for d in desc])
    return op

def get_projector(H1):
    m, n, _ = H1.shape
    from qupy.dense import Gate
    I, X, Z, Y, S, T = Gate.I, Gate.X, Gate.Z, Gate.Y, Gate.S, Gate.T
    stabs = []
    for row in H1:
        desc = tuple(tuple(bit) for bit in row)
        stabs.append(get_op(desc))
    for g in stabs:
      for h in stabs:
        assert g*h == h*g
    In = reduce(matmul, [I]*n)
    P = reduce(mul, [In + stab for stab in stabs])
    return P

def has_transversal_gate(H1, *gates):
    #result = has_transversal_S_upto_sign(H1)
    #if not result:
    #    return False
    m, n, _ = H1.shape
    P = get_projector(H1)
    L = reduce(matmul, [gates[i%len(gates)] for i in range(n)])
    return P*L == L*P


def has_transversal_S(H1):
    if not has_transversal_S_upto_sign(H1):
        return False
    from qupy.dense import Gate
    m, n, _ = H1.shape
    P = get_projector(H1)
    L = reduce(matmul, [Gate.S]*n)
    return P*L == L*P


def has_transversal_SSdag(H1):
    if not has_transversal_S_upto_sign(H1):
        return False
    from qupy.dense import Gate
    m, n, _ = H1.shape
    P = get_projector(H1)
    L = reduce(matmul, [[Gate.S, ~Gate.S][i%2] for i in range(n)])
    return P*L == L*P


def has_transversal_TTdag(H1):
    if not has_transversal_S_upto_sign(H1):
        return False
    from qupy.dense import Gate
    m, n, _ = H1.shape
    P = get_projector(H1)
    L = reduce(matmul, [[Gate.T, ~Gate.T][i%2] for i in range(n)])
    return P*L == L*P


def is_cyclic(H1):
    m, n, _ = H1.shape
    cols = [(i+1)%n for i in range(n)]
    H2 = H1[:, cols]
    return equal(H1, H2)


def is_symmetric(H1):
    m, n, _ = H1.shape
    for j in range(n-1):
        cols = list(range(n))
        cols[j:j+2] = [j+1, j]
        #print(cols)
        H2 = H1[:, cols, :]
        #print(H1.shape, H2.shape)
        if not equal(H1.view(), H2):
            return False
    return True


def main():

    accept = argv.get("accept", "has_transversal_S")
    print("accept", accept)
    accept = eval(accept)

    n = argv.get("n")
    m = argv.get("m", 1)
    if n is not None:
        count = search(n, m, accept, argv.verbose)
        print(count)
        return

    n0 = argv.get("n0", 1)
    n1 = argv.get("n1", 10)
    #for n in range(6, 7):
    #  for m in range(4, n+1):
    for n in range(n0, n1):
      for m in range(n+1):
        count = search(n, m, accept)
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



