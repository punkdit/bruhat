#!/usr/bin/env python
"""
doctrines for various code families.
"""

if 1:
    import os
    os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
    os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=1
    os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=1
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=1
    os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=1


from functools import cache, reduce
from operator import matmul, mul
import numpy

from bruhat.dev.geometry import all_codes
from bruhat.sp_pascal import i_grassmannian
from bruhat.qcode import QCode, strop
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
        H = zeros2(m, 2*n)
        H.shape = m, n, 2
        yield H
        return 

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
            #count += 1
            yield H1
    if verbose:
        print()
    #return count


def show(H1):
    if not len(H1):
        return
    H = H1.view()
    m, n, _ = H.shape
    nn = 2*n
    H.shape = m, nn
    print()
    print(strop(H))
    #if "X." in strop(H) and '.Z' in strop(H):
    #    print(has_transversal_CZ(H1))


def show_zxcat(H):
    from bruhat.unwrap import zxcat
    from qumba.qcode import QCode
    code = QCode(H)
    #print(code.longstr())
    perm = []
    for i in range(code.n):
        perm.append(i-1 if i%2 else i+1)
    dode = zxcat(code, perm)
    eode = dode.apply_H()
    if dode.equiv(eode):
        print("*", end="", flush=True)
    else:
        print(".", end="", flush=True)
    


def show_832(H):
    from bruhat.unwrap import zxcat
    from qumba.qcode import QCode
    from qumba import construct
    right = QCode(H)
    print(strop(right.H))
    print("-->")
    inner = construct.get_832()
    left = (right.n//3) * inner
    right = QCode.trivial(left.n - right.n) + right
    code = left << right
    print(strop(code.H))
    print()


    

def has_nothing(H1):
    return True


S = array2([[1,1],[0,1]])
H = array2([[0,1],[1,0]])
SH = dot2(S, H)
SHS = dot2(SH, S)
def has_transversal_S_upto_sign(H1):
    H2 = dot2(H1, S)
    return equal(H1, H2)

def has_transversal_upto_sign(H1, op):
    H2 = dot2(H1, op)
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


@cache 
def get_transversal_H(n):
    from qupy.dense import Gate
    L = reduce(matmul, [Gate.H]*n)
    return L

@cache 
def get_transversal_S(n):
    from qupy.dense import Gate
    L = reduce(matmul, [Gate.S]*n)
    return L

@cache 
def get_transversal_T(n):
    from qupy.dense import Gate
    L = reduce(matmul, [Gate.T]*n)
    return L

def has_transversal_T(H1):
    if not has_transversal_S_upto_sign(H1):
        return False
    m, n, _ = H1.shape
    P = get_projector(H1)
    L = get_transversal_T(n)
    return P*L == L*P

@cache 
def get_transversal_THT(n):
    from qupy.dense import Gate
    U = Gate.T * Gate.H * Gate.T
    L = reduce(matmul, [U]*n)
    return L

def has_transversal_THT(H1):
    m, n, _ = H1.shape
    P = get_projector(H1)
    L = get_transversal_THT(n)
    return P*L == L*P

@cache 
def get_transversal_SH(n):
    from qupy.dense import Gate
    U = Gate.S * Gate.H
    L = reduce(matmul, [U]*n)
    return L

@cache
def get_transversal_CZ(n):
    assert n%2 == 0
    from qupy.dense import Qu
    CZ = Qu((2,)*4, 'uudd')
    CZ[0, 0, 0, 0] = 1.
    CZ[0, 1, 0, 1] = 1.
    CZ[1, 0, 1, 0] = 1.
    CZ[1, 1, 1, 1] = -1.
    L = reduce(matmul, [CZ]*(n//2))
    return L


def dump(*args):
    for P in args:
        n = len(P.shape)//2
        v = P.v
        assert numpy.abs(v.imag).sum() < 1e-6
        v = v.real
        v = v.astype(int)
        v.shape = (2**n, 2**n)
        print(v)


def has_transversal_CZ(H1):
    m, n, _ = H1.shape
    if n%2:
        return False
    P = get_projector(H1)
    L = get_transversal_CZ(n)
    #dump(P, L, P*L, L*P)
    return P*L == L*P

@cache
def get_transversal_CCZ(n):
    assert n%3 == 0
    from qupy.dense import Qu
    #CCZ = Qu((2,)*6, 'ud'*3)
    CCZ = Qu((2,)*6, 'u'*3 + 'd'*3)
    for a in [0,1]:
     for b in [0,1]:
      for c in [0,1]:
        CCZ[a, b, c, a, b, c] = -1 if a==b==c==1 else 1
    L = reduce(matmul, [CCZ]*(n//3))
    return L

def has_transversal_CCZ(H1):
    m, n, _ = H1.shape
    if n%3:
        return False
    P = get_projector(H1)
    L = get_transversal_CCZ(n)
    return P*L == L*P

def has_transversal_SH(H1):
    m, n, _ = H1.shape
    if not has_transversal_upto_sign(H1, SH):
        return False
    P = get_projector(H1)
    L = get_transversal_SH(n)
    return P*L == L*P

@cache 
def get_transversal_SHS(n):
    from qupy.dense import Gate
    U = Gate.S * Gate.H * Gate.S
    L = reduce(matmul, [U]*n)
    return L

def has_transversal_SHS(H1):
    m, n, _ = H1.shape
    if not has_transversal_upto_sign(H1, SHS):
        return False
    P = get_projector(H1)
    L = get_transversal_SHS(n)
    return P*L == L*P

@cache 
def get_transversal_SHSdag(n):
    from qupy.dense import Gate
    U = Gate.S * Gate.H * ~Gate.S
    L = reduce(matmul, [U]*n)
    return L

def has_transversal_SHSdag(H1):
    m, n, _ = H1.shape
    if not has_transversal_upto_sign(H1, SHS):
        #assert P*L!=L*P
        return False
    P = get_projector(H1)
    L = get_transversal_SHSdag(n)
    return P*L == L*P


def has_transversal_clifford(H1):
    m, n, _ = H1.shape
    if not has_transversal_upto_sign(H1, S):
        return False
    if not has_transversal_upto_sign(H1, H):
        return False
    P = get_projector(H1)
    L = get_transversal_S(n)
    if P*L != L*P:
        return False
    L = get_transversal_H(n)
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


def clifford():
    from qupy.dense import Gate
    from qupy.util import mulclose

    G = mulclose([Gate.S, Gate.H])
    print(len(G))


def find_clifford():
    from qumba.qcode import QCode, dot2, eq2
    from qumba.matrix import SymplecticSpace, mulclose, Matrix
    n = 4
    space = SymplecticSpace(n)
    S = space.get_S()
    H = space.get_H()
    c1 = mulclose([S, H])
    assert len(c1) == 6 # Sp(2,2)

    gen = [space.get_S(i) for i in range(4)]
    gen += [space.get_H(i) for i in range(4)]
    G = mulclose(gen)
    print(len(G))

    CZ = space.get_CZ(2,3)
    HI = space.get_H(2)
    IH = space.get_H(3)
    SI = space.get_S(2)
    IS = space.get_S(3)
    II = space.identity()
    swap = space.get_perm([0,1,3,2])
    c2 = mulclose([HI, IH, SI, IS, CZ])
    assert len(c2) == 720 # Sp(4,2)
    c2.remove(II)

    if 0:
        code = QCode.fromstr("XXXX ZZZZ", None, "XXII ZIZI XIXI ZZII")
        print(code.longstr())
        print()
        dode = code.apply_H()
        print(dode.longstr())
        print()
        assert code.equiv(dode)
        L = code.get_logical(dode)
        #print(L)
        #print(SymplecticSpace(2).get_perm([1,0]))
        #return
    
        E = Matrix(code.get_encoder())
        print(E)
    
        print()
        print(E*swap*HI*IH)
        print()
        print(H*E)
    
        return

    tgt = SymplecticSpace(2).get_CNOT(0,1).A # yes
    #tgt = SymplecticSpace(2).get_CZ(0,1).A # no...

    found = []

    count = 0
    for H in search(4, 2, has_nothing):
        code = QCode(H)
        if code.get_distance() < 2:
            continue
        #print(code.longstr())
        E = Matrix(code.get_encoder())
        for g in [S]:
        #for g in c1:
        #for g in G:
            if g==II:
                continue
            F = g*E
            dode = QCode.from_encoder(F.A, m=code.m)
            #print(dode.longstr())
            if dode.equiv(code):
                #print(g)
                #print(strop(code.H))
                L = code.get_logical(dode)
                if eq2(L, tgt):
                    #print()
                    #print(code.longstr())
                    print("*", end="", flush=True)
                    found.append(code)
                else:
                    print(".", end="", flush=True)
        count += 1
    print()
    print(count)
    for code in found:
        print(strop(code.H))
        print()


def main():

    accept = argv.get("accept", "has_transversal_S")
    print("accept", accept)
    accept = eval(accept)

    n = argv.get("n")
    m = argv.get("m", 1)
    if n is not None:
        count = 0
        for H in search(n, m, accept):
            post(H)
            count += 1
        return

    post = argv.get("post", "lambda H:None")
    post = eval(post)

    mod = argv.get("mod", 0)

    n0 = argv.get("n0", 1)
    n1 = argv.get("n1", 10)
    #for n in range(6, 7):
    #  for m in range(4, n+1):
    for n in range(n0, n1):
      if mod!=0 and n%mod:
        continue
      print("n=%2d:"%n, end=" ", flush=True)
      for m in range(n+1):
        count = 0
        for H in search(n, m, accept):
            post(H)
            count += 1
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



