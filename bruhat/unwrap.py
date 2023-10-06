#!/usr/bin/env python

from time import time
start_time = time()
from random import shuffle, choice, randint
from functools import reduce
from operator import add

import numpy

from bruhat.action import Perm

from qumba.solve import (parse, shortstr, linear_independent, eq2, dot2, identity2,
    rank, rand2, pseudo_inverse, kernel, direct_sum, zeros2, solve2, normal_form, array2)
from qumba.qcode import QCode, SymplecticSpace, Matrix, strop, fromstr
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


def get_pairs(perm):
    pairs = []
    for (i, j) in enumerate(perm):
        if i==j:
            assert 0
        assert i!=j, "perm has fixed point"
        assert perm[j] == i
        if i < j:
            pairs.append((i, j))
    assert len(pairs)*2 == len(perm), "perm is not an involution?"
    return pairs


def wrap(code, duality):
    pairs = get_pairs(duality)

    rows = []
    for a in code.Hx:
        row = []
        for (i, j) in pairs:
            row.append(a[i])
            row.append(a[j])
        rows.append(row)
    H = array2(rows)
    code = QCode(H)
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


def unwrap_encoder(code):
    E = code.get_encoder()
    Ei = code.space.invert(E)
    space = SymplecticSpace(2*code.n)

    n, m, k = code.n, code.m, code.k
    E2 = zeros2(4*n, 4*n)
    E2[::2, ::2] = E
    E2[1::2, 1::2] = Ei.t
    E2 = Matrix(E2)
    assert space.is_symplectic(E2)
    F = space.F

    perm = list(range(4*n))
    for i in range(m):
        a, b = perm[4*i+2:4*i+4]
        perm[4*i+2:4*i+4] = b, a
    E2 = E2[:, perm]
    assert space.is_symplectic(E2)

    #HT = E2.t[:4*m, :]
    #print(strop(HT))
    #print()

    code2 = QCode.from_encoder(E2, 2*m)
    #print(code2.longstr(), code2)
    return code2


def test_513():

#    H = """
#    XZZX.
#    .XZZX
#    X.XZZ
#    ZX.XZ
#    """
#    T = """
#    XZZZX
#    XXZZZ
#    ZXXZZ
#    ZZXXZ
#    """
#
#    H = Matrix(fromstr(H))
#    T = Matrix(fromstr(T))
#    #print(QCode(H, T))
#    space = SymplecticSpace(5)
#    F = space.F
#    #print(H * F * T.t)
#
#    code = construct.get_513()
#    print(code.longstr())
#
#    n = code.n
#    perm = [(i+1)%n for i in range(n)]
#    R = space.get_perm(perm)
#
#    T, L = code.T, code.L
#    TL = T.concatenate(L)
#    nn = code.nn
#    bits = [[] for i in range(code.m)]
#    for t in numpy.ndindex((2,)*nn):
#        t = array2(t)
#        t.shape = 1, nn
#        t = Matrix(t)
#        u = (t*F*H.t)
#        if u.sum() == 1:
#            idx = u.where()[0][1]
#            bits[idx].append(t)
#
#    for bit in bits:
#      for t in bit:
#        t1 = t*R
#        u = (t*F*H.t)
#        u1 = (t1*F*H.t)
#        if u1.sum() == 1:
#            print(strop(t), strop(t1), u, u1)
#
#    return

    code = construct.get_513()
    code2 = unwrap_encoder(code)

    iso = code2.get_isomorphism( unwrap(code) )
    assert iso is not None
    print(iso)


def test_tetrahedron():

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


def fixed(f):
    return [i for i in range(len(f)) if f[i]==i]

def is_identity(f):
    for i in range(len(f)):
        if f[i] != i:
            return False
    return True

def mul(f, g):
    return [f[g[i]] for i in range(len(g))]


def test_unwrap():
    from bruhat.doctrine import search
    from bruhat.sp_pascal import i_grassmannian

    for (h,t) in [
        #("Z", "X"),
        #("XZ ZX", "ZI IZ"),
        #("XX IX", "ZI ZZ"),
        #("XI IX", "ZX XZ"), 
        ("YI IY", "ZI IZ"),
    ]:

        code = QCode.fromstr(h, t)
        #assert code.get_encoder() == code.space.get_CNOT()
        print(code.longstr())
        print(code.get_encoder())
    
        code = unwrap_encoder(code)
        space = code.space
        n = space.n
        print(code.longstr())
        E = code.get_encoder()
        print(E)
        #rhs = space.get_CNOT(0, 3)*space.get_CNOT(1, 2)*space.get_H(2)*space.get_H(3)
        #print(E==rhs)

        print()

    CX = space.get_CNOT
    H = space.get_H

    E0 = E
    word = []
    for i in range(n):
      for j in range(n):
        if j==i:
            continue
        A = space.get_CNOT(i, j)
        if str(E).count('1') > str(A*E).count('1'):
            print("get_CNOT", i, j)
            E = A*E

    print(E)

    print(CX(3,0)*CX(1,2)*H(1)*H(3) == E0)



def test():
    src = construct.toric(2, 2)
    print(src)
    #print(src.longstr())
    Ax, Az = src.Ax, src.Az
    print("Ax =")
    print(shortstr(Ax))
    print(''.join(str(i) for i in range(src.n)))
    print("Az =")
    print(shortstr(Az))
    print(''.join(str(i) for i in range(src.n)))

    duality = find_zx_duality(Ax, Az)
    autos = find_autos(Ax, Az)

    #for auto in autos:
    #    if len(fixed(auto)) == 0 and is_identity(mul(auto, auto)):
    #        print("auto:", auto, get_pairs(auto))

    codes = []
    for auto in autos:
        f = mul(duality, auto)
        if len(fixed(f)) == 0 and is_identity(mul(f, f)):
            #code = zxcat(src, f)
            #print(code.get_params())
            pairs = get_pairs(f)
            print(pairs)
            dode = wrap(src, f)
            print("H =")
            #print(strop(dode.H), dode.get_params())
            print(dode.longstr())
            codes.append(dode)
            #print(shortstr(code.H))
            #codes.append(zxcat(src, f))
            print()

    print(len(codes))

    from bruhat.equ import quotient
    lookup = quotient(codes, lambda a,b:a.is_isomorphic(b))
    found = []
    for k,v in lookup.items():
        if v in found:
            continue
        found.append(v)
        print('\t', len(v))


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





