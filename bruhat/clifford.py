#!/usr/bin/env python3

from functools import reduce
from operator import matmul, add, mul

from bruhat import element
from bruhat.vec import Space, Hom, Map
from bruhat.action import mulclose, mulclose_names, Perm, Group
from bruhat.argv import argv
from bruhat.util import factorial, partitions
from bruhat.rep import get_perms, Young

def atpow(A, n):
    return reduce(matmul, [A]*n)

def main():

    Z = element.Z
    ring = element.CyclotomicRing(8)
    x = ring.x
    r2 = x + x**7
    assert r2**2 == 2

    qubit = Space(2, ring)

    q2 = qubit @ qubit

    hom = Hom(qubit, qubit)

    I = Map.from_array([[1, 0], [0, 1]], hom)
    X = Map.from_array([[0, 1], [1, 0]], hom)
    Z = Map.from_array([[1, 0], [0, -1]], hom)
    Y = Map.from_array([[0, x**6], [x**2, 0]], hom)
    assert Y*Y == I
    assert Y == (x**2)*X*Z

    II = I@I
    XI = X@I
    IX = I@X
    XX = X@X
    ZI = Z@I
    IZ = I@Z

    assert XI*IX == XX
    T = Map.from_array([[1, 0], [0, x]], hom)
    Td = Map.from_array([[1, 0], [0, x**7]], hom)
    S = T*T
    Sd = Td*Td
    assert Z == S*S
    assert T*Td == I
    assert S*Sd == I
    assert S*Z == Z*S
    assert T*Z == Z*T

    assert (Z*X*Z) == -X
    assert ((Z@Z)*(X@X)*(Z@Z)) == X@X
    assert ((S@S)*(X@X)*(Sd@Sd)) == atpow((S*X*Sd), 2)

    SXSd = S*X*Sd
    TXTd = T*X*Td
    lhs = r2 * TXTd
    i = x**2
    assert(lhs == X + i*X*Z)

    if 0:
        names = mulclose_names([X, Z, x*I], "X Z xI".split())
        #print("SXSd:", "*".join(names[SXSd]))
    
        names = mulclose_names([XI, IX, ZI, IZ], "XI IX ZI IZ".split())
        g = atpow(SXSd, 2)
        #print("*".join(names[g]))
        
        gen = [
            X@I@I@I, I@X@I@I, I@I@X@I, I@I@I@X,
            Z@I@I@I, I@Z@I@I, I@I@Z@I, I@I@I@Z, 
            x*(I@I@I@I)]
        names = "XIII IXII IIXI IIIX ZIII IZII IIZI IIIZ wIIII".split()
        names = mulclose_names(gen, names)
        g = atpow(SXSd, 4)
        #print("*".join(names[g]))
        g = atpow(TXTd, 4)
        assert g not in names

    I8, X8, Z8 = atpow(I, 8), atpow(X, 8), atpow(Z, 8)
    T8 = atpow(T, 8)
    #print(T8.str(hide_zero=True, sep=""))
    Td8 = atpow(Td, 8)
    G = mulclose([X8, Z8])
    P = I8 + X8 + Z8 + X8*Z8
    assert P*P == 4*P
    lhs = P
    rhs = T8*P*Td8
    assert(rhs*rhs == 4*rhs)

    #T1 = reduce(matmul, [T, Td]*4)
    #T2 = reduce(matmul, [Td, T]*4)
    #print(T1*P*T2)
    #return

    Q = P@P

    RM14 = """
    1111111111111111
    .1.1.1.1.1.1.1.1
    ..11..11..11..11
    ....1111....1111
    ........11111111
    """.strip().split()
    RM14.pop(0)

    def interp(row, op):
        ops = [op if s=="1" else I for s in row]
        A = reduce(matmul, ops)
        return A

    RM14 = [row[1:] for row in RM14] # puncture
    n = len(RM14[0])
    Sx = [interp(row, X) for row in RM14]
    Sz = [interp(row, Z) for row in RM14]

    if 0:
        print("mulclose...")
        G = mulclose(Sx+Sz, maxsize=2**len(Sx+Sz)) # too slow :-(
        print(len(G))
        P = reduce(add, G)
        assert P*P == len(G)*P
    
        Tn = atpow(T, n)
        Tdn = atpow(Td, n)
        print(P == Tn*P*Tdn)

    #print(lhs.str(hide_zero=True, sep=""))
    #print()
    #print()
    #print()
    #print(rhs.str(hide_zero=True, sep=""))

    gen = [
        X@I@I@I@I@I@I@I,
        I@X@I@I@I@I@I@I,
        I@I@X@I@I@I@I@I,
        I@I@I@X@I@I@I@I,
        I@I@I@I@X@I@I@I,
        I@I@I@I@I@X@I@I,
        I@I@I@I@I@I@X@I,
        I@I@I@I@I@I@I@X,
        Z@I@I@I@I@I@I@I,
        I@Z@I@I@I@I@I@I,
        I@I@Z@I@I@I@I@I,
        I@I@I@Z@I@I@I@I,
        I@I@I@I@Z@I@I@I,
        I@I@I@I@I@Z@I@I,
        I@I@I@I@I@I@Z@I,
        I@I@I@I@I@I@I@Z,
        x*(I@I@I@I@I@I@I@I)]
    names = """
        XIIIIIII IXIIIIII IIXIIIII IIIXIIII IIIIXIII IIIIIXII IIIIIIXI IIIIIIIX 
        ZIIIIIII IZIIIIII IIZIIIII IIIZIIII IIIIZIII IIIIIZII IIIIIIZI IIIIIIIZ 
        wIIIIIIII
    """.strip().split()
    gen.pop(-1)
    names.pop(-1)

    if 0:
        g = atpow(TXTd, 8)
        for i in range(8):
            print(g == (x**i)*atpow(X*Z, 8))
        return

    if 0:
        names = mulclose_names(gen, names)
        g = atpow(TXTd, 8)
        for i in range(8):
            print(names.get((x**i)*g, "?"))
    
    #for i in [1, 2, 4]:
    #    print(atpow(TXTd, i))

    if 0:
        #print(Td*X*T)
        X4 = X@X@X@X
        T4 = T@T@T@T
        Td4 = Td@Td@Td@Td
        print(Td4 * X4 * T4)

    H = X+Z
    E = Map.from_array([[0, 1], [0, 0]], hom) # raising 
    F = Map.from_array([[0, 0], [1, 0]], hom) # lowering
    
    CNOT = Map.from_array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1],
         [0, 0, 1, 0]], hom@hom)

    SWAP = Map.from_array(
        [[1, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1]], hom@hom)

    assert SWAP * SWAP == II
    assert CNOT * CNOT == II
    
    G = mulclose([XI, IX, CNOT]) # order 8
    G = mulclose([XI, IX, CNOT, SWAP]) # order 24.. must be S_4
    G = mulclose([CNOT, SWAP])
#    print(len(G))
#    for g in G:
#        print(g)

    A = SWAP @ I
    B = I @ SWAP
    S_3 = list(mulclose([A, B]))
    assert len(S_3) == 6
#    for g in G:
#        print(g)
#    print(g.hom)

    hom = A.hom
    space = hom.src
    N = 2**3
    basis = space.get_basis()
    orbits = set()
    for v in basis:
        v1 = space.zero_vector()
        for g in S_3:
            u = g*v
            v1 = v1 + u
        orbits.add(v1)
    orbits = list(orbits)
#    for v in orbits:
#        print(v)

    HHH = H@H@H
    v = basis[7]
    #print(v)
    #print(HHH * v)

    print("OK")


if __name__ == "__main__":

    name = argv.next() or "main"
    fn = eval(name)
    fn()


