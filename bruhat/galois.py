#!/usr/bin/env python3

from random import choice

import numpy
from numpy import dot, array, empty
from numpy.linalg import eigvals, eig

from bruhat.element import PolynomialRing, Z, cyclotomic
from bruhat.argv import argv
from bruhat.action import Perm, Group


def test():
    a = (0.5)*(1+5**0.5)
    b = (0.5)*(1-5**0.5)
    
    A = array([[0., 1], [1, 1]])
    B = array([[1., -1], [-1, 0]])
    
    print(dot(A, B))
    print(dot(B, A))
    
    vals, vecs = eig(A)
    print("vals:", vals)
    print("vecs:")
    print(vecs)
    
    #print(eigvals(B))
    
    u = vecs[:, 0]
    Bu = dot(B, u) / vals[1]
    print(u, "==", Bu)


def astr(v):
    vs = empty(v.shape, dtype=object)
    for idx in numpy.ndindex(v.shape):
        vs[idx] = str(v[idx])
    s = str(vs)
    s = s.replace("'", "")
    return s


def test_quadratic():

    # here we construct a quadratic extension

    R = PolynomialRing(Z)
    X = R.x
    #p = X**2 - X - 1 # fibonacci
    p = X**2 + 2*X + 7

    S = R/p

    for i in range(-100, 100):
        assert p(i) != 0

    X1 = S.promote(X)
    # X1 + X2 == -p[1] 
    X2 = -p[1] - X1

    v = array([1, X])
    print(astr(v))

    v1 = X1*v
    v2 = X2*v

    deg = 2
    #print(astr(v1))
    A1 = empty((deg, deg), dtype=object)
    
    for i in range(deg):
      for j in range(deg):
        A1[i, j] = v1[i][j]

    print(astr(A1))

    A2 = empty((deg, deg), dtype=object)
    
    for i in range(deg):
      for j in range(deg):
        A2[i, j] = v2[i][j]

    print(astr(A2))



def test_cyclotomic():

    # here we construct a cyclotomic extension

    n = argv.get("n", 4)

    GL_1 = list(range(1, n))
    for i in range(1, n):
      for j in range(1, n):
        if (i*j)%n:
            continue
        if i in GL_1:
            GL_1.remove(i)
        if j in GL_1:
            GL_1.remove(j)

    print("non-zero-divisors:", GL_1)
    for i in GL_1:
      for j in GL_1:
        print("%2d"%((i*j)%n), end=" ")
      print()
    print()

    R = PolynomialRing(Z)
    X = R.x
    p = cyclotomic(R, n)
    print("p:", p)
    deg = p.deg

    S = R/p
    X = S.promote(X)

    v = [X**i for i in range(deg)]
    v = array(v)
    print(astr(v))

    Xv = dot(X, v)
    print(astr(Xv))

    #print(astr(v1))
    A = empty((deg, deg), dtype=object)
    
    for i in range(deg):
      for j in range(deg):
        A[i, j] = Xv[i][j]

    print(astr(A))

    #for i in GL_1[1:]:
    for i in range(1, n):
        B = A
        for j in range(i-1):
            B = dot(A, B)
        print(astr(B))


def powdot(A, n):
    B = numpy.identity(A.shape[0], dtype=A.dtype)
    for i in range(n):
        B = dot(A, B)
    return B



def fstr(v):
    vs = empty(v.shape, dtype=object)
    for idx in numpy.ndindex(v.shape):
        assert v[idx] == int(v[idx])
        vs[idx] = "%.0f"%(v[idx])
    s = str(vs)
    s = s.replace("'", "")
    #s = s.replace("-0.0000", "0.0000")
    s = s.replace("-0", "0")
    return s

def numrepr(A):
    return fstr(A) # careful !



# does not need hashable operators
def mulclose(gen, verbose=False, maxsize=None):
    ops = list(gen)
    bdy = gen
    while bdy:
        _bdy = []
        for g in bdy:
            for h in gen:
                k = dot(g, h)
                for j in ops:
                    if numpy.allclose(j, k):
                        break
                else:
                    ops.append(k)
                    _bdy.append(k)
        bdy = _bdy
        if verbose:
            print("mulclose:", len(ops))
        if maxsize and len(ops) >= maxsize:
            break
    return ops


def mulclose_fast(gen, verbose=False, maxsize=None):
    ops = list(gen)
    bdy = gen
    found = set()
    while bdy:
        _bdy = []
        for g in bdy:
            for h in gen:
                k = dot(g, h)
                s = numrepr(k)
                if s not in found:
                    found.add(s)
                    ops.append(k)
                    _bdy.append(k)
        bdy = _bdy
        if verbose:
            print("mulclose:", len(ops))
        if maxsize and len(ops) >= maxsize:
            break
    return ops


def make_group(ops):
    ops = list(ops)
    lookup = dict((numrepr(A), i) for (i, A) in enumerate(ops))
    #for A in ops:
    #    print(repr(numrepr(A)))
    items = list(range(len(ops)))
    perms = []
    for i, A in enumerate(ops):
        perm = {}
        for j, B in enumerate(ops):
            C = dot(A, B)
            k = lookup[numrepr(C)]
            perm[j] = k
        perm = Perm(perm, items)
        perms.append(perm)
    G = Group(perms, items)
    return G

def get_signature(G):
    sizes = []
    for H in G.cyclic_subgroups():
        sizes.append(len(H))
    sizes.sort()
    return sizes

    

def test():
    # 1/8 root of unity
    w = array([
        [0., 1, 0, 0],
        [0., 0, 1, 0],
        [0., 0, 0, 1],
        [-1., 0, 0, 0]])

    w7 = powdot(w, 7)
    r2 = w + w7 # root 2
    ir2 = 0.5 * r2 # 1/r2
    #print(r2)
    #print(powdot(r2, 2))

    I2 = numpy.identity(2)
    I4 = numpy.identity(4)
    I8 = numpy.identity(8)

    assert numpy.allclose(dot(r2, ir2), I4)
    
    G = mulclose([w])
    assert len(G) == 8

    Z = numpy.zeros((8, 8))
    Z[0:4, 0:4] = I4
    Z[4:8, 4:8] = -I4
    
    X = numpy.zeros((8, 8))
    X[4:8, 0:4] = I4
    X[0:4, 4:8] = I4

    G = mulclose([Z, X])
    assert len(G) == 8
    
    S = numpy.zeros((8, 8))
    S[0:4, 0:4] = I4
    S[4:6, 6:8] = I2
    S[6:8, 4:6] = -I2
    
    T = numpy.zeros((8, 8))
    T[0:4, 0:4] = I4
    T[4:8, 4:8] = w

    H = numpy.zeros((8, 8))
    H[0:4, 0:4] = ir2
    H[0:4, 4:8] = ir2
    H[4:8, 0:4] = ir2
    H[4:8, 4:8] = -ir2
    
    assert numpy.allclose(dot(H, H), I8)
    
    assert len(mulclose([X, Z])) == 8 
    assert len(mulclose([X, S])) == 32 
    assert len(mulclose([X, T])) == 128 
    assert len(mulclose([X, S, H])) == 192 

    if 0:
        G = mulclose_fast([X, T, H], maxsize=2000) # infinite group
        print(len(G))

    CZ = numpy.identity(8)
    CZ[7, 7] = -1.
    CZ[6, 6] = -1.

    CCZ = numpy.identity(8)
    CCZ[7, 7] = -1.

    G = mulclose([X, Z, CZ])
    assert len(G) == 32

    G1 = make_group(G)
    assert len(G1) == 32
    s1 = get_signature(G1)

    G2 = make_group(mulclose([X, S]))
    s2 = get_signature(G2)
    assert s1 != s2 # so these are different groups 

    assert len(mulclose([X, Z, CZ, CCZ])) == 128

    if 0:
        vals, vecs = eig(w)
        for i in range(4):
            print(vals[i], 1.j*numpy.log(vals[i]))
            print(vecs[:, i])
            print()

    print("OK")


if __name__ == "__main__":

    #test_cyclotomic()
    test()








