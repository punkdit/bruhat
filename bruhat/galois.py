#!/usr/bin/env python3

import numpy
from numpy import dot, array, empty
from numpy.linalg import eigvals, eig

from bruhat.element import PolynomialRing, Z, cyclotomic
from bruhat.argv import argv


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
    B = numpy.identity(A.shape, dtype=A.dtype)
    for i in range(n):
        B = A*B
    return B


def test():
    w = array([
        [0., 1, 0, 0],
        [0., 0, 1, 0],
        [0., 0, 0, 1],
        [-1., 0, 0, 0]])

    print(powdot(w, 7))
    return

    vals, vecs = eig(w)
    for i in range(4):
        print(vals[i], 1.j*numpy.log(vals[i]))
        print(vecs[:, i])
        print()



if __name__ == "__main__":

    #test_cyclotomic()
    test()








