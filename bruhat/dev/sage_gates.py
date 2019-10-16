#!/usr/bin/env python2

from __future__ import print_function

from bruhat.dev.sage_env import *

from bruhat.argv import argv
d = argv.get("d", 2) # qudit dimension


def dag(op):
    op = op.conjugate_transpose()
    op.set_immutable() # ARGH!
    return op

def mul(a, b):
    c = a*b
    c.set_immutable() # ARGH!
    return c

def mulclose(gen, mul=mul, verbose=False, maxsize=None):
    for A in gen:
        A.set_immutable()
    els = set(gen)
    bdy = list(els)
    changed = True
    while bdy:
        #if verbose:
        #    print "mulclose:", len(els)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = mul(A, B)
                if C not in els:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return els

if d==2:
    K = CyclotomicField(8)
    e8 = K.gen()
    e4 = e8**2
    w = e4**2 
    w2 = e4
    w3 = e8
    assert w==-1
    root_d = e8 + e8.conjugate()
    assert root_d**2 == 2

elif d==3:
    K = CyclotomicField(9)
    e9 = K.gen()
    e3 = e9**3
    w = e3
    w2 = e3
    w3 = e9

    R = PolynomialRing(K, names=('x',))
    x, = R._first_ngens(1)
    f = x**2 - d
    assert f.is_irreducible()

    #print(K.extension.__doc__)
    K = K.extension(f, "r"+str(d))
    root_d = K.gen()

else:
    K = CyclotomicField(d)
    w = K.gen()
    w2 = w
    w3 = w

    if d%4 == 3:
        R = PolynomialRing(K, names=('x',))
        #x, = R._first_ngens(1)
        x = R.gen()
        f = x**2 - d
        assert f.is_irreducible()

        #print(K.extension.__doc__)
        K = K.extension(f, "r"+str(d))
        root_d = K.gen()

    else:

        R = PolynomialRing(K, names=('x',))
        #x, = R._first_ngens(1)
        x = R.gen()
        f = x**2 - d
        assert not f.is_irreducible()
        g = f.factor()
        root_d = g[1][0][0]

assert root_d**2 == d


MS = MatrixSpace(K, d)
MS_2 = MatrixSpace(K, d**2)

I = MS.identity_matrix()
II = MS_2.identity_matrix()

wI = w2*I
wI.set_immutable()

X = MS.matrix()
Z = MS.matrix()
S = MS.matrix()
T = MS.matrix()
H = MS.matrix()
for i in range(d):
    X[i, (i+1)%d] = 1
    Z[i, i] = w**i
    S[i, i] = w2**(i**2)
    T[i, i] = w3**(i**3)
    for j in range(d):
        H[i, j] = w**(i*j) / root_d


for op in [X, Z, S, T, H]:
    op.set_immutable()
    assert op*dag(op) == I, op

assert X**0 == I
assert X**d == I
assert Z**0 == I
assert Z**d == I


def control(A):
    ms1 = A.matrix_space()
    m, n = ms1.dims()
    #I = ms1.identity_matrix()
    ms2 = MatrixSpace(K, m*d)
    CA = ms2.matrix()
    for i in range(d):
        CA[ d*i : d*(i+1), d*i : d*(i+1) ] = A**i
    return CA


def block_diag(blocks):
    #A = copy(MS_2.zero_matrix())
    A = blocks[0]
    ms1 = A.matrix_space()
    m, n = ms1.dims()
    N = len(blocks)
    ms2 = MatrixSpace(K, m*N)

    AA = ms2.matrix()
    for i in range(N):
        AA[m*i:m*(i+1), m*i:m*(i+1)] = blocks[i]
    AA.set_immutable()
    return AA


        
# -------------- two qudits ---------------

II = I.tensor_product(I)
HH = H.tensor_product(H)
IH = I.tensor_product(H)
HI = H.tensor_product(I)
CX = control(X)
CZ = control(Z)
CH = control(H)
CS = control(S)
CT = control(T)

assert HH * HH == II
#print(HH * CX * HH)
#print(IH * CX * IH)
#print(HH * CZ * HH)
#print(HH)

XI = X.tensor_product(I)
IX = I.tensor_product(X)
ZI = Z.tensor_product(I)
IZ = I.tensor_product(Z)
wII = w2*II

gen = [XI, IX, ZI, IZ, wII]
for op in gen:
    op.set_immutable()
    #print()
    #print(op)
    assert op*dag(op) == II

if d==2:
    ms = MatrixSpace(K, 4)
    CX2 = ms.matrix([
        1,0,0,0,
        0,0,0,1,
        0,0,1,0,
        0,1,0,0])
    
    SWAP = CX2*CX*CX2
    #gen = [XI, IX, CX, CX2]
    #G = mulclose(gen)
    #print(len(G))


XII = XI.tensor_product(I)
IXI = IX.tensor_product(I)
IIX = II.tensor_product(X)
ICX = I.tensor_product(CX)
CXI = CX.tensor_product(I)



