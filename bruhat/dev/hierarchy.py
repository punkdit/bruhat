#!/usr/bin/env python

from __future__ import print_function

import os
os.environ["SAGE_LOCAL"] = "/usr"
os.environ["SAGE_ROOT"] = '/usr/share/sagemath'
os.environ["SAGE_SRC"] = '/usr/share/sagemath/src'

for (k,v) in {
 'CONWAY_POLYNOMIALS_DATA_DIR': '/usr/share/sagemath/conway_polynomials',
 'DOT_SAGE': '/extra/simon/.sage',
 'ELLCURVE_DATA_DIR': '/usr/share/sagemath/ellcurves',
 'GAP_ROOT_DIR': '/usr/share/gap',
 'GRAPHS_DATA_DIR': '/usr/share/sagemath/graphs',
 'HOSTNAME': 'monic',
 'LOCAL_IDENTIFIER': 'monic.10722',
 'POLYTOPE_DATA_DIR': '/usr/share/sagemath/reflexive_polytopes',
 'PYTHON_EGG_CACHE': '/extra/simon/.sage/.python-eggs',
 'REALM': 'sage.math.washington.edu',
 'SAGE_BANNER': '',
 'SAGE_DATE': '2017-12-07',
 'SAGE_DISTFILES': '/usr/share/sagemath/upstream',
 'SAGE_DOC': '/usr/share/doc/sagemath',
 'SAGE_DOC_SRC': '/usr/share/sagemath/src/doc',
 'SAGE_DOT_GIT': '/usr/share/sagemath/.git',
 'SAGE_ETC': '/usr/etc',
 'SAGE_EXTCODE': '/usr/share/sagemath/ext',
 'SAGE_IMPORTALL': 'yes',
 'SAGE_INC': '/usr/include',
 'SAGE_LIB': '/usr/lib/python2.7/dist-packages',
 'SAGE_LOCAL': '/usr',
 'SAGE_LOGS': '/usr/share/sagemath/logs/pkgs',
 'SAGE_PKGS': '/usr/share/sagemath/build/pkgs',
 'SAGE_REPO_ANONYMOUS': 'git://trac.sagemath.org/sage.git',
 'SAGE_REPO_AUTHENTICATED': 'ssh://git@trac.sagemath.org:2222/sage.git',
 'SAGE_ROOT': '/usr/share/sagemath',
 'SAGE_SCRIPTS_DIR': '/usr/share/sagemath/bin',
 'SAGE_SHARE': '/usr/share/sagemath',
 'SAGE_SPKG_INST': '/usr/share/sagemath/installed',
 'SAGE_SRC': '/usr/share/sagemath/src',
 'SAGE_STARTUP_FILE': '/extra/simon/.sage/init.sage',
 'SAGE_URL': 'http://sage.math.washington.edu/sage/',
 'SAGE_VERSION': '8.1',
 'SINGULAR_SO': '/usr/lib/x86_64-linux-gnu/libsingular-Singular-4.1.0.so',
 #'SITE_PACKAGES': ['/usr/lib/python2.7/dist-packages'],
 'THEBE_DIR': '/usr/share/thebe',
 'TRAC_SERVER_URI': 'https://trac.sagemath.org',
 'UNAME': 'Linux'}.items():
    os.environ[k] = v


from sage.all_cmdline import *   # import sage library

#_sage_const_2 = Integer(2); _sage_const_1 = Integer(1);
#_sage_const_7 = Integer(7); _sage_const_4 = Integer(4)

#import sys
#print(sys.argv)

def inv(op):
    op = op.conjugate_transpose()
    op.set_immutable() # ARGH!
    return op

def mul(a, b):
    c = a*b
    c.set_immutable() # ARGH!
    return c

def mulclose(gen, mul, verbose=False, maxsize=None):
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

from bruhat.argv import argv
d = argv.get("d", 5)

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
        x, = R._first_ngens(1)
        f = x**2 - d
        assert f.is_irreducible()

        #print(K.extension.__doc__)
        K = K.extension(f, "r"+str(d))
        root_d = K.gen()

    else:

        R = PolynomialRing(K, names=('x',))
        x, = R._first_ngens(1)
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
    T[i, i] = w3**i
    for j in range(d):
        H[i, j] = w**(i*j) / root_d


for op in [X, Z, S, T, H]:
    op.set_immutable()
    assert op*inv(op) == I, op

assert X**0 == I
assert X**d == I
assert Z**0 == I
assert Z**d == I

C1 = mulclose([X, Z, wI], mul)
print("|C1| =", len(C1)) # == 16

C2 = mulclose([X, S, H], mul)
print("|C2| =", len(C2)) # == 192


def is_cliff(A):
    Ai = inv(A)
    for g in [X, Z]:
        h = A*g*Ai
        h.set_immutable()
        if h not in C1:
            return False
    return True

#for A in C2:
#    assert is_cliff(A)
assert S not in C1
assert is_cliff(S)

def is_third_level(A):
    Ai = inv(A)
    for g in [X, Z]:
        h = A*g*Ai
        h.set_immutable()
        if not is_cliff(h):
            return False
    return True

assert not is_cliff(T)
assert is_third_level(T)

if 1:
    C3 = set(C2)
    for a in C2:
        aT = a*T
        for b in C2:
            A = aT*b
            A.set_immutable()
            C3.add(A)
    print("|C3| =", len(C3))
    #for A in C3:
    #    assert is_third_level(A)

    def src(a):
        return set([b for b in C3 if mul(a,b) in C3])

    def tgt(a):
        return set([b for b in C3 if mul(b,a) in C3])

    srcs = []
    for b in C3:
        src_b = src(b)
        if src_b not in srcs:
            print("|src_b| = ", len(src_b))
            srcs.append(src_b)
        if d==2 and len(srcs)==4: # there is only 4 of these to be found
            break

    obs = list(srcs)
    tgts = []
    for b in C3:
        tgt_b = tgt(b)
        if tgt_b not in obs:
            obs.append(tgt_b)
        if tgt_b not in tgts:
            print("|tgt_b| = ", len(tgt_b))
            tgts.append(tgt_b)
        if d==2 and len(tgts)==4: # there is only 4 of these to be found
            break

    done = False
    while not done:
        done = True
        print("obs:", len(obs))

        obs.sort(key = len, reverse=True)
        obs1 = list(obs)
        for s in obs:
          for t in obs:
            st = s.intersection(t)
            a = ' '
            if st not in obs1:
                a = '*'
                obs1.append(st)
                done = False
            print("%4d"%(len(st)), end=a)
          print()
        obs = obs1
        obs.sort(key = len, reverse=True)

    print("obs:", len(obs), [len(ob) for ob in obs])
    print(set(C2) == obs[-1])
    for ob in obs:
        #print(len(ob), end=" ")
        iob = set([inv(a) for a in ob])
        print(iob in obs, end=" ")
        #iob = [ia for a in iob if inv(a) in ob]
        #print((len(ob), len(iob)), end=" ")

    print()




#A = MS.matrix([1, 0, 0, 0, -1, 0, 0, 0, -1])
#print(is_cliff(A)) # nope
#print(is_third_level(A)) # nope

# --------------- two qudits -------------- #

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
    assert op*inv(op) == II

C1_2 = mulclose(gen, mul)
print("|C1(2)| =", len(C1_2))

def is_cliff_2(A):
    Ai = inv(A)
    for g in [XI, IX, ZI, IZ]:
        h = A*g*Ai
        h.set_immutable()
        if h not in C1_2:
            return False
    return True

def is_third_level_2(A):
    Ai = inv(A)
    for g in [XI, IX, ZI, IZ]:
        h = A*g*Ai
        h.set_immutable()
        if not is_cliff_2(h):
            return False
    return True


def block_diag(blocks):
    A = copy(MS_2.zero_matrix())
    
    assert len(blocks)==d
    for i in range(d):
        A[d*i:d*(i+1), d*i:d*(i+1)] = blocks[i]
    
    A.set_immutable()
    return A

if 0:
    S_2 = block_diag([Z**(i**2) for i in range(d)])
    
    #print(S_2)
    print(S_2 in C1_2)
    print(is_cliff_2(S_2))
    print(is_third_level_2(S_2))


