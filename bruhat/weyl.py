#!/usr/bin/env python

"""
Represent any Weyl group as a permutation group of a root system.

See: 
    Humphreys, p63-66
    https://en.wikipedia.org/wiki/Root_system

"Root Systems and Generalized Associahedra"
Sergey Fomin and Nathan Reading
https://arxiv.org/pdf/math/0505518.pdf

"""


import sys, os
import string
from fractions import Fraction
from operator import mul

import numpy

from action import Perm, Group, mulclose
from gelim import solve, array, identity, dot, shortstr, eq, dotx, kernel
from gelim import Subspace
from argv import argv


def cross(itemss):
    if len(itemss)==0:
        yield ()
    else:
        for head in itemss[0]:
            for tail in cross(itemss[1:]):
                yield (head,)+tail


def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r
 

def choose(m, n):
    return factorial(m) // factorial(n)


def rdot(a, b):
    assert len(a)==len(b)
    r = sum(ai*bi for (ai,bi) in zip(a, b))
    return r


def rnorm2(a):
    r = sum(ai*ai for ai in a)
    return r


def rscale(v, a):
    assert Fraction(v)==v
    assert type(a) is tuple
    a = tuple(v*ai for ai in a)
    return a



def mulclose_pri(els, verbose=False, maxsize=None):
    "multiplicative closure; short words first"
    els = set(els)
    changed = True
    while changed:
        if verbose:
            print "mulclose:", len(els)
        changed = False
        _els = list(els)
        pairs = [(g, h) for g in _els for h in _els]
        pairs.sort(key = lambda (g,h) : len(g.word+h.word))

        for A, B in pairs:
            C = A*B 
            if C not in els:
                els.add(C)
                if maxsize and len(els)>=maxsize:
                    return list(els)
                changed = True
    return els 


class Weyl(object):
    def __init__(self, roots, gen, simple=None, name=None, check=True):
        assert self.__class__ != Weyl
        self.roots = roots # list of tuples
        self.gen = gen # list of Weyl group generators (Perm's)
        self.identity = Perm.identity(roots, '')
        self.simple = simple
        for g in gen:
            assert g*g == self.identity

        self.n = len(gen)
        self.name = name

        if self.simple is None:
            self.build_simple()

        assert self.simple is not None
        if check:
            self.check_simple()

    def __str__(self):
        return "%s(%d)"%(self.__class__.__name__, self.n)

    def generate(self):
        gen = self.gen
        n = len(gen)
        names = string.uppercase
        assert n<=len(names)
        for i, g in enumerate(gen):
            g.word = names[i]
            #print "%s:"%g.word,
            #print g.str()
        #print
        roots = self.roots
        weyl = mulclose_pri([self.identity]+gen)
        weyl = list(weyl)
        weyl.sort(key = lambda g : (len(g.word), g.word))
        return weyl

    def build_simple(self):
        #print "build_simple"
        gen = self.gen
        roots = self.roots
        n = len(gen)
        simple = []
        for g in gen:
          for root in roots:
            nroot = rscale(-1, root)
            if g(root) != nroot:
                continue
            assert nroot not in simple

            #print "root:", shortstr(root)
            if not simple:
                simple.append(root)
                break

            signs = [rdot(root, _root) for _root in simple]
            #print "signs:", signs
            if max(signs) == min(signs) == 0:
                simple.append(root)
            elif max(signs) > 0:
                assert min(signs) >= 0
                simple.append(nroot)
            else:
                simple.append(root)
            break

        assert len(simple) == n, "only found %s simple roots, expected %d" % (len(simple), n)
        self.simple = simple

    @classmethod
    def buildfrom_simple(cls, roots, simple, **kw):
        "use _generators from reflections of simple roots"
        for root in simple:
            assert root in roots
        n = len(roots[0])
        idxs = range(n)
        lookup = dict((root, i) for (i, root) in enumerate(roots))
        gen = []
        for alpha in simple:
            #print "alpha:", alpha
            r0 = sum(alpha[i]*alpha[i] for i in idxs)
            perm = []
            for root in roots:
                #print "    root:", root
                r = sum(alpha[i]*root[i] for i in idxs)
                _root = tuple(root[i] - 2*Fraction(r, r0)*alpha[i] for i in idxs)
                #print "    _root:", _root
                perm.append(lookup[_root])
            perm = Perm(perm, roots)
            gen.append(perm)
        return cls(roots, gen, simple, **kw)

    def check_simple(self):
        """
        A set of simple roots for a root system \Phi is
        a set of roots that form a basis for the Euclidean
        space spanned by \Phi with the special property that each
        root has components with respect to this basis that are
        either all nonnegative or all nonpositive.
        
        https://en.wikipedia.org/wiki/E8_(mathematics)#Simple_roots
        """

        simple = self.simple
        n = len(simple)
        U = array(simple)
        #print "simple:"
        #print U
        V = solve(U, identity(len(U)), check=True)
        #print "inverse:"
        #print V
        #print "UV:"
        #print dot(U, V)

        roots = self.roots
        gen = self.gen
        n = len(simple)
        for i in range(n):
            root = simple[i]
            assert gen[i](root) == rscale(-1, root)
            #print root, gen[i](root)

        for i in range(n):
            for j in range(i+1, n):
                assert rdot(simple[i], simple[j]) <= 0

        for root in roots:
            nroot = rscale(-1, root)
            for g in gen:
                if g(root) == nroot:
                    assert root in simple or rscale(-1, root) in simple

            if root in simple or nroot in simple:
                continue
            a = dot(root, V)
            #print "root:", shortstr(root, a)
            pos = (a>=0).astype(numpy.int)
            neg = (a<=0).astype(numpy.int)
            #assert pos.sum() == n or neg.sum() == n
            if pos.sum() != n and neg.sum() != n:
                assert 0, "FAIL"
        #print "check_simple: OK"

    def matrix(self, desc=None):
        "coxeter matrix"
        m = {}
        gen = self.gen
        n = len(gen)
        for i in range(n):
          for j in range(i+1, n):
            key = (i, j)
            if desc:
                key = (desc[i], desc[j])
            gi = gen[i]
            gj = gen[j]
            g = gigj = gi*gj
            k = 1
            while 1:
                if g == self.identity:
                    if k > 2:
                        m[key] = k
                    break
                g = gigj * g
                k += 1
                assert k<10
        return m

    @classmethod
    def build_A(cls, n, **kw):
        roots = []
        simple = []
        lookup = {}
        for i in range(n+1):
          for j in range(n+1):
    
            if i==j:
                continue
            root = [0]*(n+1)
            root[i] = 1
            root[j] = -1
    
            root = tuple(root)
            lookup[root] = len(roots)
            roots.append(root)
            if j==i+1:
                simple.append(root)
    
        #assert len(pos_roots) == choose(n+1, 2) # number of positive roots
        assert len(lookup) == len(roots)
        assert len(roots) == n*(n+1)
    
        gen = []
        for i in range(n):
            perm = []
            for idx, root in enumerate(roots):
                _root = list(root)
                _root[i], _root[i+1] = _root[i+1], _root[i] # swap i, i+1 coords
                jdx = lookup[tuple(_root)]
                perm.append(jdx)
            perm = Perm(perm, roots)
            gen.append(perm)
        return Weyl_A(roots, gen, simple, name="A_%d"%n, **kw)

    @classmethod
    def build_B(cls, n, k=1, **kw):
    
        roots = []
        lookup = {}

        # long roots
        for i in range(n):
          for j in range(i+1, n):
    
            for a in [-1, 1]:
             for b in [-1, 1]:
                root = [0]*n
                root[i] = a
                root[j] = b
                root = tuple(root)
                lookup[root] = len(roots)
                roots.append(root)

        # short roots
        for i in range(n):
          for a in [-1, 1]:
            root = [0]*n
            root[i] = a*k
            root = tuple(root)
            lookup[root] = len(roots)
            roots.append(root)
        assert len(lookup) == len(roots)
        assert len(roots) == 2*n**2
    
        gen = []
        for i in range(n-1):
            perm = []
            for idx, root in enumerate(roots):
                _root = list(root)
                _root[i], _root[i+1] = _root[i+1], _root[i] # swap i, i+1 coords
                jdx = lookup[tuple(_root)]
                perm.append(jdx)
            perm = Perm(perm, roots)
            gen.append(perm)

        perm = []
        for idx, root in enumerate(roots):
            _root = list(root)
            _root[n-1] = -_root[n-1]
            jdx = lookup[tuple(_root)]
            perm.append(jdx)
        perm = Perm(perm, roots)
        gen.append(perm)

        return Weyl_B(roots, gen, name="B_%d"%n, **kw)

    @classmethod
    def build_C(cls, n, **kw):
        # um,,, fix the roots here..
        return Weyl_C.build_B(n, k=2, **kw)

    @classmethod
    def build_D(cls, n, **kw):

        assert n>=2

        roots = []
        lookup = {}

        for i in range(n):
          for j in range(i+1, n):
    
            for a in [-1, 1]:
             for b in [-1, 1]:
                root = [0]*n
                root[i] = a
                root[j] = b
                root = tuple(root)
                lookup[root] = len(roots)
                roots.append(root)

        assert len(lookup) == len(roots)
        assert len(roots) == 2*n*(n-1)

        gen = []
        for i in range(n-1):
            perm = []
            for idx, root in enumerate(roots):
                _root = list(root)
                _root[i], _root[i+1] = _root[i+1], _root[i] # swap i, i+1 coords
                jdx = lookup[tuple(_root)]
                perm.append(jdx)
            perm = Perm(perm, roots)
            gen.append(perm)

        perm = []
        for idx, root in enumerate(roots):
            _root = list(root)
            _root[n-1], _root[n-2] = -_root[n-2], -_root[n-1] # swap & negate last two coords
            jdx = lookup[tuple(_root)]
            perm.append(jdx)
        perm = Perm(perm, roots)
        gen.append(perm)

        return Weyl_D(roots, gen, name="D_%d"%n, **kw)

    @classmethod
    def build_E8(cls, check=False):
        D8 = cls.build_D(8, check=check)
        half = Fraction(1, 2)
        roots = list(D8.roots)
        for signs in cross([(-1, 1)]*8):
            r = reduce(mul, signs)
            if r != 1:
                continue
            root = tuple(sign*half for sign in signs)
            roots.append(root)

        assert len(roots) == 240

        simple = []
        for i in range(6):
            root = [0]*8
            root[i] = 1
            root[i+1] = -1
            simple.append(tuple(root))

        root = [0]*8
        root[i] = 1
        root[i+1] = 1
        simple.append(tuple(root))

        root = [-half]*8
        simple.append(tuple(root))

        return Weyl_E8.buildfrom_simple(roots, simple, name="E_8", check=check)

    @classmethod
    def build_E7(cls, check=False):
        E8 = cls.build_E8(check=check)
        idxs = range(8)
        # delete one root:
        root0 = E8.simple[0]
        roots = []
        for root in E8.roots:
            if sum(root0[i]*root[i] for i in idxs)==0:
                roots.append(root)
        assert len(roots)==126
        simple = [rscale(-1, roots[0])] + [root for root in E8.simple if root in roots]
        return Weyl_E7.buildfrom_simple(roots, simple, name="E_7", check=check)

    @classmethod
    def build_E6(cls, check=False):
        E8 = cls.build_E8(check=check)
        idxs = range(8)
        # delete two roots:
        root0 = E8.simple[0]
        root1 = E8.simple[1]
        roots = []
        for root in E8.roots:
            if sum(root0[i]*root[i] for i in idxs)==0 and \
                sum(root1[i]*root[i] for i in idxs)==0:
                roots.append(root)
        assert len(roots)==72
        simple = [(0, 0, 0, -1, 0, 0, 0, 1)]
        simple += [root for root in E8.simple if root in roots]
        return Weyl_E6.buildfrom_simple(roots, simple, name="E_6", check=check)

    @classmethod
    def build_F4(cls, **kw):
        roots = []
        idxs = range(4)
        for root in cross([(-1,0,1)]*4):
            d = sum(root[i]**2 for i in idxs)
            if d==1 or d==2:
                roots.append(root)
        half = Fraction(1, 2)
        for root in cross([(-half,0,half)]*4):
            d = sum(root[i]**2 for i in idxs)
            if d==1 or d==2:
                roots.append(root)
        assert len(roots)==48
        simple = [
            (1, -1, 0, 0),
            (0, 1, -1, 0),
            (0, 0, 1, 0),
            (-half, -half, -half, -half)]
        return Weyl_F.buildfrom_simple(roots, simple, name="F_4", **kw)

    @classmethod
    def build_G2(cls, **kw):
        roots = []
        for root in cross([(-2, -1, 0, 1, 2)]*3):
            if sum(root) != 0:
                continue
            d = sum(root[i]**2 for i in range(3))
            if d==2 or d==6:
                roots.append(root)
        assert len(roots)==12
        simple = [(1, -1, 0), (-1, 2, -1)]
        return Weyl_G.buildfrom_simple(roots, simple, name="G_2", **kw)

    def prod_sk(self, sk):
        g = self.identity
        gen = self.gen
        for i in sk:
            g = g * gen[i]
        return g

    def longest_element(self):
        g = self.identity
        for sk in self.iter_sk():
            g = g * self.prod_sk(sk)
        return g


"""
Calculation of longest element:
https://arxiv.org/pdf/1108.1048v2.pdf Table 1,

Also Humphreys, "when is -I in W?":
http://people.math.umass.edu/~jeh/pub/longest.pdf

http://mathoverflow.net/questions/54926/longest-element-of-weyl-groups
"""

class Weyl_A(Weyl):
    def iter_sk(self):
        n = self.n
        for k in range(1, n+1):
            sk = []
            i = n-k
            while i<=n-1:
                sk.append(i)
                i += 1
            yield sk


class Weyl_B(Weyl):
    def iter_sk(self):
        n = self.n
        for k in range(1, n+1):
            sk = []
            i = n-k
            while i<=n-1:
                sk.append(i)
                i += 1
            i -= 1
            while i>n-k:
                i -= 1
                sk.append(i)
            #print "k=%d, sk=%s"%(k, sk)
            yield sk


class Weyl_C(Weyl_B):
    pass


class Weyl_D(Weyl):
    def iter_sk(self):
        n = self.n
        assert n>2
        yield [n-1]
        yield [n-2]
        for k in range(3, n+1):
            sk = []
            i = n-k
            while i<=n-3:
                sk.append(i)
                i += 1
            sk.append(n-1)
            i = n-2
            while i>=n-k:
                sk.append(i)
                i -= 1
            #print "k=%d, sk=%s"%(k, sk)
            yield sk


class Weyl_E6(Weyl):
    diagram = """
          3
          |
    0--1--2--4--5
    """
    def iter_sk(self):
        n = self.n
        D_5 = Weyl.build_D(5, check=False)
        for sk in D_5.iter_sk():
            yield sk
        sk = [i-1 for i in [6, 5, 3, 4, 2, 1, 3, 2, 5, 3, 4, 6, 5, 3, 2, 1]]
        yield sk


class Weyl_E7(Weyl):
    diagram = """
          4
          |
    1--2--3--5--6--0
    """
    def iter_sk(self):
        n = self.n
        D_5 = Weyl.build_D(5, check=False)
        for sk in D_5.iter_sk():
            sk = [i+1 for i in sk]
            yield sk
        sk = [6, 5, 3, 4, 2, 1, 3, 2, 5, 3, 4, 6, 5, 3, 2, 1]
        yield sk
        sk = [7, 6, 5, 3, 4, 2, 1, 3, 2, 5, 3, 4, 
            6, 5, 3, 2, 1, 7, 6, 5, 3, 4, 2, 3, 5, 6, 7]
        sk = [i%7 for i in sk] # fix numbers
        yield sk


class Weyl_E8(Weyl):
    diagram = """
          5
          |
    7--6--4--3--2--1--0
    """
    def iter_sk(self):
        tx = {1:7, 2:6, 3:4, 4:5, 5:3, 6:2, 7:1, 8:0}
        n = self.n
        D_5 = Weyl.build_D(5, check=False)
        for sk in D_5.iter_sk():
            sk = [tx[i+1] for i in sk]
            yield sk
        sk = [6, 5, 3, 4, 2, 1, 3, 2, 5, 3, 4, 6, 5, 3, 2, 1]
        sk = [tx[i] for i in sk]
        yield sk
        sk = [7, 6, 5, 3, 4, 2, 1, 3, 2, 5, 3, 4, 
            6, 5, 3, 2, 1, 7, 6, 5, 3, 4, 2, 3, 5, 6, 7]
        sk = [tx[i] for i in sk]
        yield sk
        sk = [8, 7, 6, 5, 3, 4, 2, 1, 3, 2, 5, 3, 4, 6, 5, 3, 2, 1, 7, 6, 5, 3, 4,
            2, 3, 5, 6, 7, 8, 7, 6, 5, 3, 4, 2, 1, 3, 2, 5, 3, 4, 6, 5, 3, 2, 1,
            7, 6, 5, 3, 4, 2, 3, 5, 6, 7, 8]
        sk = [tx[i] for i in sk]
        yield sk



class Weyl_F(Weyl):
    def iter_sk(self):
        B_3 = Weyl.build_B(3)
        for sk in B_3.iter_sk():
            yield sk
        sk = [i-1 for i in [4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4]]
        yield sk


class Weyl_G(Weyl):
    def iter_sk(self):
        yield [0]
        yield [1, 0, 1, 0, 1]



def representation(G, verbose=False):

    gen = G.gen
    simple = G.simple
    #print "simple:", simple

    # Simple vectors in the Euclidean (standard) basis.
    S = array(simple)
    print "S:"
    print shortstr(S)
    I = identity(G.n)

    # T maps vectors in the Euclidean basis to vectors in the Simple basis.
    Tt = solve(S, identity(G.n), check=True)
    T = Tt.transpose()
    print "T:"
    print shortstr(T)
    print '--'*10

    assert eq(dot(S, Tt), I)

    Ws = [] # weight spaces
    As = [] # As for generators

    for g in gen:

        # A is the matrix representation of g in the simple root basis.

        #print "A:", [g(s) for s in simple]
        rows = []
        for s in simple:
            r = g(s)
            r = dot(T, r)
            #print shortstr(r)
            rows.append(r)
        A = array(rows).transpose()
        assert eq(dot(A, A), identity(len(A))) # reflection
        As.append(A)
        print "A:"
        print shortstr(A)

        print "St*A*T:"
        print shortstr(dotx(S.transpose(), A, T))

        # W is the fixed space of A, expressed as vectors in the Euclidean basis.
        # These vectors (rows of W) lie on the "mirror" that is A.
        W = []
        for root in simple:
            groot = g(root)
            if groot == rscale(-1, root):
                continue

            w = array(root) + array(groot) # easy!
            W.append(w)
            assert numpy.abs(w).sum()

            #u = solve(S.transpose(), w)
            #v = dot(S.transpose(), u)
            #assert eq(v, w)

        W = array(W)
        Ws.append(W)
        print "W:"
        print shortstr(W)
        B = dotx((A-identity(len(A))), T, W.transpose())
        #print "A*T*Wt:"
        #print shortstr(dotx(A, T, W.transpose()))
        #print "(A-I)*T*Wt:"
        #print shortstr(B)
        assert numpy.abs(B).sum() == 0

        #B = dotx(S.transpose(), A, T)
        #print "B:"
        #print shortstr(B)

    print '--'*10

    Ws = [Subspace(W) for W in Ws]
    for W in Ws:
        print "W:"
        print W

    print '--'*10

    # not needed...
    # W0 = Subspace(S)
    # Ws = [W.intersect(W0) for W in Ws]

    # Each of these weights will lie on all but one of the A mirrors
    weights = []

    for i in range(G.n):

        Ws1 = list(Ws)
        Ws1.pop(i)

        W = Ws1[0]
        for j in range(1, len(Ws1)):
            W = W.intersect(Ws1[j])

        print "weight:"
        print W

        W = W.W # get the numpy array
        assert len(W)==1
        W.shape = W.shape[1],

        A = As[i]
        B = dotx(S.transpose(), A, T)
        w = W.transpose()
        Bw = dot(B, w)
        assert not eq(Bw, w)
        v = Bw - w # we will scale w so that this delta will be equal to S[i]
        v = v.transpose() # row vector
        s = S[i]
        r = None
        for j in range(len(v)):
            assert (v[j]==0) == (s[j]==0)
            if v[j]==0:
                continue
            assert r is None or r==s[j]/v[j]
            r = s[j]/v[j]

        w = r*w
        #print shortstr(w)
        v = dot(B, w) - w
        assert eq(v, s)

        weights.append(w)

    W = array(weights)
    print "W:"
    print shortstr(W)
    Wt = W.transpose()
    St = S.transpose()

    # What combination of weights gives the (simple) roots?
    U = solve(Wt, St)
    print "U:"
    print shortstr(U)

    # What combination of roots gives the weights?
    V = solve(St, Wt)
    print "V:"
    print shortstr(V)

    print '--'*10

    #for A in As:
    weyl = G.generate() # entire weyl group, sorted by word length
    print "weyl:", len(weyl)
    I = weyl[0]
    assert I.word == '' # identity

    rows = []
    for s in simple:
        r = I(s)
        r = dot(T, r)
        #print shortstr(r)
        rows.append(r)
    A0 = array(rows).transpose()

    for g in weyl:
        rows = []
        for s in simple:
            r = g(s)
            r = dot(T, r)
            #print shortstr(r)
            rows.append(r)
        A = array(rows).transpose()

        #print "A:"
        #print shortstr(A)
        print g.word or 'I'
        D = dot(A-A0, V).transpose()
        print shortstr(D)



def find_negi(G):

    roots = G.roots
    simple = G.simple
    print "roots:", len(roots)
#    print roots
#    print

    weyl = G.generate()

    print "weyl:",
    for w in weyl:
        print w.word,
    print
    print "weyl:", len(weyl)

    nroots = [rscale(-1, root) for root in roots]
    n = len(roots)
    N = len(weyl)

    for i in range(N):
      for j in range(i+1, N):
        X = weyl[i]
        Z = weyl[j]
        XZ = X*Z
        ZX = Z*X

        for root in simple:
            if XZ(root) != rscale(-1, ZX(root)):
                break
        else:
            if len(X.word)==len(Z.word) or 0:
                x, z = X.word, Z.word
                print "%s*%s = -%s*%s" % (x, z, z, x)
                return


    return

    for g in weyl:
        for i in range(n):
            if g(roots[i]) != nroots[i]:
                break
        else:
            print "%s = -I" % g.word



def test_monoid(G):
    "build the Coxeter-Bruhat monoid, represented as a monoid of functions."
    "this representation is not faithful."

    roots = G.roots
    print "roots:", len(roots)

    weyl = G.generate()
    print "weyl:",
    for w in weyl:
        print w.word,
    print
    print "weyl:", len(weyl)

    r0 = roots[0]
    bdy = set([r0])
    seen = set()
    identity = dict((r, r) for r in roots)
    perms = [dict(identity) for g in gen]
    while bdy:

        _bdy = set()
        seen.update(bdy)
        for r0 in bdy:
          for i in range(n):
            g = gen[i]
            perm = perms[i]
            r1 = g*r0
            assert perm.get(r0) == r0
            if r1 not in seen:
                perm[r0] = r1
                _bdy.add(r1)

        bdy = _bdy

    gen = [Perm(perms[i], roots, gen[i].word) for i in range(len(perms))]
    identity = Perm(identity, roots)
    for g in gen:
        print g.str()
        assert g*g == g
    
    monoid = mulclose_pri([identity]+gen)
    print "monoid:", len(monoid)
    #monoid = mulclose(monoid)
    #print "monoid:", len(monoid)

    desc = "ABCDEFG"[:n]

    def txlate(word):
        "monoid to weyl"
        g = identity
        for c in word:
            g = g * G.gen[desc.index(c)]
        return g

    monoid = list(monoid)
    monoid.sort(key = lambda g : (len(g.word), g.word))
    for g in monoid:
        tgt = list(set(g.perm.values()))
        g = txlate(g.word)
        print "%6s"%g.word, len(tgt)
    print

    return

    m = G.matrix(desc)

    from bruhat import BruhatMonoid
    monoid = BruhatMonoid(desc, m, bruhat=True, build=True)
    #for w in monoid.words:
    #    print w,
    #print

    def txlate(word):
        "_translate to function monoid"
        g = identity
        for c in word:
            g = g * gen[desc.index(c)]
        return g

    lookup = {}
    for w0 in monoid.words:
        g = txlate(w0)
        w1 = lookup.get(g)
        if w1 is not None:
            print w0, "=", w1
        else:
            lookup[g] = w0
            print w0

    for w0 in monoid.words:
      for w1 in monoid.words:
        w2 = monoid.mul[w0, w1]
        if txlate(w0)*txlate(w1) == txlate(w2):
            pass
        else:
            print "%r*%r = %r" % (w0, w1, w2),
            print " ****************** FAIL"



def test(n):


    G = Weyl.build_A(n)

#    print "%d roots" % len(G.roots)
#    print G.roots
#    for g in G.gen:
#        print [g(root) for root in G.roots]

    gen = G.gen
    for i in range(n):
      for j in range(n):
        gi = gen[i]
        gj = gen[j]
        if i==j:
            assert gi*gj == G.identity
        elif abs(i-j)==1:
            assert gi*gj != G.identity
            assert (gi*gj)**2 != G.identity
            assert (gi*gj)**3 == G.identity
        else:
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity

    if n < 5:
        assert len(mulclose(G.gen)) == factorial(n+1)

    # ---------------------------------------------------------

    G = Weyl.build_B(n)

    gen = G.gen
    for g in gen:
        assert g*g == G.identity

    for i in range(n-1):
      for j in range(i+1, n-1):
        gi = gen[i]
        gj = gen[j]
        if abs(i-j)==1:
            assert gi*gj != G.identity
            assert (gi*gj)**2 != G.identity
            assert (gi*gj)**3 == G.identity
        else:
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity
        if i < n-1:
            gj = gen[n-1]
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity

    if n>2:
        gi = gen[n-2]
        gj = gen[n-1]
        assert (gi*gj) != G.identity
        assert (gi*gj)**2 != G.identity
        assert (gi*gj)**3 != G.identity
        assert (gi*gj)**4 == G.identity

    if n < 5:
        assert len(mulclose(G.gen)) == (2**n)*factorial(n)

    # ---------------------------------------------------------
    G = Weyl.build_C(n)

    # ---------------------------------------------------------

    if n<3:
        return # <---------------------- return

    G = Weyl.build_D(n)

    gen = G.gen
    for i in range(n-1):
        gi = gen[i]
        for j in range(i+1, n-1):
            gj = gen[j]
            if abs(i-j)==1:
                assert gi*gj != G.identity
                assert (gi*gj)**2 != G.identity
                assert (gi*gj)**3 == G.identity
            else:
                assert gi*gj != G.identity
                assert (gi*gj)**2 == G.identity

        gj = gen[n-1]
        if i < n-3 or i==n-2:
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity
        elif i==n-3:
            assert gi*gj != G.identity
            assert (gi*gj)**2 != G.identity
            assert (gi*gj)**3 == G.identity

    if n < 5:
        assert len(mulclose(G.gen)) == (2**(n-1))*factorial(n)


def test_longest_element():

    G = Weyl.build_A(2)
    g = G.longest_element()
    A, B = G.gen
    assert g == A*B*A

    G = Weyl.build_A(3)
    g = G.longest_element()
    A, B, C = G.gen
    assert g == B*C*B*A*B*C

    w = A*C*B
    assert g == (w**2)
    return

    while w != g:
        w = w*A*C*C
        print "."

    G = Weyl.build_B(2)
    g = G.longest_element()
    A, B = G.gen
    assert g == B*A*B*A

    G = Weyl.build_B(3)
    g = G.longest_element()
    A, B, C = G.gen
    assert g == C*A*B*C*A*B*C*A*B

    G = Weyl.build_D(4)
    g = G.longest_element()
    A, B, C, D = G.gen
    assert g == B*A*D*B*A*D*C*B*D*A*B*C

    G = Weyl.build_D(6)
    g = G.longest_element()
    for root in G.roots:
        assert g(root) == rscale(-1, root)

    G = Weyl.build_E7()
    g = G.longest_element()
    for root in G.roots:
        assert g(root) == rscale(-1, root)

    G = Weyl.build_E8()
    g = G.longest_element()
    for root in G.roots:
        assert g(root) == rscale(-1, root)

    G = Weyl.build_F4()
    g = G.longest_element()
    for root in G.roots:
        assert g(root) == rscale(-1, root)

    G = Weyl.build_G2()
    g = G.longest_element()
    for root in G.roots:
        assert g(root) == rscale(-1, root)

    print "OK"
    

def main():

    n = argv.n
    if n:
        test(n)

    elif argv.test:
        for n in range(2, 6):
            test(n)
        test_longest_element()


    if argv.test_longest_element:
        test_longest_element()


    G = None

    if argv.E_8:
        G = Weyl.build_E8()
        assert G.matrix() == {
            (0, 1):3, (1, 2):3, (2, 3):3, (3, 4):3, (4, 5):3, (4, 6):3, (6, 7):3}

    if argv.E_7:
        G = Weyl.build_E7()
        assert G.matrix() == {
            (0, 6): 3, (1, 2): 3, (2, 3): 3, (3, 4): 3, (3, 5): 3, (5, 6): 3}

    if argv.E_6:
        G = Weyl.build_E6()
        assert G.matrix() == {
            (0, 1): 3, (1, 2): 3, (2, 3): 3, (2, 4): 3, (4, 5): 3}
        #items = mulclose(G.gen, verbose=True)
        #print "|E6|=", len(items)

    if argv.F_4:
        G = Weyl.build_F4()
        assert G.matrix() == {(0, 1): 3, (1, 2): 4, (2, 3): 3}
        #items = mulclose(G.gen, verbose=True)
        #print "|F4|=", len(items) # == 1152

    if argv.G_2:
        G = Weyl.build_G2()
        assert G.matrix() == {(0, 1): 6}
        assert len(mulclose(G.gen)) == 12 # dihedral group D_6

    for arg in argv:
        if len(arg) < 3 or arg[1]!="_":
            continue
        try:
            n = int(arg[2:])
        except:
            continue

        print "constructing %s"%arg
        if arg.startswith("A"):
            G = Weyl.build_A(n)
        if arg.startswith("B"):
            G = Weyl.build_B(n)
        if arg.startswith("C"):
            G = Weyl.build_C(n)
        if arg.startswith("D"):
            G = Weyl.build_D(n)

        print "roots:", len(G.roots)
        if argv.order:
            print "order:", len(mulclose(G.gen))

    if G is None:
        return

    if argv.longest_element:
        g = G.longest_element()
        print g

    if argv.monoid:
        test_monoid(G)

    if argv.representation:
        representation(G)

    if argv.find_negi:
        find_negi(G)


if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:

        main()






