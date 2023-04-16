#!/usr/bin/env python3

"""

An attempt to implement the ideas in [GH12] & [Hei21].

[GH12] 
"The Weil representation in characteristic two"
Shamgar Gurevich, Ronny Hadani
https://www.maths.ed.ac.uk/~v1ranick/papers/gurevichhadani2.pdf

[Hei21]
"On stabiliser techniques and their _application
to simulation and certification of quantum devices"
Markus Heinrich
https://kups.ub.uni-koeln.de/50465/1/dissertation_heinrich.pdf

"""


from time import time
start_time = time()

from random import shuffle, seed, choice, randint
#seed(0) # ?
from functools import reduce
from operator import mul

import numpy
from numpy import dot, alltrue, zeros, array, identity, kron, concatenate

from bruhat.action import mulclose, mulclose_hom, Perm
from bruhat.util import cross
from bruhat.argv import argv


int_scalar = numpy.int8

def zeros2(m, n):
    return numpy.zeros((m, n), dtype=int_scalar)


_sy_cache = {}
def symplectic_form(n):
    if n in _sy_cache:
        return _sy_cache[n]
    F = numpy.zeros((2*n, 2*n), dtype=int_scalar)
    I = numpy.identity(n, dtype=int_scalar)
    F[:n, n:] = I
    F[n:, :n] = I
    _sy_cache[n] = F
    return F


class Vector(object):

    def __init__(self, v):
        n = len(v)
        assert n%2 == 0
        n //= 2
        self.n = n
        self.v = array(v, dtype=int_scalar) % 2

    def __str__(self):
        return str(self.v)
    __repr__ = __str__

    def __lt__(u, v):
        return u.v.tobytes() < v.v.tobytes()

    def __eq__(u, v):
        assert u.n == v.n
        return alltrue(u.v == v.v)

    def __hash__(u):
        key = u.v.tobytes()
        return hash(key)

    @classmethod
    def promote(cls, v):
        if isinstance(v, cls):
            return v
        return cls(v)

    def sy(u, v):
        assert u.n == v.n
        n = u.n
        F = symplectic_form(n)
        a = dot(u.v, dot(F, v.v)) % 2
        return a

    @property
    def X(u):
        return u.v[:u.n]

    @property
    def Z(u):
        return u.v[u.n:]

#    def beta(u, v):
#        assert u.n == v.n
#        n = u.n
#        d = dot(u.Z, v.X) % 2
#        return 2*d

    def beta(u, v):
        # [Hei21] 3.70 (this is different to beta in [GH12])
        rhs = (gamma(u+v) - gamma(u) - gamma(v) + 2*(dot(u.Z, v.X))) % 4
        return rhs

    def gamma(u):
        # [Hei21] (3.69)
        d = dot(u.Z, u.X) % 4
        return d

    def __add__(u, v):
        assert u.n == v.n
        w = (u.v + v.v) % 2
        return Vector(w)

    def __matmul__(u, v):
        w = concatenate((u.X, v.X, u.Z, v.Z))
        return Vector(w)


beta = lambda u, v : u.beta(v)
sy = lambda u, v : u.sy(v)
gamma = lambda u : u.gamma()

_vec_cache = {}
def vecspace(n):
    if n in _vec_cache:
        return _vec_cache[n]
    items = [(0, 1)]*n
    space = []
    for v in cross(items):
        space.append(Vector(v))
    _vec_cache[n] = space
    return space


class Heisenberg(object):

    def __init__(self, v, phase=0):
        assert 0<=phase<4 # Z/4
        v = Vector.promote(v)
        self.v = v
        self.n = v.n
        self.phase = phase

    def __str__(self):
        return "%s[%s]"%(self.phase, self.v)
    __repr__ = __str__

    @property
    def X(u):
        return u.v.X

    @property
    def Z(u):
        return u.v.Z

    def __eq__(u, v):
        assert u.n == v.n
        return u.phase == v.phase and u.v == v.v

    def __hash__(u):
        key = (u.v.v.tobytes(), u.phase)
        return hash(key)

    def __neg__(u):
        return Heisenberg(u.v, (u.phase+2)%4)
        
    def __mul__(u, v):
        # [Hei21] (3.64)
        assert u.n == v.n
        phase = (beta(u.v, v.v) + u.phase + v.phase) % 4
        w = u.v + v.v
        return Heisenberg(w, phase)

    def __matmul__(u, v):
        phase = (u.phase + v.phase) % 4
        v = concatenate((u.X, v.X, u.Z, v.Z))
        return Heisenberg(v, phase)

    def __pow__(self, count):
        assert count >= 0
        if count == 0:
            return Heisenberg(Vector([0]*self.n*2), 0)
        if count == 1:
            return self
        A = self
        while count > 1:
            A = self * A
            count -= 1
        return A


class Symplectic(object):
    def __init__(self, A):
        A = array(A, dtype=int_scalar)
        m, n = A.shape
        assert m==n
        assert n%2 == 0
        n //= 2
        self.A = A
        self.n = n
        #self.check()

    def check(A):
        n, A = A.n, A.A
        F = symplectic_form(n)
        B = dot(dot(A.transpose(), F), A) % 2
        assert alltrue(B==F)

    def __str__(self):
        return str(self.A)

    def __eq__(A, B):
        assert A.n == B.n
        return numpy.alltrue(A.A==B.A)

    def __hash__(A):
        return hash(A.A.tobytes())

    @classmethod
    def identity(cls, n):
        A = identity(2*n, dtype=int_scalar)
        return Symplectic(A)

    @classmethod
    def h_gate(cls, n, idx=0):
        A = zeros2(2*n, 2*n)
        for i in range(2*n):
            if i==idx:
                A[i, n+i] = 1
            elif i==n+idx:
                A[i, i-n] = 1
            else:
                A[i, i] = 1
        return Symplectic(A)

    @classmethod
    def cx_gate(cls, n, src=0, tgt=1):
        assert src!=tgt
        A = zeros2(2*n, 2*n)
        for i in range(2*n):
            A[i, i] = 1
        A[src, tgt] = 1
        A[tgt+n, src+n] = 1
        return Symplectic(A)

    @classmethod
    def cz_gate(cls, n, src=0, tgt=1):
        CN = cls.cx_gate(n, src, tgt)
        H = cls.h_gate(n, tgt)
        CZ = H * CN * H
        return CZ

    @classmethod
    def s_gate(cls, n, i=0):
        A = cls.identity(n).A
        assert 0<=i<n
        A[i, i+n] = 1
        return Symplectic(A)

    def inverse(self):
        A = self.A
        n = self.n
        F = symplectic_form(n)
        Fi = F # an involution
        Ai = dot(Fi, dot(A.transpose()), F) % 2
        return Symplectic(Ai)

    def __call__(A, v):
        assert isinstance(v, Vector), type(v)
        u = dot(A.A, v.v) % 2
        return Vector(u)

    def __mul__(A, B):
        assert isinstance(B, Symplectic)
        C = dot(A.A, B.A) % 2
        return Symplectic(C)

    def __pow__(self, count):
        if count < 0:
            return self.inverse()**(-count) # recurse
        if count == 0:
            return Symplectic.identity(self.n)
        if count == 1:
            return self
        A = self
        while count > 1:
            A = self * A
            count -= 1
        return A


class ASymplectic(object):
    def __init__(self, g, alpha):
        assert isinstance(g, Symplectic)
        self.n = g.n
        self.g = g
        self.alpha = alpha
        assert len(alpha) == 2**(2*g.n)
        #self.check()

    def __hash__(g):
        n, g, alpha = g.n, g.g, g.alpha
        keys = list(alpha.keys())
        keys.sort()
        vals = tuple(alpha[k] for k in keys)
        return hash((g.A.tobytes(), vals))

    def __eq__(g, h):
        assert isinstance(h, ASymplectic)
        assert g.n == h.n
        return g.g == h.g and g.alpha == h.alpha

    def __str__(self):
        return "%s %s"%(self.g, self.alpha)
        
    def check(self):
        n, g, alpha = self.n, self.g, self.alpha
        space = vecspace(2*n)
        for v in space:
          for w in space:
            lhs = (beta(g(v), g(w)) - beta(v, w)) % 4
            rhs = (alpha[v+w] - alpha[v] - alpha[w]) % 4
            assert lhs == rhs
            #assert alpha[v+w] == alpha[v] + alpha[w]

    def __call__(g, v):
        assert isinstance(v, Heisenberg)
        assert g.n == v.n
        w = g.g(v.v)
        phase = (v.phase + g.alpha[v.v]) % 4
        return Heisenberg(w, phase)

    def __mul__(g, h):
        # eq. 3.73 in [Hei21], sec 1.3.1 in [GH12]
        assert isinstance(h, ASymplectic)
        assert g.n == h.n
        gh = g.g*h.g
        alpha = {}
        for v in vecspace(2*g.n):
            alpha[v] = (g.alpha[h.g(v)] + h.alpha[v]) % 4
        return ASymplectic(gh, alpha)



def find_alpha(g):
    # alpha is determined by its values on the basis elements
    # (no it's not linear!) but here we search through all functions
    # (dict's) space -> 2Z/4Z . There are always N solutions .
    n = g.n
    space = vecspace(2*n)
    lhs = {}
    pairs = [(v, w) for v in space for w in space]
    for (v, w) in pairs:
        lhs[v, w] = (beta(g(v), g(w)) - beta(v, w)) % 4

    N = 2**(2*n)
    for value in cross([(0,1)]*N):
        alpha = {}
        for idx, v in enumerate(space):
            alpha[v] = 2*value[idx]
        #print(alpha)

        for (v, w) in pairs:
            #if alpha[v] + alpha[w] != alpha[v+w]:
            #    break
            #lhs = (beta(g(v), g(w)) - beta(v, w)) % 4
            rhs = (alpha[v+w] - alpha[v] - alpha[w]) % 4
            if lhs[v, w] != rhs:
                break
        else:
            yield alpha



def get_order(g):
    n = 1
    a = g
    while a*a != a: # identity
        a = a*g
        n += 1
    return n

def get_identity(G):
    I = [g for g in G if g*g==g]
    assert len(I) == 1
    return I[0]


def find_pairgen(gens, G=None):
    if G is None:
        G = mulclose(gens)
    print("find_pairgen(%d, %d)"%(len(gens), len(G)))
    N = len(gens)
    if N>2:
        pairs = [(gens[i], gens[j]) for i in range(N) for j in range(i+1, N)]
        for pair in pairs:
            H = mulclose(pair)
            if len(H) == len(G):
                print("found pair")
                return pair # <----------- return
            print(len(H), end=".", flush=True)
    return gens

def get_perms(gens, G):
    items = list(range(len(G)))
    lookup = dict((g, i) for (i, g) in enumerate(G))
    perms = []
    for g in gens:
        i = lookup[g]
        perm = [lookup[g*h] for h in G]
        perm = dict(enumerate(perm))
        perm = Perm(perm, items)
        perms.append(perm)
    return perms


def main():

    space = vecspace(4)

    # test beta is a 2-cocycle
    for g in space:
      for h in space:
        for k in space:
            lhs = (beta(g, h) - beta(h, k))%4
            rhs = (beta(g, h+k) - beta(g+h, k))%4
            assert lhs == rhs # [Hei21] eq. 3.65

    for u in space:
      for v in space:
        assert beta(u, v) % 2 == sy(u, v)
        assert (beta(u, v) - beta(v, u))%4 == 2*sy(u, v) # [Hei21] eq. 3.67
        rhs = (gamma(u+v) - gamma(u) - gamma(v) + 2*(dot(u.Z, v.X))) % 4
        assert beta(u, v) == rhs

    # ------------------------------------------------

    # wtf is going on here, i don't know...

    I = Heisenberg([0, 0])
    w = Heisenberg([0, 0], 1)
    X = Heisenberg([1, 0])
    Z = Heisenberg([0, 1])
    Y = Heisenberg([1, 1], 2)  # ?
    XZ = Heisenberg([1, 1], 1) # ?
    ZX = Heisenberg([1, 1], 3) # ?

    assert X*X == I
    assert Z*Z == I
    assert XZ != ZX
    assert XZ == -ZX
    assert X*Z == XZ
    assert Z*X == ZX
    assert Y == w*X*Z

    assert w != I
    assert w**2 != I
    assert w**3 != I
    assert w**4 == I

    XI = X@I
    ZI = Z@I
    IX = I@X
    IZ = I@Z
    XZ = X@Z
    ZX = Z@X
    wII = w@I

#    assert XI*IX == X@X
#    assert ZI*XI == -XI*ZI
#    assert XZ * ZX == ZX * XZ

    G = mulclose([w, X, Y, Z])
    assert w in G
    assert len(G) == 16, len(G)

    gens = [wII, XI, ZI, IX, IZ]
    HW = mulclose(gens)
    assert len(HW) == 64

    #rels = get_presentation(gens, HW)
    #return

    # ------------------------------------------------

    n = 1
    S = Symplectic.s_gate(n)
    H = Symplectic.h_gate(n)

    Sp2 = mulclose([S, H])
    assert len(Sp2) == 6

    # ------------------------------------------------

    n = 2
    SI = Symplectic.s_gate(n)
    HI = Symplectic.h_gate(n)
    IS = Symplectic.s_gate(n, 1)
    IH = Symplectic.h_gate(n, 1)
    CZ = Symplectic.cz_gate(n)

    Sp4 = mulclose([SI, HI, IS, IH, CZ])
    assert len(Sp4) == 720, len(Sp4)

    # ------------------------------------------------

    space = vecspace(2*n)
    assert len(space) == 2**(2*n)

    Sp = list(Sp4)
    shuffle(Sp)

    tgt = []
    gens = []
    for g in Sp:
        #print(g)
        for alpha in find_alpha(g):
            #print(alpha)
            op = ASymplectic(g, alpha)
            gens.append(op)
            tgt.append(g)
            break
        if len(gens) > 5:
            break

    for g in gens:
      for h in gens:
        for a in HW:
            assert (g*h)(a) == g(h(a))

    hom = mulclose_hom(gens, tgt)
    ASp = list(hom.keys())
    shuffle(ASp)
    for _ in range(64):
        g, h = choice(ASp), choice(ASp)
        assert hom[g]*hom[h] == hom[g*h]
    assert len(ASp) == 11520, "whoops, try again?"
    print("|ASp| =", len(ASp))

    I0 = get_identity(Sp)
    I = get_identity(ASp)
    section = {t:s for (s,t) in hom.items()}
    section[I0] = I
    is_identity = lambda g : g*g==g
    A = []
    for g in ASp:
        if g.g == I0:
            A.append(g)
    def get_inverse(G, g):
        for h in G:
            if g*h == I:
                return h
        assert 0
                
    print("%s >--> %s -->> %s"%(len(A), len(ASp), len(Sp)))

    # See: https://math.stanford.edu/~conrad/210BPage/handouts/Cohomology&Extensions.pdf

    decompose = {} # ASp --> A x Sp
    for a in A:
      for x in Sp:
        g = a*section[x]
        assert g not in decompose
        decompose[g] = (a, x)
    print("decompose:", len(decompose))

    f = lambda x : section[x]
    def act(x, b):
        fx = f(x)
        fxi = get_inverse(ASp, fx)
        return fx*b*fxi

    for _ in range(10):
        g = choice(ASp)
        h = choice(ASp)
        a, x = decompose[g]
        b, y = decompose[h]
        fx = f(x)
        #xi = get_inverse(Sp, x)
        lhs = g*h
        rhs = a * act(x, b)*(f(x) * f(y) * get_inverse(ASp, f(x*y)))*f(x*y)
        assert lhs==rhs
        print(lhs)
        print("ok")

    return

    #for g in Sp:
    for _ in range(10):
        g = choice(Sp)
        x = section[g]
        xi = get_inverse(ASp, x)
        for b in A:
            assert x*b*xi in A

    if 0:
        # Success! ASp is isomorphic to aff_sy constructed in dev/clifford.gap
        from bruhat.oeqc import gap_code
        gens = find_pairgen(gens, ASp)
        perms = get_perms(gens, ASp)
        print(str(perms)[:200])
        print(gap_code(perms)[:200])
    
        f = open("gap_ASp.gap", 'w')
        print(gap_code(perms), file=f)
        f.close()

    # ------------------------------------------------

    if 0:
        found = []
        for g in Sp:
            total = 0
            for v in space:
              for w in space:
                lhs = (beta(g(v), g(w)) - beta(v, w)) % 4
                total += lhs
            #print(total, end=" ")
            if total == 0:
                found.append(g)
        assert len(mulclose(found)) == 6 # hmmm
    


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

    print("OK: finished in %.3f _seconds"%(time() - start_time))
    print()

