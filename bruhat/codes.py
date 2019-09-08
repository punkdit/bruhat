#!/usr/bin/env python3

"""
Classical binary linear codes.
Build Reed-Muller codes.
Construct triorthogonal matrices.
MacWilliams identities.
"""

from random import random, randint
from functools import reduce
from operator import mul

import numpy

from bruhat.solve import array2, zeros2, dot2, shortstr, rank, find_kernel, span
from bruhat.solve import linear_independent, parse, pseudo_inverse, eq2, rand2
from bruhat.action import mulclose
from bruhat.comm import Poly
from bruhat.argv import argv
from bruhat.util import choose, cross


class Code(object):
    """
        binary linear code, as defined by a generator matrix.
    """
    def __init__(self, G, H=None, d=None, desc="", check=True):
        assert len(G.shape)==2
        self.G = G.copy()
        self.k, self.n = G.shape
        self.d = d
        self.desc = desc

        if H is None:
            H = list(find_kernel(G))
            H = array2(H)
        if H.shape == (0,):
            H.shape = (0, self.n)
        self.H = H.copy()

        if check:
            self.check()

    def check(self):
        G, H = self.G, self.H
        assert rank(G) == len(G)
        A = dot2(H, G.transpose())
        assert A.sum()==0

    def __str__(self):
        desc = ', "%s"'%self.desc if self.desc else ""
        return "Code([[%s, %s, %s]]%s)" % (self.n, self.k, self.d, desc)

    def dump(self):
        G, H = self.G, self.H
        print("G =")
        print(shortstr(G))
    
        print("H =")
        print(shortstr(H))

    def is_selfdual(self):
        #G = self.G
        #x = dot2(G, G.transpose())
        #return x.sum()==0
        return self.eq(self.get_dual())

    def get_dual(self):
        return Code(self.H, self.G)

    def eq(self, other):
        "Two codes are equal if their generating matrices have the same span."
        G1, G2 = self.G, other.G
        if len(G1) != len(G2):
            return False
        A = dot2(self.H, other.G.transpose())
        B = dot2(other.H, self.G.transpose())
        assert (A.sum()==0) == (B.sum()==0)
        return A.sum() == 0

    def get_distance(self):
        G = self.G
        d = None
        for v in span(G):
            w = v.sum()
            if w==0:
                continue
            if d is None or w<d:
                d = w
        if self.d is None:
            self.d = d
        return d

    def puncture(self, i):
        assert 0<=i<self.n
        G = self.G
        A = G[:, :i]
        B = G[:, i+1:]
        G = numpy.concatenate((A, B), axis=1)
        G = linear_independent(G, check=True)
        return Code(G)

    def get_even(self):
        G = self.G
        rows = [row for row in G if row.sum()%2==0]
        G = array2(rows)
        return Code(G)

    def is_morthogonal(self, m):
        return is_morthogonal(self.G, m)

    def is_triorthogonal(self):
        return is_morthogonal(self.G, 3)

    def weight_enum(self, A, B):
        G = self.G
        the_op = None
        for v in span(G):
            op = A if v[0]==0 else B
            for vi in v[1:]:
                op = op * (A if vi==0 else B)
            the_op = op if the_op is None else the_op + op
        return the_op

    def tensor_enum(self, A, B):
        G = self.G
        the_op = None
        for v in span(G):
            op = A if v[0]==0 else B
            for vi in v[1:]:
                op = op @ (A if vi==0 else B)
            the_op = op if the_op is None else the_op + op
        return the_op



def is_morthogonal(G, m):
    k = len(G)
    if m==1:
        for v in G:
            if v.sum()%2 != 0:
                return False
        return True
    if m>2 and not is_morthogonal(G, m-1):
        return False
    items = list(range(k))
    for idxs in choose(items, m):
        v = G[idxs[0]]
        for idx in idxs[1:]:
            v = v * G[idx]
        if v.sum()%2 != 0:
            return False
    return True


def strong_morthogonal(G, m):
    k = len(G)
    assert m>=1
    if m==1:
        for v in G:
            if v.sum()%2 != 0:
                return False
        return True
    if not strong_morthogonal(G, m-1):
        return False
    items = list(range(k))
    for idxs in choose(items, m):
        v = G[idxs[0]]
        for idx in idxs[1:]:
            v = v * G[idx]
        if v.sum()%2 != 0:
            return False
    return True


def reed_muller(r, m, puncture=False):
    "Build Reed-Muller code"

    assert 0<=r<=m, "r=%s, m=%d"%(r, m)

    n = 2**m # length

    one = array2([1]*n)
    basis = [one]

    vs = [[] for i in range(m)]
    for i in range(2**m):
        for j in range(m):
            vs[j].append(i%2)
            i >>= 1
        assert i==0

    vs = [array2(v) for v in vs]

    for k in range(r):
        for items in choose(vs, k+1):
            v = one
            #print(items)
            for u in items:
                v = v*u
            basis.append(v)
        
    G = numpy.array(basis)

    code = Code(G, d=2**(m-r), desc="reed_muller(%d, %d)"%(r, m))

    if puncture:
        code = code.puncture(0)

    return code




class Tensor(object):

    """ Some kind of graded ring element... I*I*I + X*X*X etc.
    """

    zero = 0
    one = 1
    def __init__(self, items, grade=None):
        # map key -> coeff, key is ("A", "B") etc.
        assert items or (grade is not None)
        keys = list(items.keys())
        keys.sort()
        self.items = {}
        nz = []
        for key in keys:
            assert grade is None or grade==len(key)
            grade = len(key)
            v = int(items[key])
            if v != 0:
                self.items[key] = v # uniquify
                nz.append(key)
        self.keys = nz
        self.grade = grade

    def __add__(self, other):
        assert self.grade == other.grade
        items = dict(self.items)
        for (k, v) in other.items.items():
            items[k] = items.get(k, self.zero) + v
        return Tensor(items, self.grade)

    def __sub__(self, other):
        assert self.grade == other.grade
        items = dict(self.items)
        for (k, v) in other.items.items():
            items[k] = items.get(k, self.zero) - v
        return Tensor(items, self.grade)

    def __mul__(self, other):
        items = {}
        for (k1, v1) in self.items.items():
          for (k2, v2) in other.items.items():
            k = k1+k2
            assert k not in items
            items[k] = v1*v2
        return Tensor(items, self.grade+other.grade)
    tensor = __mul__

    def __rmul__(self, r):
        items = {}
        for (k, v) in self.items.items():
            items[k] = r*v
        return Tensor(items, self.grade)

    def subs(self, rename):
        the_op = Tensor({}, self.grade) # zero
        one = self.one
        for (k, v) in self.items.items():
            final = None
            for ki in k:
                op = rename.get(ki, Tensor({ki : one}))
                if final is None:
                    final = op
                else:
                    final = final * op # tensor
            the_op = the_op + v*final
        return the_op

    def __str__(self):
        ss = []
        for k in self.keys:
            v = self.items[k]
            s = ''.join(str(ki) for ki in k)
            if v == 1:
                pass
            elif v == -1:
                s = "-"+s
            else:
                s = str(v)+"*"+s
            ss.append(s)
        ss = '+'.join(ss) or "0"
        ss = ss.replace("+-", "-")
        return ss

    def __repr__(self):
        return "Tensor(%s)"%(self.items)

    def __eq__(self, other):
        return self.items == other.items

    def __ne__(self, other):
        return self.items != other.items

    def __hash__(self):
        return hash((str(self), self.grade))



def test():

    for m in range(2, 7):
      for r in range(0, m+1):
        code = reed_muller(r, m)
        assert code.n == 2**m
        k = 1
        for i in range(1, r+1):
            k += len(list(choose(list(range(m)), i)))
        assert code.k == k
        if code.k < 12:
            assert code.get_distance() == 2**(m-r)

        if 0<=r<=m-1:
            dual = code.get_dual()
            code1 = reed_muller(m-r-1, m)
            #print(m-r-1 == r, dual.eq(code), code)
            assert code1.eq(dual)

    test_hamming()

    print("OK")


def genus_enum1():
    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)
    puncture = argv.puncture
    code = reed_muller(r, m, puncture)
    G = code.G
    print(shortstr(G))
    m, n = G.shape

    poly = lambda cs : Poly(cs, 2, "x_0 x_1".split())
    #x_1 = poly({(0, 1) : 1})
    #x_0 = poly({(1, 0) : 1})
    #xs = [x0, x1]
    cs = {}
    for v0 in span(G):
        exp = [0, 0]
        #vv = numpy.array([2*v0, v1])
        for i in range(n):
            exp[v0[i]] += 1
        exp = tuple(exp)
        cs[exp] = cs.get(exp, 0) + 1
    p = poly(cs)
    print(p)
    
    Z = numpy.array([[1, 0], [0, 1]])
    q = p.transform(Z)
    print("invariant under Z", q==p)

    S2 = numpy.array([[0, 1], [-1, 0]])
    q = p.transform(S2)
    print("invariant under CS2", q==p)



def genus_enum2():
    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)
    puncture = argv.puncture
    code = reed_muller(r, m, puncture)
    G = code.G
    print(shortstr(G))
    m, n = G.shape

    poly = lambda cs : Poly(cs, 4, "x_{00} x_{01} x_{10} x_{11}".split())
#    x_11 = poly({(0, 0, 0, 1) : 1})
#    x_10 = poly({(0, 0, 1, 0) : 1})
#    x_01 = poly({(0, 1, 0, 0) : 1})
#    x_00 = poly({(1, 0, 0, 0) : 1})
#    xs = [[x_00, x_01], [x_10, x_11]]
    cs = {}
    for v0 in span(G):
        print(".",end='',flush=True)
        for v1 in span(G):
            exp = [0, 0, 0, 0]
            #vv = numpy.array([2*v0, v1])
            for i in range(n):
                exp[2*v0[i] + v1[i]] += 1
            exp = tuple(exp)
            cs[exp] = cs.get(exp, 0) + 1
        #break
    print()
    p = poly(cs)
    print(p)
    
    CZ = numpy.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, -1]])
    q = p.transform(CZ)
    print("invariant under CZ", q==p)

    CS2 = numpy.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, -1, 0]])
    q = p.transform(CS2)
    print("invariant under CS2", q==p)

    T2 = numpy.array([
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [-1, 0, 0, 0]])
    q = p.transform(T2)
    print("invariant under T2", q==p)

    A = numpy.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 1, 0, 0],
        [-1, 0, 0, 0]])
    q = p.transform(A)
    print("invariant under A", q==p)



def genus_enum3():
    r = argv.get("r", 1) # degree
    m = argv.get("m", 4)
    puncture = argv.puncture
    code = reed_muller(r, m, puncture)
    G = code.G
    print(shortstr(G))
    m, n = G.shape

    rank = 8
    poly = lambda cs : Poly(cs, rank, 
        "x_{000} x_{001} x_{010} x_{011} x_{100} x_{101} x_{110} x_{111}".split())
    cs = {}
    for v0 in span(G):
        print(".",end='',flush=True)
        for v1 in span(G):
          for v2 in span(G):
            exp = [0, 0, 0, 0, 0, 0, 0, 0]
            #vv = numpy.array([2*v0, v1])
            for i in range(n):
                exp[4*v0[i] + 2*v1[i] + v2[i]] += 1
            exp = tuple(exp)
            cs[exp] = cs.get(exp, 0) + 1
        #break
    print()
    p = poly(cs)
    print(p)
    print()
    
    CCZ = numpy.array([
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, -1]])
    q = p.transform(CCZ)
    print("invariant under CCZ", q==p)

    CS2 = numpy.array([
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, -1, 0]])
    q = p.transform(CS2)
    print("invariant under CS2", q==p)

    T3 = numpy.array([
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, -1, 0, 0, 0]])
    q = p.transform(T3)
    print("invariant under T3", q==p)

    A = numpy.array([
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, -1, 0, 0, 0]])
    q = p.transform(A)
    print("invariant under A", q==p)




def test_tri_rm():
    r = argv.get("r", 1) # degree
    m = argv.get("m", 4)
    puncture = argv.puncture
    code = reed_muller(r, m, puncture)
    G = code.G
    rows = [row for row in G if row.sum()%2==0]
    G = array2(rows)
    k = len(rows)
    print(G.shape)
    for i in range(1, 10):
        # statistical test for m-orthogonality
        for trials in range(10000):
            vecs = [rows[randint(0, k-1)] for j in range(i)]
            v = vecs[0].copy()
            for w in vecs[1:]:
                v = v*w
            if v.sum()%2 != 0:
                break
        else:
            print("(%d)"%i, end=" ", flush=True)
            continue
        break
    print()


def test_rm():
    params = [(r, m) for m in range(2, 8) for r in range(1, m)]
    r = argv.get("r", None) # degree
    m = argv.get("m", None)
    if r is not None and m is not None:
        params = [(r, m)]
    
    for (r, m) in params:
        #code = reed_muller(r, m)
#      for code in [ reed_muller(r, m), reed_muller(r, m).puncture(0) ]:
      for code in [reed_muller(r, m)]:
        if argv.puncture:
            print(code, end=" ", flush=True)
            code = code.puncture(0)
            code = code.get_even()
            if argv.puncture==2:
                code = code.puncture(0)
                code = code.get_even()
            G = code.G
            k, n = G.shape
            #code = Code(G)
            #d = code.get_distance()
            d = "."
            print("puncture [%d, %d, %s]" % (n, k, d), end=" ", flush=True)
        else:
            G = code.G
            print(code, end=" ", flush=True)
        i = 1
        while i<8:
            if (is_morthogonal(G, i)):
                print("(%d)"%i, end="", flush=True)
                i += 1
            else:
                break
            if i > code.k:
                print("*", end="")
                break
        print()
        if argv.show:
            print(G.shape)
            print(shortstr(G))
            print(dot2(G, G.transpose()).sum())
        if len(G) >= 14:
            continue
        A = array2(list(span(G)))
        for i in [1, 2, 3]:
            assert strong_morthogonal(G, i) == strong_morthogonal(A, i)

def gen():
    r = argv.get("r", None) # degree
    m = argv.get("m", None)

    if r is not None and m is not None:
        code = reed_muller(r, m)
    
        #print(code)
        #print("d =", code.get_distance())
        #code.dump()
    
        #code = code.puncture(3)
    
        #print(code)
        code = code.puncture(0)
        print(code)
        for g in code.G:
            print(shortstr(g), g.sum())
        print()
        #code.dump()
        #print("d =", code.get_distance())
    
        return

    for m in range(2, 8):
      for r in range(0, m+1):
        code = reed_muller(r, m)
        print(code, end=" ")
        if code.is_selfdual():
            print("is_selfdual", end=" ")
        if code.is_morthogonal(2):
            print("is_biorthogonal", end=" ")
        if code.is_morthogonal(3):
            print("is_triorthogonal", end=" ")
        if dot2(code.H, code.H.transpose()).sum()==0:
            print("***", end=" ")
        p = code.puncture(0)
        if p.is_morthogonal(3):
            print("puncture.is_triorthogonal", end=" ")
        if p.is_selfdual():
            print("puncture.is_selfdual", end=" ")
        if dot2(p.H, p.H.transpose()).sum()==0:
            print("***", end=" ")
        print()

        if p.is_triorthogonal() and p.k < 20:
            G = p.G
            #print(shortstr(G))
            A = list(span(G))
            A = array2(A)
            print(is_morthogonal(A, 3))


def test_triorth():
    code = reed_muller(1, 5)
    code = code.puncture(0)
    code.dump()

    print(code.is_triorthogonal())
    A = array2(list(span(code.G)))
    print(is_morthogonal(A, 2))
    #print(shortstr(A))

    k = len(A)

    for i in range(k):
      for j in range(i+1, k):
        u = A[i]
        v = A[j]
        x = (u*v).sum() % 2
        if x == 0:
            continue
        #print(shortstr(u))
        #print(shortstr(v))
        #print()

    for a in range(k):
      for b in range(a+1, k):
       for c in range(b+1, k):
        u = A[a]
        v = A[b]
        w = A[c]
        x = (u*v*w).sum() % 2
        #if x:
            #print(a, b, c)


def test_hamming():
    G = parse("""
    1....111
    .1..1.11
    ..1.11.1
    ...1111.
    """)

    assert strong_morthogonal(G, 2)


def main():
    I = Tensor({"I" : 1})
    X = Tensor({"X" : 1})
    Y = Tensor({"Y" : 1})
    Z = Tensor({"Z" : 1})

    II = I*I
    XI = X*I
    IX = I*X
    XX = X*X
    assert II+II == 2*II

    assert X*(XI + IX) == X*X*I + X*I*X

    assert ((I-Y)*I + I*(I-Y)) == 2*I*I - I*Y - Y*I
    assert (XI + IX).subs({"X": I-Y}) == ((I-Y)*I + I*(I-Y))

    A = Tensor({"A":1})
    B = Tensor({"B":1})
    p = A*A*A + B*B*A + B*A*B + A*B*B
    q = A*A*A + B*B*B
    p1 = p.subs({"A": A+B, "B": A-B})
    assert p1 == 4*A*A*A + 4*B*B*B

    verbose = argv.verbose

    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)

    code = reed_muller(r, m)
    code.dump()

    p = code.weight_enum(A, B)
    if verbose:
        print("code:")
        print(p)

    dual = code.get_dual()
    q = dual.weight_enum(A, B)
    if verbose:
        print("dual:")
        print(q)
    print("p==q:", p==q)
    print("code.is_selfdual:", code.is_selfdual())

    #r = p.subs({"A": A+B, "B": A-B})
    r = code.weight_enum(A+B, A-B)
    if verbose:
        print("P(A+B, A-B)")
        print(r)
    coeff = 2**len(code.G)
    print("MacWilliams:", r == coeff*q)


    print("OK")


def test_dual():

    from vec import Space, Hom, Map

    import element

    ring = element.Z
    #ring = element.Q
    one = ring.one

    space = Space(2, ring)
    hom = Hom(space, space)

    I = Map.from_array([[1, 0], [0, 1]], hom)
    X = Map.from_array([[0, 1], [1, 0]], hom)
    Z = Map.from_array([[1, 0], [0, -1]], hom)

    assert X+X == (2*X)
    assert X*Z == -Z*X

    assert (X@X) * (Z@Z) == (Z@Z) * (X@X)
    assert (I@X@X) * (I@Z@Z) == (I@Z@Z) * (I@X@X)

    if argv.code == "repitition":
        G = parse("111")
        H = parse("""
        11.
        .11""")
    elif argv.code == "steane":
        G = parse("""
        1111...
        .1.11.1
        ..11.11
        """)
    else:
        return

    code = Code(G)
    if argv.dual:
        code = code.get_dual()

    dual = code.get_dual()

    code.dump()

    #W = lambda A, B : (A@A@A + B@B@B)
    #WD = lambda A, B : (A@A@A + B@B@A + B@A@B + A@B@B)
    W = code.tensor_enum
    WD = dual.tensor_enum

    a = 2**len(code.G)
    b = 2**len(dual.G)
    A = WD(I, X)
    B = W(I, Z)
    assert A*B == B*A
    AA = W(I+X, I-X)
    BB = WD(I+Z, I-Z)
    assert AA == a*A
    assert BB == b*B
    assert AA*BB == a*b*A*B
    #print(W(I+X, I-X))
    #print(WD(I, X))
    #print(W(I, Z))

    src = Space(2, ring)
    tgt = Space(3, ring)
    hom = Hom(src, tgt)

    A = Map.from_array([
        [1, 2],
        [5, 6],
        [7, 8]], hom)

    B = Map.from_array([
        [7, 2],
        [6, 3],
        [5, 4]], hom)

    assert W(A+B, A-B) == a*WD(A, B)
    assert WD(A+B, A-B) == b*W(A, B)




def search():
    # Bravyi, Haah, 1209.2426v1 sec IX.
    # https://arxiv.org/pdf/1209.2426.pdf

    verbose = argv.get("verbose")
    m = argv.get("m", 6) # _number of rows
    k = argv.get("k", None) # _number of odd-weight rows

    # these are the variables N_x
    xs = list(cross([(0, 1)]*m))

    maxweight = argv.maxweight
    minweight = argv.get("minweight", 1)

    xs = [x for x in xs if minweight <= sum(x)]
    if maxweight:
        xs = [x for x in xs if sum(x) <= maxweight]

    N = len(xs)

    lhs = []
    rhs = []

    # bi-orthogonality
    for a in range(m):
      for b in range(a+1, m):
        v = zeros2(N)
        for i, x in enumerate(xs):
            if x[a] == x[b] == 1:
                v[i] = 1
        if v.sum():
            lhs.append(v)
            rhs.append(0)

    # tri-orthogonality
    for a in range(m):
      for b in range(a+1, m):
       for c in range(b+1, m):
        v = zeros2(N)
        for i, x in enumerate(xs):
            if x[a] == x[b] == x[c] == 1:
                v[i] = 1
        if v.sum():
            lhs.append(v)
            rhs.append(0)

#    # dissallow columns with weight <= 1
#    for i, x in enumerate(xs):
#        if sum(x)<=1:
#            v = zeros2(N)
#            v[i] = 1
#            lhs.append(v)
#            rhs.append(0)

    if k is not None:
      # constrain to k _number of odd-weight rows
      assert 0<=k<m
      for a in range(m):
        v = zeros2(N)
        for i, x in enumerate(xs):
          if x[a] == 1:
            v[i] = 1
        lhs.append(v)
        if a<k:
            rhs.append(1)
        else:
            rhs.append(0)

    A = array2(lhs)
    rhs = array2(rhs)
    #print(shortstr(A))

    B = pseudo_inverse(A)
    soln = dot2(B, rhs)
    if not eq2(dot2(A, soln), rhs):
        print("no solution")
        return
    if verbose:
        print("soln:")
        print(shortstr(soln))

    soln.shape = (N, 1)
    rhs.shape = A.shape[0], 1

    K = array2(list(find_kernel(A)))
    #print(K)
    #print( dot2(A, K.transpose()))
    #sols = []
    #for v in span(K):
    best = None
    density = 1.0
    size = 99*N
    trials = argv.get("trials", 1024)
    count = 0
    for trial in range(trials):
        u = rand2(len(K), 1)
        v = dot2(K.transpose(), u)
        #print(v)
        v = (v+soln)%2
        assert eq2(dot2(A, v), rhs)

        if v.sum() > size:
            continue
        size = v.sum()

        Gt = []
        for i, x in enumerate(xs):
            if v[i]:
                Gt.append(x)
        if not Gt:
            continue
        Gt = array2(Gt)
        G = Gt.transpose()
        assert is_morthogonal(G, 3)
        if G.shape[1]<m:
            continue

        if 0 in G.sum(1):
            continue

        if argv.strong_morthogonal and not strong_morthogonal(G, 3):
            continue

        #print(shortstr(G))
#        for g in G:
#            print(shortstr(g), g.sum())
#        print()

        _density = float(G.sum()) / (G.shape[0]*G.shape[1])
        #if best is None or _density < density:
        if best is None or G.shape[1] <= size:
            best = G
            size = G.shape[1]
            density = _density

        if 0:
            #sols.append(G)
            Gx = even_rows(G)
            assert is_morthogonal(Gx, 3)
            if len(Gx)==0:
                continue
            GGx = array2(list(span(Gx)))
            assert is_morthogonal(GGx, 3)

        count += 1

    print("found %d solutions" % count)
    if best is None:
        return

    G = best
    #print(shortstr(G))
    for g in G:
        print(shortstr(g), g.sum())
    print()
    print("density:", density)
    print("shape:", G.shape)

    G = linear_independent(G)
    A = list(span(G))
    print(strong_morthogonal(A, 1))
    print(strong_morthogonal(A, 2))
    print(strong_morthogonal(A, 3))
    
    #print(shortstr(dot2(G, G.transpose())))

    if 0:
        B = pseudo_inverse(A)
        v = dot2(B, rhs)
        print("B:")
        print(shortstr(B))
        print("v:")
        print(shortstr(v))
        assert eq2(dot2(B, v), rhs) 


def build_toric(l=3, allgen=True):
    keys = []
    keymap = {}
    for i in range(l):
      for j in range(l):
        for k in (0, 1):
            m = len(keys)
            keys.append((i, j, k))
            for di in (-l, 0, l):
              for dj in (-l, 0, l):
                keymap[i+di, j+dj, k] = m

    if l>2:
        assert keys[keymap[2, 1, 0]] == (2, 1, 0)

    if allgen:
        m = l**2 # rows (constraints)
    else:
        m = l**2-1 # rows (constraints)
    n = len(keys) # cols (bits)
    assert n == 2*(l**2)

    Lx = zeros2(2, n)
    Lz = zeros2(2, n)
    Hx = zeros2(m, n)
    Tz = zeros2(m, n)
    Hz = zeros2(m, n)
    Tx = zeros2(m, n)

    for i in range(l):
        Lx[0, keymap[i, l-1, 1]] = 1
        Lx[1, keymap[l-1, i, 0]] = 1
        Lz[0, keymap[0, i, 1]] = 1
        Lz[1, keymap[i, 0, 0]] = 1

    row = 0
    xmap = {}
    for i in range(l):
      for j in range(l):
        if (i, j)==(0, 0) and not allgen:
            continue
        Hx[row, keymap[i, j, 0]] = 1
        Hx[row, keymap[i, j, 1]] = 1
        Hx[row, keymap[i-1, j, 0]] = 1
        Hx[row, keymap[i, j-1, 1]] = 1
        xmap[i, j] = row
        i1 = i
        while i1>0:
            Tz[row, keymap[i1-1, j, 0]] = 1
            i1 -= 1
        j1 = j
        while j1>0:
            Tz[row, keymap[i1, j1-1, 1]] = 1
            j1 -= 1
        row += 1

    row = 0
    zmap = {}
    for i in range(l):
      for j in range(l):
        if i==l-1 and j==l-1 and not allgen:
            continue
        Hz[row, keymap[i, j, 0]] = 1
        Hz[row, keymap[i, j, 1]] = 1
        Hz[row, keymap[i+1, j, 1]] = 1
        Hz[row, keymap[i, j+1, 0]] = 1
        zmap[i, j] = row
        i1 = i
        while i1<l-1:
            Tx[row, keymap[i1+1, j, 1]] = 1
            i1 += 1
        j1 = j
        while j1<l-1:
            Tx[row, keymap[i1, j1+1, 0]] = 1
            j1 += 1
        row += 1

    return Hx


def search_extend():
    # Extend the checks of a random code to make it triorthogonal.
    # Based on the search function above.

    verbose = argv.get("verbose")

    m = argv.get("m", 6)
    n = argv.get("n", m+2)
    k = argv.get("k") # odd _numbered rows ( logical operators)
    code = argv.get("code", "rand")

    if code == "rand":
        while 1:
            G0 = rand2(m, n)
            counts = G0.sum(0)
            if min(counts)==2 and rank(G0) == m:
                cols = set()
                for i in range(n):
                    cols.add(tuple(G0[:, i]))
                if len(cols) == n: # no repeated cols
                    break

    elif code == "toric":
        G0 = parse("""
        11.11...
        .111..1.
        1...11.1
        """) # l=2 toric code X logops + X stabs 

        l = argv.get("l", 3)
        G0 = build_toric(l)

        m, n = G0.shape
    else:
        return

    code = Code(G0, check=False)
    print(shortstr(G0))
    print("is_triorthogonal:", code.is_triorthogonal())

    # these are the variables N_x
    xs = list(cross([(0, 1)]*m))
    N = len(xs)

    lookup = {}
    for i, x in enumerate(xs):
        lookup[x] = i

    lhs = []
    rhs = []

    taken = set()
    for i in range(n):
        x = G0[:, i]
        idx = lookup[tuple(x)]
        assert idx not in taken
        taken.add(idx)

    if verbose:
        for idx in range(N):
            print(idx, xs[idx], "*" if idx in taken else "")

    for idx in taken:
        v = zeros2(N)
        v[idx] = 1
        lhs.append(v)
        rhs.append(1)

    # bi-orthogonality
    for a in range(m):
      for b in range(a+1, m):
        v = zeros2(N)
        for i, x in enumerate(xs):
            if x[a] == x[b] == 1:
                v[i] += 1
        assert v.sum()
        lhs.append(v)
        rhs.append(0)

    # tri-orthogonality
    for a in range(m):
      for b in range(a+1, m):
       for c in range(b+1, m):
        v = zeros2(N)
        for i, x in enumerate(xs):
            if x[a] == x[b] == x[c] == 1:
                v[i] += 1
        assert v.sum()
        lhs.append(v)
        rhs.append(0)

    # dissallow columns with weight <= 1
    for i, x in enumerate(xs):
        if sum(x)<=1:
            v = zeros2(N)
            v[i] = 1
            lhs.append(v)
            rhs.append(0)

    if k is not None:
      # constrain to k _number of odd-weight rows
      assert 0<=k<m
      for a in range(m):
        v = zeros2(N)
        for i, x in enumerate(xs):
          if x[a] == 1:
            v[i] = 1
        lhs.append(v)
        if a<k:
            rhs.append(1)
        else:
            rhs.append(0)

    A = array2(lhs)
    rhs = array2(rhs)

    if verbose:
        print("lhs:")
        print(shortstr(A))
    
        print("rhs:")
        print(shortstr(rhs))

    B = pseudo_inverse(A)
    soln = dot2(B, rhs)
    if not eq2(dot2(A, soln), rhs):
        print("no solution")
        return
    if verbose:
        print("soln:")
        print(shortstr(soln))

    soln.shape = (N, 1)
    rhs.shape = A.shape[0], 1

    K = array2(list(find_kernel(A)))

    best = None
    density = 1.0
    size = 9999*n
    trials = argv.get("trials", 1024)
    count = 0
    for trial in range(trials):
        u = rand2(len(K), 1)
        v = dot2(K.transpose(), u)
        #print(v)
        assert dot2(A, v).sum()==0
        #if v.sum() != n:
        #    continue
        assert v[0]==0
        v = (v+soln)%2
        assert eq2(dot2(A, v), rhs)

        Gt = list(G0.transpose())
        for i, x in enumerate(xs):
            if v[i] and not i in taken:
                Gt.append(x)
        if not Gt:
            continue
        Gt = array2(Gt)
        G = Gt.transpose()
        if verbose:
            print("G0")
            print(shortstr(G0))
            print("solution:")
            print(shortstr(G))
        assert is_morthogonal(G, 3)
        if G.shape[1]<m:
            continue

        if 0 in G.sum(1):
            continue

        #print(shortstr(G))
#        for g in G:
#            print(shortstr(g), g.sum())
#        print()

        _density = float(G.sum()) / (G.shape[0]*G.shape[1])
        #if best is None or _density < density:
        if best is None or G.shape[1] < size:
            best = G
            density = _density
            size = G.shape[1]

        if 0:
            #sols.append(G)
            Gx = even_rows(G)
            assert is_morthogonal(Gx, 3)
            if len(Gx)==0:
                continue
            GGx = array2(list(span(Gx)))
            assert is_morthogonal(GGx, 3)

        count += 1

    print("found %d solutions" % count)

    G = best
    #print(shortstr(G))
    for g in G:
        print(shortstr(g), g.sum())
    print()
    print("density:", density)

    #print(dot2(G, G.transpose())) # yes it's self-dual
    

def search_selfdual():

    verbose = argv.get("verbose")
    m = argv.get("m", 6) # _number of rows
    k = argv.get("k", None) # _number of odd-weight rows


    maxweight = argv.get("maxweight", m)
    minweight = argv.get("minweight", 1)

    # these are the variables N_x
    print("building xs...")

    if 0:
        xs = cross([(0, 1)]*m)
        xs = [x for x in xs if minweight <= sum(x) <= maxweight]
    
        prune = argv.get("prune", 0.5)
        xs = [x for x in xs if random() < prune]

    xs = []
    N = argv.get("N", m*100)
    colweight = argv.get("colweight", maxweight)
    assert colweight <= m
    for i in range(N):
        x = [0]*m
        total = 0
        while total < colweight:
            idx = randint(0, m-1)
            if x[idx] == 0:
                x[idx] = 1
                total += 1
        xs.append(x)

    N = len(xs)

    lhs = []
    rhs = []

    # bi-orthogonality
    for a in range(m):
      for b in range(a+1, m):
        v = zeros2(N)
        for i, x in enumerate(xs):
            if x[a] == x[b] == 1:
                v[i] = 1
        if v.sum():
            lhs.append(v)
            rhs.append(0)

    k = 0 # all rows must have even weight
    # constrain to k _number of odd-weight rows
    assert 0<=k<m
    for a in range(m):
      v = zeros2(N)
      for i, x in enumerate(xs):
        if x[a] == 1:
          v[i] = 1
      lhs.append(v)
      if a<k:
          rhs.append(1)
      else:
          rhs.append(0)

    logops = argv.logops

    A = array2(lhs)
    rhs = array2(rhs)
    #print(shortstr(A))

    print("solve...")
    B = pseudo_inverse(A)
    soln = dot2(B, rhs)
    if not eq2(dot2(A, soln), rhs):
        print("no solution")
        return

    if verbose:
        print("soln:")
        print(shortstr(soln))

    soln.shape = (N, 1)
    rhs.shape = A.shape[0], 1

    K = array2(list(find_kernel(A)))
    print("kernel:", K.shape)
    if len(K)==0:
        return
    #print(K)
    #print( dot2(A, K.transpose()))
    #sols = []
    #for v in span(K):
    best = None
    density = 1.0
    size = 99*N
    trials = argv.get("trials", 1024)
    count = 0
    print("trials...")
    for trial in range(trials):
        u = rand2(len(K), 1)
        v = dot2(K.transpose(), u)
        #print(v)
        v = (v+soln)%2
        assert eq2(dot2(A, v), rhs)

        if v.sum() >= size:
            continue

        if v.sum() < m:
            continue

        if v.sum():
            print(v.sum(), end=" ", flush=True)

        size = v.sum()

        if logops is not None and size != 2*m+logops:
            continue

        Gt = []
        for i, x in enumerate(xs):
            if v[i]:
                Gt.append(x)

        Gt = array2(Gt)
        G = Gt.transpose()
        if dot2(G, Gt).sum() != 0:
            # not self-dual
            print(shortstr(dot2(G, Gt)))
            assert 0
            return

        #if G.shape[1]<m:
        #    continue

        if 0 in G.sum(1):
            print(".", end="", flush=True)
            continue

        #print(shortstr(G))
#        for g in G:
#            print(shortstr(g), g.sum())
#        print()

        _density = float(G.sum()) / (G.shape[0]*G.shape[1])
        #if best is None or _density < density:
        if best is None or G.shape[1] <= size:
            best = G
            size = G.shape[1]
            density = _density

        if 0:
            #sols.append(G)
            Gx = even_rows(G)
            assert is_morthogonal(Gx, 3)
            if len(Gx)==0:
                continue
            GGx = array2(list(span(Gx)))
            assert is_morthogonal(GGx, 3)

        count += 1

    print("found %d solutions" % count)
    if best is None:
        return

    G = best
    #print(shortstr(G))
    f = open("selfdual.ldpc", "w")
    for spec in ["Hx =", "Hz ="]:
        print(spec, file=f)
        for g in G:
            print(shortstr(g), file=f)
    f.close()

    print()
    print("density:", density)
    print("shape:", G.shape)
    

    if 0:
        B = pseudo_inverse(A)
        v = dot2(B, rhs)
        print("B:")
        print(shortstr(B))
        print("v:")
        print(shortstr(v))
        assert eq2(dot2(B, v), rhs) 


def even_rows(G):
    Gx = []
    for u in G:
        #print(shortstr(u), u.sum()%2)
        parity = u.sum()%2
        if parity==0:
            Gx.append(u)
    Gx = array2(Gx)
    return Gx


def get_code():

    name = argv.get("code")
    code = None
    if name == "toric":
        G = parse("""
        1.1.....
        .1...1..
        11.11...
        .111..1.
        1...11.1
        """) # l=2 toric code X logops + X stabs 

    elif name == "toric3":
        G = parse("""
        .....1.....1.....1
        ............1.1.1.
        .111..........1...
        ...111..........1.
        1.....11...1......
        ..1....111........
        ....1....111......
        ......1.....11...1
        ........1....111..
        ..........1....111
        """) # l=3 toric code X logops + X stabs 

    elif name == "steane":
        G = parse("""
        1111111
        1111...
        .1.11.1
        ..11.11
        """)

    elif name == "haah":
        # triorthogonal matrix
        G = parse("""
        1111111.......
        .......1111111
        1.1.1.11.1.1.1
        .11..11.11..11
        ...1111...1111
        """)
        G = parse("""
        1.1.1.11.1.1.1
        .11..11.11..11
        ...1111...1111
        """)

        G = parse("""
        ..............1111111111111111
        ......11111111........11111111
        ..1111....1111....1111....1111
        .1..11..11..11..11..11..11..11
        11.1.1.1.1.1.1.1.1.1.1.1.1.1.1
        """)

    elif name == "rm":
        r = argv.get("r", 1)
        m = argv.get("m", 5)
        code = reed_muller(r, m)
        if argv.puncture:
            code = code.puncture(0)

    else:
        return None

    if code is None:
        code = Code(G)

    if argv.dual:
        code = code.get_dual()

    return code


def triortho():
    code = get_code()

    code.dump()
    print(code)

    Gx = []
    for u in code.G:
        print(shortstr(u), u.sum()%2)
        parity = u.sum()%2
        if parity==0:
            Gx.append(u)
    Gx = array2(Gx)

    print("is_triorthogonal:", code.is_triorthogonal())

    A = array2(list(span(Gx)))
    print("span(Gx) is_morthogonal(2):", is_morthogonal(A, 2))
    print("span(Gx) is_morthogonal(3):", is_morthogonal(A, 3))

    return

    G = code.G

#    A = array2(list(span(G)))
#    poly = {}
#    for v in A:
#        w = v.sum()
#        poly[w] = poly.get(w, 0) + 1
#    print(poly)

    k, n = G.shape

    if 0:
        from comm import Poly
        a = Poly({(1,0):1})
        b = Poly({(0,1):1})
        poly = Poly.zero(2)
        for v in span(G):
            w = v.sum()
            term = Poly({(n-w,0) : 1}) * Poly({(0,w) : 1})
            poly = poly + term
        print(poly)

    # print higher genus weight enumerator
    genus = argv.get("genus", 1)
    assert 1<=genus<=4
    N = 2**genus
    idxs = list(cross([(0,1)]*genus))

    cs = {} # _coefficients : map exponent to coeff
    for vs in cross([list(span(G)) for _ in range(genus)]):
        key = [0]*N
        for i in range(n):
            ii = tuple(v[i] for v in vs)
            idx = idxs.index(ii)
            key[idx] += 1
        key = tuple(key)
        cs[key] = cs.get(key, 0) + 1
    #print(cs)
    keys = list(cs.keys())
    keys.sort()
    print(idxs)
    for key in keys:
        print(key, cs[key])
            

#def test_full_weight():
#
#    from vec import Space, Hom, Map
#
#    import element
#
#    ring = element.Z
#    #ring = element.Q
#    one = ring.one
#
#    space = Space(2, ring)
#    hom = Hom(space, space)
#
#    I = Map.from_array([[1, 0], [0, 1]], hom)
#    X = Map.from_array([[0, 1], [1, 0]], hom)
#    Z = Map.from_array([[1, 0], [0, -1]], hom)
#    Y = Z*X
#
#    stab = [I@I, X@X]
#    for g in stab:
#      for h in stab:
#        assert g*h == h*g
#    G = mulclose(stab)
#    assert -(I@I) not in G
#    print("|G| =", len(G))
#
#    fwe = lambda I, Y, X, Z : I@I + X@X
#    fwe_d = lambda I, Y, X, Z : I@I + X@X + Z@Z + Y@Y
#
#    P = fwe_d(I, Y, X, Z)
#    print("P=")
#    print(P)
#    print("P*P=")
#    print(P*P)
#    #assert (P*P) == len(G)*P
#
#    x, y, z, t = I, Y, X, Z
#    Q = fwe(x+y+z+t, x+y-z-t, x-y+z-t, x-y-z+t)
#    print(Q)
#    assert Q == len(G) * P # FAIL
#
#    print("OK")


class Op(object):
    def __init__(self, desc, sign=1):
        self.desc = desc
        for s in desc:
            assert s in "IXZY"
        self.sign = sign
        self.n = len(desc)

    def __mul__(self, other):
        assert self.n == other.n
        sign = self.sign * other.sign
        res = []
        for a, b in zip(self.desc, other.desc):
            ab = a+b
            if a==b:
                res.append("I")
            elif a=="I":
                res.append(b)
            elif b=="I":
                res.append(a)
            elif ab=="XZ":
                res.append("Y") # Y == XZ
            elif ab=="XY":
                res.append("Z") # Y == XZ
            elif ab=="YX":
                res.append("Z")
                sign *= -1
            elif ab=="ZX":
                res.append("Y")
                sign *= -1
            elif ab=="ZY":
                res.append("X")
                sign *= -1
            elif ab=="YZ":
                res.append("X")
            else:
                assert 0
        op = Op("".join(res), sign)
        return op

    def __neg__(self):
        return Op(self.desc, -self.sign)

    def __str__(self):
        return self.desc
    __repr__ = __str__

    def __hash__(self):
        return hash((self.desc, self.sign))

    def __eq__(self, other):
        return (self.desc, self.sign) == (other.desc, other.sign)

    def __ne__(self, other):
        return (self.desc, self.sign) == (other.desc, other.sign)

    def build(self, ns):
        A = None
        for g in self.desc:
            op = ns[g]
            if A is None:
                A = op
            else:
                A = A@op # tensor
        return self.sign*A


def test_full_weight():


    I = Op("I")
    X = Op("X")
    Z = Op("Z")
    assert I*X == X
    assert X*X == I
    assert Z*X == -X*Z

    def build():
        if argv.steane:
            zops = "1111111 1001011 0101101 0010111".split()
            xops = "1001011 0101101 0010111".split()
        elif argv.toric:
            zops = "111..1.. 1.11...1 .1..111.".split()
            #xops = "1.1..... .1...1.. 11.11... .111..1. 1...11.1".split()
            xops = "11.11... .111..1. 1...11.1".split()
        else:
            return

        stab = []
        for s in xops:
            s = s.replace("0", "I")
            s = s.replace(".", "I")
            gx = s.replace("1", "X")
            stab.append(Op("".join(gx)))
        for s in zops:
            s = s.replace("0", "I")
            s = s.replace(".", "I")
            gz = s.replace("1", "Z")
            stab.append(Op("".join(gz)))

        for g in stab:
          for h in stab:
            assert g*h == h*g
        G = mulclose(stab)
        print("|G| =", len(G))
    
        return G

    def shortstr(P):
        s = str(P)
        s = s.replace("0", ".")
        s = s.replace(" ", "")
        return s

    G = build()
    #print(G)

    from vec import Space, Hom, Map
    import element
    ring = element.Z
    space = Space(2, ring)
    hom = Hom(space, space)
    I = Map.from_array([[1, 0], [0, 1]], hom)
    X = Map.from_array([[0, 1], [1, 0]], hom)
    Z = Map.from_array([[1, 0], [0, -1]], hom)
    Y = X*Z

    ns = {"I":I, "X":X, "Z":Z, "Y":Y}
    ops = [g.build(ns) for g in G]

    # P is projector onto span(logical |0>, logical |1>)
    P = ops[0]
    for op in ops[1:]:
        P = P+op
    assert P*P == len(G)*P

    x, y, z, t = I, Y, X, Z
    ns = {"I":x+y+z+t, "Y":x+y-z-t, "X":x-y+z-t, "Z":x-y-z+t}
    ops = [g.build(ns) for g in G]

    # Q is projector onto logical |0>
    Q = ops[0]
    for op in ops[1:]:
        Q = Q+op
    Q = Q/1024
    #assert Q*Q == 8*Q

    import elim
    def astr(A):
        s = str(elim.shortstr(A))
        s = s.replace(" 0 ", " . ")
        s = s.replace(" ", "")
        return s

    Q = elim.row_reduce(ring, Q.to_array(), truncate=True)
    print("Q:")
    print(astr(Q))

    #s = astr(Q)
    #for i, c in enumerate(s):
    #    if c=="1":
    #        print(bin(i-2))

    P = elim.row_reduce(ring, P.to_array(), truncate=True)
    print("P:")
    print(astr(P))

    return

    x, y, z, t = I, Y, X, Z
    I, Y, X, Z = x+y+z+t, x+y-z-t, x-y+z-t, x-y-z+t



    P = fwe_d(I, Y, X, Z)
    print("P=")
    print(P)
    print("P*P=")
    print(P*P)
    #assert (P*P) == len(G)*P

    x, y, z, t = I, Y, X, Z
    Q = fwe(x+y+z+t, x+y-z-t, x-y+z-t, x-y-z+t)
    print(Q)
    assert Q == len(G) * P

    print("OK")



if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        from util import set_seed
        set_seed(_seed)

    name = argv.next()
    if name is None:
        test()
    else:
        print(name+"()")
        fn = eval(name)
        fn()
        print("OK.")



