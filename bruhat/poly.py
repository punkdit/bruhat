#!/usr/bin/env python3

"""
Multivariate commutative polynomials.
"""

from functools import reduce
import operator
from random import randint, shuffle, seed

import numpy

from bruhat.argv import argv
from bruhat.util import cross, factorial, choose, determinant
from bruhat.theta import divisors
from bruhat.element import Fraction, Q
from bruhat import series


#def is_scalar(x):
#    return isinstance(x, int)

def tpl_zip(left, right):
    i = j = 0
    while 1:
        if i < len(left) and j < len(right):
            kl, vl = left[i]
            kr, vr = right[j]
            if kl < kr:
                yield (kl, vl, 0)
                i += 1
            elif kl == kr:
                yield (kl, vl, vr)
                i += 1
                j += 1
            else:
                assert kl > kr
                yield (kr, 0, vr)
                j += 1
        elif i < len(left):
            kl, vl = left[i]
            yield (kl, vl, 0)
            i += 1
        elif j < len(right):
            kr, vr = right[j]
            yield (kr, 0, vr)
            j += 1
        else:
            break


def tpl_add(left, right):
    result = tuple((k, l+r) for (k, l, r) in tpl_zip(left, right) if l+r)
    #assert result == tpl_add_old(left, right)
    return result

def tpl_sub(left, right):
    result = tuple((k, l-r) for (k, l, r) in tpl_zip(left, right) if l-r)
    return result

def tpl_union(left, right):
    result = tuple((k, max(l,r)) for (k, l, r) in tpl_zip(left, right))
    return result

def tpl_disjoint(left, right):
    for (k, l, r) in tpl_zip(left, right):
        if l and r:
            return False
    return True

def tpl_ge(left, right):
    for (k, l, r) in tpl_zip(left, right):
        if l<r:
            return False
    return True

def tpl_gt(left, right):
    "total order on monomials"
    if left==right:
        return False
    l = sum([l[1] for l in left], 0)
    r = sum([r[1] for r in right], 0)
    if l != r:
        return l>r
    for (k, l, r) in tpl_zip(left, right):
        if l != r:
            return l>r
    assert 0 # should not get here...



def tpl_degree(tpl):
    d = reduce(operator.add, [v for (k,v) in tpl], 0)
    #print("tpl_degree", tpl, "=", d)
    vals = [k for (k,v) in tpl]
    vals.sort()
    return d, tuple(vals)

def tex_str(v):
    if isinstance(v, Fraction):
        if v.bot == 1:
            return str(v.top)
        else:
            return r"\frac{%s}{%s}" % (v.top, v.bot)
    else:
        return str(v)

def shortstr(v):
    return str(v)


class Poly(object):

    def __init__(self, cs, ring):
        if isinstance(cs, str):
            cs = {((cs,1),):ring.one}
        coefs = {}
        keys = []
        degree = 0
        head = ()
        for key, value in cs.items():
            if value == 0:
                continue
            key = list(key)
            key.sort()
            key = tuple(key)
            if tpl_gt(key, head):
                head = key
            dkey = dict(key)
            assert len(dkey) == len(key)
            deg = 0
            for v in dkey.values():
                assert v, "got some zero exponents: %s"%repr(dkey)
                deg += v
            degree = max(degree, deg)
            coefs[key] = ring.promote(value)
            keys.append(key)
        #assert len(coefs) == len(cs), (coefs, cs)
        keys.sort()
        self.keys = keys
        self.cs = coefs
        self.degree = degree
        self.head = head
        self.ring = ring

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, idx):
        key = self.keys[idx]
        val = self.cs[key]
        return Poly({key : val}, self.ring)

    def get(self, key):
        return self.cs.get(key, self.ring.zero)

    def terms(self):
        cs = self.cs
        for key in self.keys:
            val = self.cs[key]
            yield Poly({key : val}, self.ring)

    def is_const(self):
        cs = self.cs
        return not cs or list(cs.keys()) == [(),]

    def get_const(self):
        return self.cs.get((), self.ring.zero)

    def get_vars(self):
        vs = set()
        for key in self.keys:
            for item in key:
                vs.add(item[0])
        vs = list(vs)
        vs.sort()
        return vs

    @classmethod
    def promote(cls, item, ring):
        if isinstance(item, Poly):
            return item
        #assert is_scalar(item)
        item = ring.promote(item)
        return Poly({() : item}, ring)

    def get_zero(self):
        return Poly({():self.ring.zero}, self.ring)

    def __eq__(self, other):
        other = self.promote(other, self.ring)
        return self.cs == other.cs

    def __ne__(self, other):
        other = self.promote(other, self.ring)
        return self.cs != other.cs

    def __hash__(self):
        cs = list(self.cs.items())
        cs.sort(key = lambda k_v:k_v[0])
        cs = tuple(cs)
        return hash(cs)

    def __add__(self, other):
        other = self.promote(other, self.ring)
        cs = dict(self.cs)
        ring = self.ring
        for key, value in list(other.cs.items()):
            cs[key] = cs.get(key, 0) + value
        return Poly(cs, ring)
    __radd__ = __add__

    def __sub__(self, other):
        other = self.promote(other, self.ring)
        cs = dict(self.cs)
        ring = self.ring
        for key, value in list(other.cs.items()):
            cs[key] = cs.get(key, 0) - value
        return Poly(cs, ring)

    def __neg__(self):
        cs = {}
        ring = self.ring
        for key, value in list(self.cs.items()):
            cs[key] = -value
        return Poly(cs, ring)

    def __rmul__(self, r):
        #assert is_scalar(r)
        ring = self.ring
        r = ring.promote(r)
        cs = {}
        for key, value in list(self.cs.items()):
            cs[key] = r * value
        return Poly(cs, ring)

    def __mul__(self, other):
        ring = self.ring
        value = ring.promote(other)
        #if is_scalar(other):
        if value is not None:
            return self.__rmul__(value)
        cs = {}
        for k1, v1 in list(self.cs.items()):
          for k2, v2 in list(other.cs.items()):
            k = tpl_add(k1, k2)
            cs[k] = cs.get(k, 0) + v1*v2
        return Poly(cs, ring)

    def reduce1(P, Q):
        "return R, S such that Q=R*P+S "
#        print("reduce1(%s, %s)"%(P, Q))
        m0 = P.head

        best = None
        #for m1 in Q.cs:
        for m1 in Q.keys:
            if tpl_ge(m1, m0):
              if best is None or tpl_gt(m1, best):
                best = m1
        if best is None:
            return P.get_zero(), Q

#        print("reduce1: m0=%s, best=%s"%(m0, best))
        m2 = tpl_sub(best, m0)
#        print("reduce1: best-m0 = ", m2)
        ring = P.ring
        m2 = Poly({m2:ring.one}, ring)
#        print("--->", Q.get(best), P.get(m0), Q.get(best)/P.get(m0))
        R = (Q.get(best)/P.get(m0)) * m2
        S = Q - R * P
        assert S.degree <= Q.degree
#        print("\tR=%s, S=%s"%(R, S))
#        print()
        assert S.get(best) == 0
        return R, S

    def reduce(P, Q):
        "return R, S such that Q=R*P+S "
#        print("reduce(%s, %s)"%(P, Q))
        if P==0:
            return 0, Q
        R, S = P.reduce1(Q)
        assert Q == R*P + S
        while R and S:
            R1, S1 = P.reduce1(S)
            S = S1
            R = R + R1
            assert Q == R*P + S
            if R1 == 0:
                break
        return R, S

    def __div__(self, other):
        other = self.promote(other, self.ring)
        R, S = other.reduce(self)
        if S!=0:
            raise Exception("%s not divisible by %s" % (self, other))
        return R
    __truediv__ = __div__

    def normalized(P):
        "divide by leading coefficient"
        head = P.head
        c = P.get(head)
        if c != P.ring.one:
            P = P / c
        return P

    def critical_pair(P, Q):
        P = P.normalized()
        Q = Q.normalized()
        if tpl_disjoint(P.head, Q.head):
            return P.ring.zero
        m = tpl_union(P.head, Q.head)
        mi = tpl_sub(m, P.head)
        mj = tpl_sub(m, Q.head)
        ring = P.ring
        mi = Poly({mi:1}, ring)
        mj = Poly({mj:1}, ring)
        S = mi * P - mj * Q
        return S

    def diff(self, var):
        "_differentiate wrt var"
        if isinstance(var, Poly):
            v = str(var)
            assert var == Poly(v, self.ring)
            var = v
        #print("diff", self, ",", var)
        cs = {}
        for k, coeff in self.cs.items():
            #print(k, coeff)
            match = (var, 0)
            rest = []
            for item in k:
                if item[0] == var:
                    match = item
                else:
                    rest.append(item)
            #print("match:", match[0], match[1])
            deg = match[1]
            if deg == 0:
                continue
            if deg > 1:
                rest.append((var, deg-1))
                coeff *= deg
            rest = tuple(rest)
            cs[rest] = coeff
        return Poly(cs, self.ring)

    def __pow__(self, n):
        ring = self.ring
        if n==0:
            return Poly({():ring.one}, ring)
        assert n>0, repr(n)
        p = self
        for i in range(n-1):
            p = self*p
        return p

    def str(self, shortstr=tex_str, POW="^", MUL="*", OPEN="{", CLOSE="}"):
        items = list(self.cs.items())
        if not items:
            return "0"
        items.sort(key = (lambda kv:tpl_degree(kv[0])))
        #print(items)
        ring = self.ring
        terms = []
        for (k, v) in items:
            assert v != 0, self.cs
            ss = []
            for name, exp in k:
                assert exp>0
                #if len(name)>1:
                #    name += " "
                if exp == 1:
                    ss.append(name)
                else:
                    if exp>9:
                        ss.append("%s%s%s%s%s"%(name, POW, OPEN, exp, CLOSE))
                    else:
                        ss.append("%s%s%s"%(name, POW, exp))
            s = '*'.join(ss)
            if s:
                if v==ring.one:
                    terms.append(s)
                else:
                    terms.append("%s%s%s" % (shortstr(v), MUL, s))
            else:
                if v==ring.one:
                    terms.append("1")
                else:
                    terms.append(shortstr(v))
        s = " + ".join(terms)
        #s = s.replace("-1"+MUL, "-") # oof, careful with this!
        s = s.replace("-1*", "-")
        s = s.replace("+ -", "- ")
        s = s.replace(" "+POW, POW)
        if s and s[-1] == " ":
            s = s[:-1]
        return s

    def __str__(self):
        s = self.str(shortstr)
        return s
    __repr__ = __str__

    def python_str(self):
        s = self.str(shortstr, "**", "*", "", "")
        s = s.replace("{", "")
        s = s.replace("}", "")
        return s

    def otherstr(self, mul=""):
        cs = self.cs
        keys = set()
        for key in self.keys:
            for k,v in key:
                keys.add(k)
        keys = list(keys)
        keys.sort()
        n = len(keys)
        ss = []
        for key in self.keys:
            value = cs[key]
            items = [0] * n
            for k, v in key:
                items[keys.index(k)] = v
            if value==1:
                s = "%s" % (tuple(items),)
            else:
                s = "%s%s%s" % (value, mul, tuple(items))
            s = s.replace(' ', '')
            ss.append(s)
        return ' + '.join(ss) or "0"

    def flatstr(self, mul=""):
        cs = self.cs
        keys = set()
        for key in self.keys:
            for k,v in key:
                keys.add(k)
        keys = list(keys)
        keys.sort()
        n = len(keys)
        #print(keys)
        entries = []
        for key in self.keys:
            value = cs[key]
            items = [0] * n
            for k, v in key:
                items[keys.index(k)] = v
            entries.append((items, value))
        def sortfunc(entry):
            _items, value = entry
            items = list(_items)
            items.sort()
            nz = [i for i in _items if i]
            return (items, value, nz)
        entries.sort(key = sortfunc)

        ss = []
        for items, value in entries:
            if value==1:
                s = "%s" % (tuple(items),)
            else:
                s = "%s%s%s" % (value, mul, tuple(items))
            s = s.replace(' ', '')
            ss.append(s)
        return ' + '.join(ss) or "0"

    def qstr(self, name="q"):
        keys = [()] + [((name, i),) for i in range(1, self.degree+1)]
        items = [self.cs.get(key) for key in keys]
        if max(items)>9:
            s = ",".join([str(i) for i in items])
        else:
            s = "".join([str(i) for i in items])
        return s

    def substitute(self, ns): # try not to use this :P
        s = self.python_str()
        assert type(ns) is dict, repr(ns)
        p = eval(s, ns) # argh, boo, lame!
        return p

    def __call__(self, **kw):
        #return self.substitute(kw)
        ring = self.ring
        #val = ring.zero
        val = self
        for (k, v) in kw.items():
            p = Poly(k, ring)
            val = (p-v).reduce(val)[1]
        if val.is_const():
            val = val.get_const()
        return val

    @classmethod
    def random(cls, vars, ring, degree=3, terms=3):
        cs = {}
        rank = len(vars)
        for i in range(terms):
            key = [0]*rank
            remain = degree
            for j in range(rank):
                key[j] = randint(0, remain)
                remain -= key[j]
            shuffle(key)
            key = tuple((vars[i], key[i]) for i in range(rank) if key[i])
            value = randint(-2, 2)
            cs[key] = value
        P = cls(cs, ring)
        return P


def reduce_many(polys, P):
    seen = set([P])
    while 1:
#        print("reduce_many", len(seen))
        for Q in polys:
#            print(Q, "reduce", P)
            Q1, P1 = Q.reduce(P)
#            print("\t", Q1, "||", P1)
            P = P1
        if P in seen:
            break
#        print("\t", P)
        seen.add(P)
        #assert len(seen) < 5
    return P


def grobner(polys, reduce=True, verbose=False):
    "Build a grobner basis from polys."

    assert polys
    ring = polys[0].ring
    polys = [p.normalized() for p in polys if p!=0]
    if verbose:
        print("grobner(%s)"%(polys,))

    basis = set(polys)
    pairs = []
    for i in range(len(polys)):
      for j in range(i+1, len(polys)):
        pairs.append((polys[i], polys[j]))

    seen = set()
    while pairs:
        pairs.sort(key = lambda pair:len(str(pair)))
#        print("grobner", len(pairs))
        P, Q = pairs.pop(0)
        if (P, Q) in seen:
            continue
        seen.add((P, Q))
        S1 = P.critical_pair(Q)
        if S1==0:
            continue

        S = reduce_many(basis, S1)
        if S==0:
            continue

        S = S.normalized()
        if verbose:
            print("grobner: P=%s, Q=%s, S(P,Q)=%s --> %s"%(P, Q, S1, S))

        if S in basis:
            continue

        for P in basis:
            if (P, S) not in seen:
                pairs.append((P, S))
        basis.add(S)


    if reduce:
        new = set()
        for p in basis:
            if new:
                p = reduce_many(new, p)
            if p != 0:
                new.add(p)
        basis = new

#        for p in basis:
#          for q in basis:
#            if q is p:
#                continue
#            r = p.reduce(q)
#            if r==0:
#                assert 0, (p, q, r)

    basis = list(basis)
    basis.sort(key = str)

    return basis



def test_sl2():

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    a, b, c, d, e, f, g, h = [Poly(c, ring) for c in 'abcdefgh']

    def repr_2d(a, b, c, d):
        A = numpy.empty((2, 2), dtype=object)
        A[:] = [[a, b], [c, d]]
        #A = A.transpose()
        return A

    def repr_3d(a, b, c, d):
        A = numpy.empty((3, 3), dtype=object)
        A[:] = [[a*a, a*c, c*c], [2*a*b, a*d+b*c, 2*c*d], [b*b, b*d, d*d]]
        A = A.transpose()
        return A

    #a, b, c, d, e, f, g, h = [1, 1, 1, 2, 1, 2, 0, 1]

    dot = numpy.dot
    A2 = repr_2d(a, b, c, d)
    B2 = repr_2d(e, f, g, h)
    #print(A2)
    #print(B2)
    AB2 = dot(A2, B2)
    #print(AB2)
    a0, b0, c0, d0 = AB2.flat

    A3 = repr_3d(a, b, c, d)
    B3 = repr_3d(e, f, g, h)
    AB3 = dot(A3, B3)
    #print(A3)
    #print(B3)
    #print(AB3)
    #print(repr_3d(a0, b0, c0, d0))
    assert numpy.alltrue(AB3 == repr_3d(a0, b0, c0, d0))


def test():

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)
    xx = Poly({(("x", 2),) : 1}, ring)
    xy = Poly({(("x", 1), ("y", 1)) : 1}, ring)

    assert str(one) == "1"
    assert str(one+one) == "2"
    assert zero*x == zero
    assert x+zero == x
    assert one*x == x
    assert (one + one) * x == x+x
    assert 2*x == x+x
    assert x*x == xx
    assert x*y == y*x == xy
    assert one.degree == 0
    assert x.degree == 1
    assert xy.degree == 2

    p = (x+y+1)**3
    #print(p)
    #return
    assert reduce(operator.add, p.terms()) == p

    assert eval(p.python_str(), locals()) == p

    a = zero
    b = zero
    for i in range(5):
        a += Poly("a_%d"%i, ring)
        b += Poly("b_%d"%i, ring)

    p = x**2 * y + x*y - 17
    q = x**2 * y**2 + x*y**2 - 3
    _, R = p.reduce(q)
    assert R == 17*y - 3

    div, rem = y.reduce(x*y+y)
    assert rem == 0
    assert div == x+1
    assert (x*y + y) / y == (x + 1)
    assert (x**2 - y**2) / (x+y) == (x-y)

    # -------------------------------------
    # 

    for trial in range(100):
        p = Poly.random(list('xyz'), ring)
        q = Poly.random(list('xyz'), ring)
        #print p, "||", q
        r, s = p.reduce(q)

        if p==0:
            continue
        pq = p*q
        r, s = p.reduce(pq)
        assert s==0
        assert r==q

    #return

    # -------------------------------------
    # 

    p = x**2*y-x**2
    q = x*y**2-y**2
    assert p.critical_pair(q) == x*y**2 - x**2*y
    assert p.critical_pair(q) == -q.critical_pair(p)
    assert p.critical_pair(p) == 0
    assert q.critical_pair(q) == 0


    # -------------------------------------
    # 

    x1 = Poly("x1", ring)
    x2 = Poly("x2", ring)
    p = x1**2*x2-x1**2
    q = x1*x2**2-x2**2
    polys = [p, q]
    basis = grobner(polys, verbose=True)
    print(basis)

    # -------------------------------------
    # 

    p = x**2*y-x**2
    q = x*y**2-y**2
    polys = [p, q]
    basis = grobner(polys)
    #print(basis)

    # yikes, easy for this to go off the rails!
    for i in range(12):
        polys = [Poly.random(list('xyz'), ring) for _ in range(3)]
        print(polys)
        basis = grobner(polys)
        print(len(basis))

    polys = [x+y+z, x*y+x*z+y*z, x*y*z]
    basis = grobner(polys)
    #for p in basis:
    #    print(p)
    #print(basis)



dim_sl3 = lambda a, b : (a+1)*(b+1)*(a+b+2)//2

def test_hilbert_sl3():
    # ---------------------------------------------------------------
    # Here we work out the coefficients of a Hilbert polynomial
    # given by the rational function top/bot.
    # These coefficients give the dimension of irreps of SL(3).
    # See also test_graded_sl3 below.

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)

    top = one - x*y
    bot = ((one-x)**3) * ((one-y)**3)

    def diff(top, bot, var):
        top, bot = top.diff(var)*bot - bot.diff(var)*top, bot*bot
        return top, bot

    fracs = {}
    fracs[0,0] = (top, bot)

    N = 3
    for i in range(N):
      for j in range(N):
        top, bot = fracs[i, j]
        top, bot = diff(top, bot, 'x')
        fracs[i, j+1] = top, bot

        top, bot = fracs[i, j]
        top, bot = diff(top, bot, 'y')
        fracs[i+1, j] = top, bot
        print(".", end="", flush=True)
    print()

    for i in range(N+1):
      for j in range(N+1):
        if (i, j) not in fracs:
            continue
        top, bot = fracs[i, j]
        t = top.get_const()
        b = bot.get_const()
        val = t//b//factorial(i)//factorial(j)
        assert val == dim_sl3(i, j)
        print(val, end=" ")
      print()


def all_monomials(vs, deg, ring):
    n = len(vs)
    assert n>0
    if n==1:
        v = vs[0]
        yield v**deg
        return

    items = list(range(deg+1))
    for idxs in cross((items,)*(n-1)):
        idxs = list(idxs)
        remain = deg - sum(idxs)
        if remain < 0:
            continue
        idxs.append(remain)
        p = Poly({():1}, ring)
        assert p==1
        for idx, v in zip(idxs, vs):
            p = p * v**idx
        yield p
            


def test_graded_sl3():
    # ---------------------------------------------------------------
    # slightly more explicit calculation than test_hilbert above

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)
    u = Poly("u", ring)
    v = Poly("v", ring)
    w = Poly("w", ring)

    rel = x*u + y*v + z*w
    rels = [rel]
    rels = grobner(rels)

    for a in range(4):
      for b in range(4):
        gens = []
        for p in all_monomials([x, y, z], a, ring):
          for q in all_monomials([u, v, w], b, ring):
            rem = p*q
            #gens.append(pq)
            for rel in rels:
                div, rem = rel.reduce(rem)
                #print(pq, div, rem)
            gens.append(rem)
    
        basis = grobner(gens)
        assert len(basis) == dim_sl3(a, b)
        print(len(basis), end=' ', flush=True)
      print()




def test_graded_sl4():
    # See: Miller & Sturmfels, p276

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    poly = lambda v : Poly(v, ring)

    p1 = poly("p1")
    p2 = poly("p2")
    p3 = poly("p3")
    p4 = poly("p4")
    p12 = poly("p12")
    p13 = poly("p13")
    p14 = poly("p14")
    p23 = poly("p23")
    p24 = poly("p24")
    p34 = poly("p34")
    p123 = poly("p123")
    p124 = poly("p124")
    p134 = poly("p134")
    p234 = poly("p234")

    rels = [
        p23*p1 - p13*p2 + p12*p3,       p24*p1 - p14*p2 + p12*p4,
        p34*p1 - p14*p3 + p13*p4,       p34*p2 - p24*p3 + p23*p4,
        p14*p23 - p13*p24 + p12*p34,    p234*p1 - p134*p2 + p124*p3 - p123*p4,
        p134*p12 - p124*p13 + p123*p14, p234*p12 - p124*p23 + p123*p24,
        p234*p13 - p134*p23 + p123*p34, p234*p14 - p134*p24 + p124*p34,
    ]
    rels = grobner(rels, verbose=True)
    print("rels:", rels)
    print()

    grades = [
        [p1, p2, p3, p4],
        [p12, p13, p14, p23, p24, p34],
        [p123, p124, p134, p234],
    ]
    multi = argv.get("multi")
    n = 5 if multi is None else sum(multi)+1
    n = argv.get("n", n)
    for g0 in range(n):
     for g1 in range(n):
      for g2 in range(n):
        if multi is not None and (g0, g1, g2)!=multi:
            #print(".  ", end='')
            continue
        elif g0+g1+g2 > n-1:
            print(".  ", end='')
            continue
        gens = []
        for m0 in all_monomials(grades[0], g0, ring):
         for m1 in all_monomials(grades[1], g1, ring):
          for m2 in all_monomials(grades[2], g2, ring):
            m = m0*m1*m2
            #for rel in rels:
            #    div, m = rel.reduce(m)
            m = reduce_many(rels, m)
            if m != 0:
                gens.append(m)
            
        print(len(gens), end=':', flush=True)
        basis = grobner(gens)
        lhs = len(basis)
        rhs = (g0+1)*(g1+1)*(g2+1)*(g0+g1+2)*(g1+g2+2)*(g0+g1+g2+3)//12
        assert lhs==rhs, ("%s != %s"%(lhs, rhs))
        print(len(basis), end=' ', flush=True)

#        basis.sort(key=str)
#        heads = {}
#        for p in basis:
#            print(p.head, p)
#            heads[p.head] = p
#        print(len(heads))
#        return

      print()
     print()



def test_plucker():
    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)

    rows, cols = argv.get("rows", 2), argv.get("cols", 4)
    U = numpy.empty((rows, cols), dtype=object)
    for i in range(rows):
      for j in range(cols):
        U[i, j] = Poly("x[%d,%d]"%(i, j), ring)

    print(U)
    COLS = list(range(cols))
    w = {} # the plucker coordinates
    for idxs in choose(COLS, rows):
        V = U[:, idxs]
        #print(V)
        a = determinant(V)
        w[idxs] = a
        #print(idxs, a)

    if (rows, cols) == (2, 4):
        assert w[0,1]*w[2,3]-w[0,2]*w[1,3]+w[0,3]*w[1,2] == 0

    for idxs in choose(COLS, rows-1):
      for jdxs in choose(COLS, rows+1):
        if len(idxs) and idxs[-1] >= jdxs[0]:
            continue
        #print(idxs, jdxs)
        sign = ring.one
        rel = ring.zero
        for l in range(rows+1):
            ldxs = idxs+(jdxs[l],)
            rdxs = jdxs[:l] + jdxs[l+1:]
            rel += sign*w[ldxs]*w[rdxs]
            sign *= -1
        assert rel==0


def test_plucker_flag():
    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)

    n = argv.get("n", 4)
    U = numpy.empty((n, n), dtype=object)
    for i in range(n):
      for j in range(n):
        U[i, j] = Poly("x[%d,%d]"%(i, j), ring)

    print(U)

    N = list(range(n))
    w = {} # the plucker coordinates
    for k in range(1, n):
      for idxs in choose(N, k):
        V = U[:k, idxs]
        #print(V)
        a = determinant(V)
        if k==1:
            w[idxs[0]] = a
        else:
            w[idxs] = a
        #print(idxs, a)

    assert n==4
    p1 = w[0]
    p2 = w[1]
    p3 = w[2]
    p4 = w[3]
    p12 = w[0,1]
    p13 = w[0,2]
    p14 = w[0,3]
    p23 = w[1,2]
    p24 = w[1,3]
    p34 = w[2,3]
    p123 = w[0,1,2]
    p124 = w[0,1,3]
    p134 = w[0,2,3]
    p234 = w[1,2,3]

    for rel in [
        p23*p1 - p13*p2 + p12*p3,       p24*p1 - p14*p2 + p12*p4,
        p34*p1 - p14*p3 + p13*p4,       p34*p2 - p24*p3 + p23*p4,
        p14*p23 - p13*p24 + p12*p34,    p234*p1 - p134*p2 + p124*p3 - p123*p4,
        p134*p12 - p124*p13 + p123*p14, p234*p12 - p124*p23 + p123*p24,
        p234*p13 - p134*p23 + p123*p34, p234*p14 - p134*p24 + p124*p34,
    ]:
        assert rel == 0

    return

    for idxs in choose(N, rows-1):
      for jdxs in choose(N, rows+1):
        if len(idxs) and idxs[-1] >= jdxs[0]:
            continue
        print(idxs, jdxs)
        sign = ring.one
        rel = ring.zero
        for l in range(rows+1):
            ldxs = idxs+(jdxs[l],)
            rdxs = jdxs[:l] + jdxs[l+1:]
            rel += sign*w[ldxs]*w[rdxs]
            sign *= -1
        print(rel)


class Formal(series.Series):
    """
        Formal power series with coefficients in a polynomial ring with infinitely many variables.
    """
    def __init__(self, name, ring, subs={}):
        series.Series.__init__(self, ring)
        self.name = name
        self.subs = dict(subs)

    def get_name(self, idx):
        name = self.name
        if idx < 10:
            name = "%s_%d"%(name, idx)
        else:
            name = "%s_{%d}"%(name, idx)
        return name

    def getitem(self, idx):
        name = self.get_name(idx)
        p = Poly(name, Q)
        p = self.subs.get(name, p)
        return p


class FormalExp(Formal):
    def getitem(self, idx):
        name = self.get_name(idx)
        p = Poly(name, Q)
        p = self.subs.get(name, p)
        #p = Poly({((name, 1),) : Q.one//factorial(idx)}, Q)
        p = Q.one//factorial(idx) * p
        return p


# could package this up in another Series subclass...
def solve(f, g, fg, ring):
    # solve fg=1 for g

    soln = {}
    sf = Formal(f.name, ring)
    sg = Formal(g.name, ring)
    i = 1
    while 1:
        p = fg[i]
        if p.degree == 0:
            i += 1
            continue
        soln[sf[i].python_str()] = sf[i]
        #print("soln:", soln)
        #print("p =", p)
        b_i = sg[i]
        #print(p.cs)
        c = p.cs[((str(b_i), 1),)]
        #print("c =", c)
        rhs = b_i - (1/c)*p
        #print("eval", rhs, soln)
        rhs = eval(rhs.python_str(), dict(soln))
        soln[b_i.python_str()] = rhs

        yield b_i, rhs
        i += 1



def main():

    zero = Poly({}, Q)
    one = Poly({():Q.one}, Q)

    ring = type("ARing", (object,), {})
    ring.zero = zero
    ring.one = one

    # -----------------------------------
    # Multiplicative inverse 

    print("\n\n# Multiplicative inverse ")

    f = FormalExp("a", ring, {"a_0":1})
    g = FormalExp("b", ring, {"b_0":1})

    print("f =", f)
    print("g =", g)

    fg = f*g
    for i in range(5):
        print(fg[i])
    print()

    N = argv.get("N", 6)
    items = solve(f, g, fg, ring)
    for i in range(N):
        lhs, rhs = items.__next__()
        print(r"    %s &= %s \\" % (lhs.str(), rhs.str()))
        #print(lhs, "=", rhs)

    # -----------------------------------
    # Compositional inverse

    print("\n\n# Compositional inverse")

    f = Formal("a", ring, {"a_0":0, "a_1":1})
    g = Formal("b", ring, {"b_0":0, "b_1":1})

    print("f =", f)
    print("g =", g)

    fg = f(g)
    for i in range(5):
        print(fg[i])
    print()
    
    items = solve(f, g, fg, ring)
    for i in range(5):
        lhs, rhs = items.__next__()
        print(r"    %s &= %s \\" % (lhs.str(), rhs.str()))
        #print(lhs, "=", rhs)

    # -----------------------------------
    # Dirichlet inverse

    print("\n\n# Dirichlet inverse")

    f = Formal("a", ring, {"a_0":0, "a_1":1})
    g = Formal("b", ring, {"b_0":0, "b_1":1})

    print("f =", f)
    print("g =", g)

    fg = f.dirichlet(g)
    for i in range(5):
        print(fg[i])
    print()
    
    items = solve(f, g, fg, ring)
    for i in range(37):
        lhs, rhs = items.__next__()
        print(r"    %s &= %s \\" % (lhs.str(), rhs.str()))
        #print(lhs, "=", rhs)




if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        seed(_seed)

    profile = argv.profile
    fn = argv.next()

    if profile:
        import cProfile as profile
        profile.run("test()")

    elif fn is None:
        test()

    else:
        fn = eval(fn)
        fn()

    print("OK")

