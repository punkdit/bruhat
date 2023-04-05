#!/usr/bin/env python3

"""
Multivariate commutative polynomials.
"""

import os
from functools import reduce
import operator
from random import randint, shuffle, seed

import numpy

from bruhat.argv import argv
from bruhat.util import cross, factorial, choose, determinant
from bruhat.theta import divisors
if argv.fast:
    from bruhat._element import Fraction, Q, Ring
else:
    from bruhat.element import Fraction, Q, Ring
from bruhat import series
#from bruhat.chain import Space, Lin


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
        self._diff_cache = {}
        self._subs_cache = {}

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

    def __rsub__(self, other):
        "other - self"
        return other + self.__neg__()

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

    def diff(self, var, count=1):
        "_differentiate wrt var"
        assert count>=0
        if count == 0:
            return self
        if count > 1:
            p = self
            while count:
                p = p.diff(var) # recurse
                count -= 1
            return p
        if isinstance(var, Poly):
            v = str(var)
            assert var == Poly(v, self.ring)
            var = v
        _diff_cache = self._diff_cache
        if var in _diff_cache:
            return _diff_cache[var]
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
        result = Poly(cs, self.ring)
        _diff_cache[var] = result
        return result

    def partial(self, **kw):
        p = self
        for k, v in kw.items():
            p = p.diff(k, v)
        return p

    def __pow__(self, n):
        ring = self.ring
        if n==0:
            return Poly({():ring.one}, ring)
        assert n>0, repr(n)
        p = self
        for i in range(n-1):
            p = self*p
        return p

    def str(self, shortstr=tex_str, POW="^", MUL="", OPEN="{", CLOSE="}"):
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
                    terms.append("%s%s%s" % (shortstr(v), "*", s))
            else:
                if v==ring.one:
                    terms.append("1")
                else:
                    terms.append(shortstr(v))
        s = " + ".join(terms)
        #s = s.replace("-1"+MUL, "-") # oof, careful with this!
        s = s.replace("-1*", "-")
        assert "*" not in POW, "whoops"
        s = s.replace("*", MUL)
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

    def py_substitute(self, ns): # about 4 times slower than substitute()
        assert type(ns) is dict, repr(ns)
        vs = self.get_vars()
        for v in vs:
            if v not in ns:
                ns[v] = Poly(v, self.ring)      # cumtime from profile:
        pystr = self.python_str()               # <-- 11 seconds
        co = compile(pystr, "<string>", 'eval') # <-- 10 seconds
        val = eval(co, ns)                      # <--  3 seconds
        return val

    def substitute(self, items):
        "items: tuple of (str, val)"
        assert type(items) is tuple
        _subs_cache = self._subs_cache
        result = _subs_cache.get(items)
        if result is not None:
            return result
        ns = dict(items)
        cs = self.cs
        ring = self.ring
        result = ring.zero
        for (k,term) in cs.items():
            for (ps,exp) in k:
                p = ns.get(ps)
                if p is None:
                    p = Poly({((ps,1),):ring.one}, ring)
                    ns[ps] = p
                if exp != 1:
                    p = p**exp
                term *= p
            result += term
        _subs_cache[items] = result
        return result

    def __call__(self, **kw):
        #return self.py_substitute(kw)
        kw = tuple(kw.items())
        return self.substitute(kw)

    def __slow_call__(self, **kw):
        ring = self.ring
        val = self
        for (k, v) in kw.items():
            p = Poly(k, ring)
            val = (p-v).reduce(val)[1] # slow !!!!!
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

    p = (1-x)**5
    assert p.diff(x, 0) == 1 - 5*x + 10*x**2 - 10*x**3 + 5*x**4 - x**5 == p
    assert p.diff(x, 1) == -5 + 20*x - 30*x**2 + 20*x**3 - 5*x**4
    assert p.diff(x, 2) == 20 - 60*x + 60*x**2 - 20*x**3
    assert p.diff(x, 3) == -60 + 120*x - 60*x**2


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
    basis = grobner(polys, verbose=False)
    #print(basis)

    # -------------------------------------
    # 

    a = (x+1)*(y-1)*z
    assert a() == a
    assert a(z=1) == (x+1)*(y-1)
    assert a(y=1) == 0

    a = (x+1)*(y-1)*z
    values = [a(x=i, y=2, z=k) 
        for i in [-2, 5]
        for k in [3, 4]]
    assert values == [-3, -4, 18, 24]


    return # <----------- return

#    # -------------------------------------
#    # 
#
#    p = x**2*y-x**2
#    q = x*y**2-y**2
#    polys = [p, q]
#    basis = grobner(polys)
#    #print(basis)
#
#    # yikes, easy for this to go off the rails!
#    for i in range(12):
#        polys = [Poly.random(list('xyz'), ring) for _ in range(3)]
#        print(polys)
#        basis = grobner(polys)
#        print(len(basis))
#
#    polys = [x+y+z, x*y+x*z+y*z, x*y*z]
#    basis = grobner(polys)
#    #for p in basis:
#    #    print(p)
#    #print(basis)



def test_genus_2():
    # output some povray code
    # https://math.stackexchange.com/a/349914
    # https://www.povray.org/documentation/3.7.0/r3_4.html#r3_4_5_3_5

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)

    p =  ((x**2+y**2)**2 - x**2 + y**2)**2 + z**2 - one/100
    #print(p)

    print("polynomial { %d" % (p.degree,), end="")
    for key in p.keys:
        print(",")
        x, y, z = 0, 0, 0
        for item in key:
            if item[0]=='x':
                x = item[1]
            if item[0]=='y':
                y = item[1]
            if item[0]=='z':
                z = item[1]
        print("    xyz(%s,%s,%s):%s" % (x, y, z, p.cs[key]), end="")

    print("\n    sturm\n}")
    


# ----------------------------------------------------------------------
#
#


class Formal(series.Series):
    """
        Formal power series with coefficients in a polynomial ring 
        with infinitely many variables.
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

    # https://golem.ph.utexas.edu/category/2018/01/more_secrets_of_the_associahed.html

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

