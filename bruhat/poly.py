#!/usr/bin/env python3

"""
Multivariate commutative polynomials.
"""

from functools import reduce
import operator
from random import randint, shuffle, seed

from bruhat.argv import argv
from bruhat.util import cross, factorial
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


#def tpl_add_old(tj, tk):
#    tj = dict(tj)
#    tk = dict(tk)
#    for k,v in tk.items():
#        tj[k] = tj.get(k, 0) + v
#        assert tj[k] > 0
#    tj = list(tj.items())
#    tj.sort()
#    tj = tuple(tj)
#    return tj
                

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

def tpl_ge(left, right):
    for (k, l, r) in tpl_zip(left, right):
        if l<r:
            return False
    return True

def tpl_compare(left, right):
    l = sum([l[1] for l in left], 0)
    r = sum([r[1] for r in right], 0)
    if l > r:
        return True
    if l == r and left>right: # ?
        return True
    return False


#def tpl_div(top, bot): # not used
#    top = dict(top)
#    bot = dict(bot)
#    result = {}
#    for k,v in bot.items():
#        w = top.get(k, 0)
#        if w < v:
#            return None
#        if v < w:
#            result[k] = w-v
#    result = list(result.items())
#    result.sort()
#    return result


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
            if tpl_compare(key, head):
                head = key
            dkey = dict(key)
            assert len(dkey) == len(key)
            deg = 0
            for v in dkey.values():
                assert v, "got some zero exponents: %s"%repr(dkey)
                deg += v
            degree = max(degree, deg)
            coefs[key] = value
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

    def get_const(self):
        return self.cs.get((), self.ring.zero)

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
        m0 = P.head
        for m1 in list(Q.cs.keys()):
            if tpl_ge(m1, m0):
                break
        else:
            return P.get_zero(), Q

        m2 = tpl_sub(m1, m0)
        #print("reduce: m1-m0 = ", m1, m0, m2)
        ring = P.ring
        m2 = Poly({m2:ring.one}, ring)
        R = (Q.get(m1)//P.get(m0)) * m2
        S = Q - R * P
        assert S.degree <= Q.degree
        return R, S

    def reduce(P, Q):
        "return R, S such that Q=R*P+S "
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
        return Poly(cs, ring)

    def __pow__(self, n):
        ring = self.ring
        if n==0:
            return Poly({():ring.one}, ring)
        assert n>0
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

    def substitute(self, ns):
        s = self.python_str()
        assert type(ns) is dict, repr(ns)
        p = eval(s, ns)
        return p

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
        print("reduce_many", len(seen))
        for Q in polys:
            print(Q, "reduce", P)
            Q1, P1 = Q.reduce(P)
            print("\t", Q1, "||", P1)
            P = P1
        if P in seen:
            break
        print("\t", P)
        seen.add(P)
        #assert len(seen) < 5
    return P


def grobner(polys):
    "Build a grobner basis from polys."

    assert polys
    ring = polys[0].ring
    polys = [p for p in polys if p]

    basis = set(polys)
    pairs = set()
    for i in range(len(polys)):
      for j in range(i+1, len(polys)):
        pairs.add((polys[i], polys[j]))

    while pairs:
        print("grobner", len(pairs))
        P, Q = pairs.pop()
        S = P.critical_pair(Q)
        if S==0:
            continue

        #write("(")
        S = reduce_many(basis, S)
        #write(")")
        if S==0:
            continue

        S = S.normalized()

        for P in basis:
            pairs.add((P, S))
        basis.add(S)
        #write('.')

    return basis






def test():

    global ring

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

    assert (x*y + y) / y == (x + 1)
    assert (x**2 - y**2) / (x+y) == (x-y)

    # -------------------------------------
    # 

    for trial in range(20):
        p = Poly.random(list('xyz'), ring)
        q = Poly.random(list('xyz'), ring)
        #print p, "||", q
        r, s = p.reduce(q)

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

    p = x**2*y-x**2
    q = x*y**2-y**2
    polys = [p, q]
    basis = grobner(polys)
    #print basis

    for i in range(10):
        polys = [Poly.random(list('xyz'), ring) for _ in range(2)]
        basis = grobner(polys)
        #print len(basis)

    polys = [x+y+z, x*y+x*z+y*z, x*y*z]
    basis = grobner(polys)
    #for p in basis:
    #    print(p)
    #print(basis)


    print("OK")


def test_hilbert():
    # ---------------------------------------------------------------
    # Here we work out the coefficients of a Hilbert polynomial
    # given by the rational function top/bot.
    # These coefficients give the dimension of irreps of SL(3).

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
        print(val, end=" ")
      print()



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

    fn = argv.next()
    if fn is None:
        test()
    else:
        fn = eval(fn)
        fn()


