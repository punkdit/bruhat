#!/usr/bin/env python3

from functools import reduce
import operator

from bruhat.argv import argv
from bruhat.util import cross, factorial
from bruhat.theta import divisors
from bruhat.element import Fraction, Q
from bruhat import series


#def is_scalar(x):
#    return isinstance(x, int)


def tpl_add(tj, tk):
    tj = dict(tj)
    tk = dict(tk)
    for k,v in tk.items():
        tj[k] = tj.get(k, 0) + v
        assert tj[k] > 0
    tj = list(tj.items())
    tj.sort()
    tj = tuple(tj)
    return tj
                

def tpl_div(top, bot): # not used
    top = dict(top)
    bot = dict(bot)
    result = {}
    for k,v in bot.items():
        w = top.get(k, 0)
        if w < v:
            return None
        if v < w:
            result[k] = w-v
    result = list(result.items())
    result.sort()
    return result


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
        for key, value in cs.items():
            if value == 0:
                continue
            key = list(key)
            key.sort()
            key = tuple(key)
            dkey = dict(key)
            assert len(dkey) == len(key)
            deg = 0
            for v in dkey.values():
                assert v
                deg += v
            degree = max(degree, deg)
            coefs[key] = value
            keys.append(key)
        #assert len(coefs) == len(cs), (coefs, cs)
        keys.sort()
        self.keys = keys
        self.cs = coefs
        self.degree = degree
        self.ring = ring

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, idx):
        key = self.keys[idx]
        val = self.cs[key]
        return Poly({key : val}, self.ring)

    def terms(self):
        cs = self.cs
        for key in self.keys:
            val = self.cs[key]
            yield Poly({key : val}, self.ring)

    @classmethod
    def promote(cls, item, ring):
        if isinstance(item, Poly):
            return item
        #assert is_scalar(item)
        item = ring.promote(item)
        return Poly({() : item}, ring)

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
        #s = s.replace("-1"+MUL, "-")
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
        return s

#    def subs(self, vals): # argh... too hard
#        print("subs", vals, "->", self)
#        if not vals:
#            return self
#        if len(self)>1:
#            # recurse
#            return reduce(operator.add, [p.subs(vals) for p in self.terms()])
#        print("subs", vals, "->", self.cs)
#        cs = list(self.cs.items())
#        assert len(cs) == 1
#        key, value = cs[0]
#        k = tpl_div(key, bot)
#        if k is not None:
#            key = k




        


def get_a(i, a_1=False):
    if i==0:
        return 0 # a_0 == 0
    if i==1 and not a_1:
        return 1 # a_1 == 1
    if i < 10:
        return Poly("a_%d"%i, ring)
    else:
        return Poly("a_{%d}"%i, ring)


def get_b(i, b_1=False):
    if i==0:
        return 0 # b_0 == 0
    if i==1 and not b_1:
        return 1 # b_1 == 1
    if i < 10:
        return Poly("b_%d"%i, ring)
    else:
        return Poly("b_{%d}"%i, ring)


def mul():
    zero = Poly({}, ring)
    #yield zero

    n = 1
    while 1:

        #print("mul: n=%s"%n)
        v = zero
        for j in range(n+1):
            k = n-j
            v += get_a(j) * get_b(k)
            #print(j, k, v)
        yield v

        n += 1


def pow_b(n, i):
    assert n>0
    #print("pow_b(%d, %d)"%(n, i))
    if n==1:
        return get_b(i)
    zero = Poly({}, ring)
    v = zero
    items = list(range(1, i+1))
    for idxs in cross([items]*(n-1)):
        total = sum(idxs)
        if total >= i:
            continue
        #idxs.append(i - total)
        last = i-total
        #cs = dict((("b_%d"%idx, 1), 1) for idx in idxs)
        p = get_b(last)
        for idx in idxs:
            p *= get_b(idx)
        v += p
    return v


def dirichlet(linear_term=False):
    zero = Poly({}, ring)
    n = 1
    while 1:

        v = zero
        for j in divisors(n):
            k = n//j
            v += get_a(j, linear_term) * get_b(k, linear_term)
        yield v

        n += 1


def compose():
    #print("compose")
    zero = Poly({}, ring)
    #yield zero

    i = 1
    while 1:
        n = 1
        q = zero
        while 1:
            v = pow_b(n, i)
            if v==0:
                break
            p = get_a(n)
            #print(p*v, end= ", ")
            q += p*v
            n += 1
        #print()
        i += 1
        yield q
    #print()


def test():
    global ring

    ring = Q

    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
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

#    for i in range(2, 6):
#        for n in range(1, 6):
#            vs = pow_b(n, i)
#            #print("pow_b(%d, %d) = %s" % (n, i, vs))
#            print(vs, end= ", ")
#        print()
#    print()

#    gen = dirichlet(True)
#    for i in range(1, 13):
#        s = str(gen.__next__())
#        s = s.replace("*", " ")
#        print(r"(%s) x^{%d} \\"%(s, i))
#    return

    # Lagrange inversion...

    name = argv.next()
    assert name and name in "mul compose dirichlet".split()
    gen = eval(name)

    idx = 2
    soln = {}
    items = gen()
    N = argv.get("N", 5)
    for i in range(1, N+1):
        p = items.__next__()
        if p.degree == 0:
            continue
        #print("p =", p)
        b_i = get_b(idx, False)
        a_i = get_a(idx, False)
        soln[a_i.python_str()] = a_i
        rhs = b_i - p
        #print("eval", rhs, soln)
        rhs = eval(rhs.python_str(), soln)
        soln[b_i.python_str()] = rhs
        idx += 1

        print(r"    %s &= %s \\" % (b_i.str(), rhs.str()))



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



def main():

    zero = Poly({}, Q)
    one = Poly({():Q.one}, Q)

    ring = type("ARing", (object,), {})
    ring.zero = zero
    ring.one = one


    f = FormalExp("a", ring)
    g = FormalExp("b", ring)

    fg = f*g
    for i in range(5):
        print(fg[i])
    print()
    
    g = FormalExp("b", ring, {"b_0":0})
    print("g =", g)
    fg = f(g)
    for i in range(5):
        print(fg[i])
    print()
    


if __name__ == "__main__":

    if argv.test:
        test()
    else:
        main()


