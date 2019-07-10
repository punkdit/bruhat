#!/usr/bin/env python3

from functools import reduce
import operator

from bruhat.argv import argv
from bruhat.util import cross
from bruhat.theta import divisors
from bruhat import series


def is_scalar(x):
    return isinstance(x, int)


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


class Poly(object):

    def __init__(self, cs):
        if isinstance(cs, str):
            cs = {((cs,1),):1}
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

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, idx):
        key = self.keys[idx]
        val = self.cs[key]
        return Poly({key : val})

    def terms(self):
        cs = self.cs
        for key in self.keys:
            val = self.cs[key]
            yield Poly({key : val})

    @classmethod
    def promote(cls, item):
        if isinstance(item, Poly):
            return item
        assert is_scalar(item)
        return Poly( {() : item} )

    def __eq__(self, other):
        other = self.promote(other)
        return self.cs == other.cs

    def __ne__(self, other):
        other = self.promote(other)
        return self.cs != other.cs

    def __hash__(self):
        cs = list(self.cs.items())
        cs.sort(key = lambda k_v:k_v[0])
        cs = tuple(cs)
        return hash(cs)

    def __add__(self, other):
        other = self.promote(other)
        cs = dict(self.cs)
        for key, value in list(other.cs.items()):
            cs[key] = cs.get(key, 0) + value
        return Poly(cs)
    __radd__ = __add__

    def __sub__(self, other):
        other = self.promote(other)
        cs = dict(self.cs)
        for key, value in list(other.cs.items()):
            cs[key] = cs.get(key, 0) - value
        return Poly(cs)

    def __neg__(self):
        cs = {}
        for key, value in list(self.cs.items()):
            cs[key] = -value
        return Poly(cs)

    def __rmul__(self, r):
        assert is_scalar(r)
        cs = {}
        for key, value in list(self.cs.items()):
            cs[key] = r * value
        return Poly(cs)

    def __mul__(self, other):
        if is_scalar(other):
            return self.__rmul__(other)
        cs = {}
        for k1, v1 in list(self.cs.items()):
          for k2, v2 in list(other.cs.items()):
            k = tpl_add(k1, k2)
            cs[k] = cs.get(k, 0) + v1*v2
        return Poly(cs)

    def __pow__(self, n):
        if n==0:
            return Poly({():1})
        assert n>0
        p = self
        for i in range(n-1):
            p = self*p
        return p

    def __str__(self):
        items = list(self.cs.items())
        if not items:
            return "0"
        items.sort(key = (lambda kv:tpl_degree(kv[0])))
        #print(items)
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
                        ss.append("%s^{%s}"%(name, exp))
                    else:
                        ss.append("%s^%s"%(name, exp))
            s = '*'.join(ss)
            if s:
                if v==1:
                    terms.append(s)
                else:
                    terms.append("%s*%s" % (v, s))
            else:
                if v==1:
                    terms.append("1")
                else:
                    terms.append("%s" % v)
        s = " + ".join(terms)
        s = s.replace("-1*", "-")
        s = s.replace("+ -", "- ")
        s = s.replace(" ^", "^")
        if s and s[-1] == " ":
            s = s[:-1]
        return s
    __repr__ = __str__

    def python_str(self):
        s = self.__str__()
        s = s.replace("^", "**")
        s = s.replace("{", "")
        s = s.replace("}", "")
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
        return Poly("a_%d"%i)
    else:
        return Poly("a_{%d}"%i)


def get_b(i, b_1=False):
    if i==0:
        return 0 # b_0 == 0
    if i==1 and not b_1:
        return 1 # b_1 == 1
    if i < 10:
        return Poly("b_%d"%i)
    else:
        return Poly("b_{%d}"%i)


def mul():
    zero = Poly({})
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
    zero = Poly({})
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
    zero = Poly({})
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
    zero = Poly({})
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

    zero = Poly({})
    one = Poly({():1})
    x = Poly("x")
    y = Poly("y")
    xx = Poly({(("x", 2),) : 1})
    xy = Poly({(("x", 1), ("y", 1)) : 1})

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
        a += Poly("a_%d"%i)
        b += Poly("b_%d"%i)

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

        s = "%s = %s" % (b_i, rhs)
        s = s.replace("*", " ")
        s = s.replace("=", "&=")
        print(r"    %s \\"%s)
        #print()


def main():

    zero = Poly({})
    one = Poly({():1})

    ring = object()
    ring.zero = zero
    ring.one = one

    


if __name__ == "__main__":

    test()


