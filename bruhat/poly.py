#!/usr/bin/env python3

from bruhat.util import cross
from bruhat.theta import divisors


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
                


class Poly(object):

    def __init__(self, cs):
        if isinstance(cs, str):
            cs = {((cs,1),):1}
        coefs = {}
        for key, value in cs.items():
            if value == 0:
                continue
            key = list(key)
            key.sort()
            key = tuple(key)
            dkey = dict(key)
            assert len(dkey) == len(key)
            for v in dkey.values():
                assert v
            coefs[key] = value
        #assert len(coefs) == len(cs), (coefs, cs)
        self.cs = coefs

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
        items.sort()
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
                    ss.append("%s^%s"%(name, exp))
            s = ' '.join(ss)
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
        s = s.replace("+ -", "- ")
        s = s.replace(" ^", "^")
        if s and s[-1] == " ":
            s = s[:-1]
        return s

    __repr__ = __str__


def get_a(i):
    if i==0:
        return 0 # a_0 == 0
    if i==1:
        return 1 # a_1 == 1
    return Poly("a_%d"%i)


def get_b(i):
    if i==0:
        return 0 # b_0 == 0
    if i==1:
        return 1 # b_1 == 1
    return Poly("b_%d"%i)


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


def dirichlet():
    zero = Poly({})
    n = 1
    while 1:

        v = zero
        for j in divisors(n):
            k = n//j
            v += get_a(j) * get_b(k)
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

    a = zero
    b = zero
    for i in range(5):
        a += Poly("a_%d"%i)
        b += Poly("b_%d"%i)

    for i in range(2, 6):
        for n in range(1, 6):
            vs = pow_b(n, i)
            #print("pow_b(%d, %d) = %s" % (n, i, vs))
            print(vs, end= ", ")
        print()
    print()


    #gen = mul()
    gen = compose()
    #gen = dirichlet()
    for i in range(1, 6):
        p = gen.__next__()
        print("compose(%d) = %s"%(i, p))
        print()



if __name__ == "__main__":

    test()


