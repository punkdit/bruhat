#!/usr/bin/env python3



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
        assert len(coefs) == len(cs)
        self.cs = cs

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
        items.sort()
        terms = []
        for (k, v) in items:
            ss = []
            for name, exp in k:
                assert exp>0
                if len(name)>1:
                    name += " "
                if exp == 1:
                    ss.append(name)
                else:
                    ss.append("%s^%s"%(name, exp))
            s = ''.join(ss)
            assert v != 0
            if v==1:
                terms.append(s)
            else:
                terms.append( "%s*%s" % (v, s))
        s = "+ ".join(terms)
        s = s.replace("+ -", "- ")
        return s

    __repr__ = __str__



def test():

    zero = Poly({})
    one = Poly({():1})
    x = Poly("x")
    y = Poly("y")
    xx = Poly({(("x", 2),) : 1})
    xy = Poly({(("x", 1), ("y", 1)) : 1})

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

    print(a, b)
    print(a*b)



if __name__ == "__main__":

    test()


