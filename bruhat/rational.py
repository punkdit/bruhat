#!/usr/bin/env python3

"""
Find power series of rational functions
"""

from bruhat.util import factorial, cross
from bruhat.element import Q
from bruhat.poly import Poly


def choose(m, n):
    assert m>=n>=0
    value = factorial(m) // (factorial(n) * factorial(m-n))
    assert value > 0
    return value


def lstr(items):
    if isinstance(items, list):
        s = "[%s]"%(', '.join(str(c) for c in items))
    elif isinstance(items, dict):
        s = "{%s}"%(', '.join("%s:%s"%(k, str(v)) for k,v in items.items()))
    else:
        s = str(items)
    return s


def tpl_range(tpl): # inclusive !
    if len(tpl) == 0:
        yield () # ?
        return 
    if len(tpl) == 1:
        k = tpl[0]
        for i in range(k+1): # inclusive !
            yield (i,)
        return
    rest = tpl[1:]
    k = tpl[0]
    for idxs in tpl_range(rest):
        for i in range(k+1):
            yield (i,)+idxs


def tpl_sub(idxs, jdxs):
    return tuple(i-j for (i,j) in zip(idxs, jdxs))

def tpl_choose(idxs, jdxs):
    assert len(idxs) == len(jdxs)
    r = 1
    for i in range(len(idxs)):
        r *= choose(idxs[i], jdxs[i])
    return r

def tpl_factorial(idxs):
    r = 1
    for i in idxs:
        r *= factorial(i)
    return r


assert list(tpl_range((1,2,1))) == [
    (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 2, 0), 
    (1, 2, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1), (0, 2, 1), (1, 2, 1)]


class Rational(object):
    def __init__(self, base, p, q, vs=None):
        # f(x) == p(x) / q(x)
        # q(x)*f(x) == p(x)
        p = Poly.promote(p, base)
        q = Poly.promote(q, base)
    
        zero, one = base.zero, base.one
        
        if vs is None:
            vs = p.get_vars() + q.get_vars()
            vs = list(set(vs))
            vs.sort()
        n = len(vs)
        assert n>0
    
        #print("divpoly(%s, %s, %s)" % (p, q, vs))
    
        fs = {}
        vzero = {v:zero for v in vs}
    
        top = p(**vzero)
        bot = q(**vzero)
        #print("top=%s, bot=%s" % (top, bot))
        f0 = top/bot
        #print("f0:", lstr(f0))
        fs[(0,)*n] = f0

        self.base = base
        self.p = p
        self.q = q
        self.fs = fs
        self.vzero = vzero
        self.vs = vs
        self.bot = bot
    
    def pump(self):
        base = self.base
        p = self.p
        q = self.q
        fs = self.fs
        vzero = self.vzero
        vs = self.vs
        bot = self.bot
        n = len(vs)
    
        deg = 0
        for tpl in fs.keys():
            deg = max(deg, sum(tpl))
        deg += 1
        for idxs in tpl_range((deg,)*(n-1)):
            total = sum(idxs) 
            if total > deg:
                continue
            idxs += (deg-total,)
            #print(idxs)
            top = p.partial(**{vs[i]:idxs[i] for i in range(n)})(**vzero) # eye's hurting?
            #print(top)
            for jdxs in tpl_range(idxs):
                if jdxs == idxs:
                    continue
                kdxs = tpl_sub(idxs, jdxs)
                q0 = q.partial(**{vs[i]:kdxs[i] for i in range(n)})(**vzero)
                top -= tpl_choose(idxs, jdxs) * fs[jdxs] * q0
            #assert top.is_const()
            f0 = top / bot
            assert fs.get(idxs) is None or fs.get(idxs) == f0
            fs[idxs] = f0
            #print(idxs, f0)
            #print("\t", idxs)
    
    def __getitem__(self, idx):
        fs = self.fs
        vs = self.vs
        if type(idx) is int:
            idx = idx,
        assert len(idx) == len(vs)
        for i in idx:
            assert 0<=i
        #print("__getitem__", idx)
        while fs.get(idx) is None:
            self.pump()
        #coeffs = {idx : fs[idx] // tpl_factorial(idx) for idx in fs.keys()}
        #return coeffs
        val = fs[idx] // tpl_factorial(idx)
        return val





def test():

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)
    xx = Poly({(("x", 2),) : 1}, ring)
    xy = Poly({(("x", 1), ("y", 1)) : 1}, ring)

    f = Rational(ring, one, (1-x)**2)
    cs = [f[i] for i in range(5)]
    assert cs == [1, 2, 3, 4, 5]

    f = Rational(ring, 1-x**2, (1-x)**6)
    cs = [f[i] for i in range(5)]
    assert cs == [1, 6, 20, 50, 105]

    f = Rational(ring, 1-x*y, (1-x)**3 * (1-y)**3)
    cs = [[f[i,j] for i in range(4)] for j in range(4)]
    assert cs == [[1, 3, 6, 10], [3, 8, 15, 24], [6, 15, 27, 42], [10, 24, 42, 64]]


    J = Poly("J", ring)
    L = Poly("L", ring)
    Li = Poly("Li", ring)
    vs = [J, L, Li]
    f = Rational(ring, one, (1-J*L)*(1-J*Li), "J L Li".split())
    for i in range(4):
      for j in range(4):
        for k in range(4):
            print(f[i, j, k], end=' ')
        print()
      print()



if __name__ == "__main__":

    test()

