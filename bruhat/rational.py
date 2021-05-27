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


def divpoly(base, p, q, count=5):
    # f(x) == p(x) / q(x)
    # q(x)*f(x) == p(x)
    p = Poly.promote(p, base)
    q = Poly.promote(q, base)

    if q.is_const():
        q = q.get_const()
        result = (1/q)*p
        if result.is_const():
            result = result.get_const() # demote
        return result # <------------------- return

    zero, one = base.zero, base.one
    
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

    for deg in range(count):
        for idxs in tpl_range((count,)*(n-1)):
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

    coeffs = {idx : fs[idx] // tpl_factorial(idx) for idx in fs.keys()}
    return coeffs





def test():

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)
    xx = Poly({(("x", 2),) : 1}, ring)
    xy = Poly({(("x", 1), ("y", 1)) : 1}, ring)

    cs = divpoly(ring, one, (1-x)**2)
    print(lstr(cs))

    #assert cs == [1, 2, 3, 4, 5]
    cs = divpoly(ring, 1-x**2, (1-x)**6)
    print(lstr(cs))
    #assert cs == [1, 6, 20, 50, 105]

    cs = divpoly(ring, 1-x*y, (1-x)**3 * (1-y)**3)
    print(lstr(cs))


if __name__ == "__main__":

    test()

