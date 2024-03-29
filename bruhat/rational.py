#!/usr/bin/env python3

"""
Find power series of rational functions
"""

from functools import lru_cache

from bruhat.util import factorial, cross
from bruhat.argv import argv
if argv.fast:
    from bruhat._element import Q
else:
    from bruhat.element import Q
from bruhat.poly import Poly


@lru_cache(maxsize=None)
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
    "Rational function p/q where p and q are Poly's over a base ring."
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
        #assert n>0
    
        #print("Rational(%s, %s, %s)" % (p, q, vs))
    
        fs = {}
        vzero = tuple((v,zero) for v in vs)
    
        top = p.substitute(vzero)
        bot = q.substitute(vzero)
        if bot == 0:
            fs = None
        else:
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

    def promote(self, other):
        if isinstance(other, Rational):
            return other
        #if isinstance(other, Poly):
        #    return other
        other = Poly.promote(other, self.base)
        other = Rational(self.base, other, self.base.one)
        return other

    def __eq__(lhs, rhs):
        return lhs.p*rhs.q == rhs.p*lhs.q
    
    def __ne__(lhs, rhs):
        return lhs.p*rhs.q != rhs.p*lhs.q

    def __str__(self):
        return "(%s)/(%s)" % (self.p, self.q)
    __repr__ = __str__

    def __add__(self, other):
        other = self.promote(other)
        assert self.base == other.base
        atop, abot = self.p, self.q
        btop, bbot = other.p, other.q
        top = atop * bbot + btop * abot
        bot = abot * bbot
        return Rational(self.base, top, bot)

    __radd__ = __add__ # commutative

    def __sub__(self, other):
        other = self.promote(other)
        assert self.base == other.base
        atop, abot = self.p, self.q
        btop, bbot = other.p, other.q
        top = atop * bbot - btop * abot
        bot = abot * bbot
        return Rational(self.base, top, bot)

    def __rsub__(self, other):
        # other - self
        other = self.promote(other)
        return self.__neg__().__add__(other)

    def __neg__(self):
        return Rational(self.base, -self.p, self.q)

    def __mul__(self, other):
        other = self.promote(other)
        assert self.base == other.base
        return Rational(self.base, self.p*other.p, self.q*other.q)

    __rmul__ = __mul__ # commutative

    def __truediv__(self, other):
        other = self.promote(other)
        assert self.base == other.base
        return Rational(self.base, self.p*other.q, self.q*other.p)

    def __rtruediv__(self, other):
        # other / self
        other = self.promote(other)
        assert self.base == other.base
        return Rational(self.base, other.p*self.q, other.q*self.p)

    def __call__(self, **kw):
        p, q = self.p, self.q
        kw = tuple(kw.items())
        p = p.substitute(kw)
        q = q.substitute(kw)
        return Rational(self.base, p, q)
    
    def _pump(self, verbose=False):
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
            top = p.partial(**{vs[i]:idxs[i] for i in range(n)})
            if verbose:
                print("_pump: idxs=%s"%(idxs,), end=" ")
            top = top.substitute(vzero)
            for jdxs in tpl_range(idxs):
                if jdxs == idxs or fs[jdxs]==0:
                    continue
                if verbose:
                    print("*", end=" ")
                kdxs = tpl_sub(idxs, jdxs)
                q0 = q.partial(**{vs[i]:kdxs[i] for i in range(n)})
                #print("\t", fs[jdxs])
                q0 = q0.substitute(vzero)
                top -= tpl_choose(idxs, jdxs) * fs[jdxs] * q0
            if verbose:
                print()
            #assert top.is_const()
            f0 = top / bot
            assert fs.get(idxs) is None # or fs.get(idxs) == f0
            fs[idxs] = f0
            #print(idxs, f0)
            #print("\t", idxs)
    
    def __getitem__(self, idx):
        fs = self.fs
        if fs is None:
            raise Exception("rational expression %s has no derivatives at zero."%self)
        vs = self.vs
        if type(idx) is int:
            idx = idx,
        assert len(idx) == len(vs)
        for i in idx:
            assert 0<=i
        #print("__getitem__", idx)
        while fs.get(idx) is None:
            self._pump()
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

    # course generating function for SL(2) irreps
    f = Rational(ring, one, (1-x)**2)
    cs = [f[i] for i in range(5)]
    assert cs == [1, 2, 3, 4, 5]

    f = Rational(ring, 1-x**2, (1-x)**6)
    cs = [f[i] for i in range(5)]
    assert cs == [1, 6, 20, 50, 105]

    # course generating function for SL(3) irreps
    f = Rational(ring, 1-x*y, (1-x)**3 * (1-y)**3)
    cs = [[f[i,j] for i in range(4)] for j in range(4)]
    assert cs == [[1, 3, 6, 10], [3, 8, 15, 24], [6, 15, 27, 42], [10, 24, 42, 64]]

    # fine generating function for SL(2) irreps
    J = Poly("J", ring)
    L = Poly("L", ring)
    Li = Poly("Li", ring)
    vs = [J, L, Li]
    f = Rational(ring, one, (1-J*L)*(1-J*Li), "J L Li".split())
    assert f(L=one, Li=one) == Rational(ring, one, (1-J)**2)

    if 0:
        promote = lambda p : Rational(ring, p, one)
    
        r_one = Rational(ring, one, one)
        r_J = promote(J)
        r_L = promote(L)
        r_Li = r_one / r_L
        f = r_one / (1-r_J*r_L)*(1-r_J*r_Li)
        print(f)
        for i in range(4):
          for j in range(4):
            print(f[i, j], end=' ')
          print()
    
        return

    for i in range(4):
      for j in range(4):
        for k in range(4):
            print(f[i, j, k], end=' ', flush=True)
        print()
      print()

    # fine generating function for SL(3) irreps
    J = Poly("J", ring)
    K = Poly("K", ring)
    L = Poly("L", ring)
    Li = Poly("Li", ring)
    M = Poly("M", ring)
    Mi = Poly("Mi", ring)
    vs = [J, K, L, Li, M, Mi]
    f = Rational(ring, 
        1-J*K, 
        (1-J*L)*(1-J*M*Li)*(1-J*Mi)*(1-K*Li)*(1-K*L*Mi)*(1-K*M), 
        "J K L Li M Mi".split())
    assert f(L=one, Li=one, M=one, Mi=one) \
        == Rational(ring, 1-J*K, (1-J)**3 * (1-K)**3)

    if 0:
        f._pump(verbose=True)
        f._pump(verbose=True)
        f._pump(verbose=True)
    
        return

    def sample(ji, ki, li, mi):
        N = 3
        if li>=0 and mi>=0:
            val = sum(f[ji, ki, l+li, l,    m+mi, m   ] for l in range(N) for m in range(N))
        elif li>=0 and mi<0:
            val = sum(f[ji, ki, l+li, l,    m,    m-mi] for l in range(N) for m in range(N))
        elif li<0 and mi<0:
            val = sum(f[ji, ki, l,    l-li, m,    m-mi] for l in range(N) for m in range(N))
        elif li<0 and mi>=0:
            val = sum(f[ji, ki, l,    l-li, m+mi, m   ] for l in range(N) for m in range(N))
        return val

    assert sample(0, 0, 0, 0) == 1
    for (li, mi) in [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1), (0, 0), (0, 1),
        (1, -1), (1, 0), (1, 1),
    ]:
        print(sample(1, 0, li, mi))

    for li in range(-2, 3):
      for mi in range(-2, 3):
        print(sample(1, 0, li, mi), end=" ", flush=True)
      print()

    #assert sum(f[1, 1, i, i, j, j] for i in range(4) for j in range(4)) == 2
    #assert sum(f[1, 1, i+1, i, j+1, j] for i in range(4) for j in range(4)) == 1
    #assert sum(f[1, 1, i+1, i, j, j] for i in range(4) for j in range(4)) == 0




if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("test()")
    else:
        test()


