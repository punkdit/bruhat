#!/usr/bin/env python3

"""
Univariate formal power series over a ring: R[[x]].
"""

from bruhat.element import Q
from bruhat.util import factorial, cross, factorize
from bruhat.theta import divisors


class Series(object):
    prec = 6

    def __init__(self, ring):
        self.ring = ring
        self._cache = {}

    def __str__(self):
        return "Series(%s, ...)" % (
            ', '.join(str(value) for value in [self[i] for i in range(self.prec)]))

    def getitem(self, idx):
        return self.ring.zero

    def __getitem__(self, idx):
        assert isinstance(idx, int)
        assert idx >= 0
        value = self._cache.get(idx)
        if value is None:
            value = self.getitem(idx)
            self._cache[idx] = value
        return value

    @classmethod
    def promote(cls, item, ring):
        if isinstance(item, Series):
            return item
        return RMul(item, One(ring))

    def __add__(a, b):
        b = Series.promote(b, a.ring)
        return Add(a, b)

    def __neg__(a):
        return RMul(-self.ring.one, a)

    def __sub__(a, b):
        b = Series.promote(b, a.ring)
        return Add(a, RMul(-a.ring.one, b))

    def __rmul__(a, r):
        return RMul(r, a)

    def __mul__(a, b):
        b = Series.promote(b, a.ring)
        return Mul(a, b)

    def __pow__(a, n):
        assert isinstance(n, int)
        if n==0:
            return One(a.ring) # what about 0**0 ?
        if n==1:
            return a
        if n==2:
            return a*a
        return Pow(a, n)

    def __call__(a, b):
        b = Series.promote(b, a.ring)
        return Compose(a, b)

    def eq(a, b, idx=10):
        b = Series.promote(b, a.ring)
        for i in range(idx):
            if a[i] != b[i]:
                return False
        return True

    def dirichlet(a, b):
        b = Series.promote(b, a.ring)
        return Dirichlet(a, b)

    def is_multiplicative(D):
        # just test a few terms here !
        if D[0] != 0:
            return False
        if D[1] != 1:
            return False
        for (m, n) in [(2,3), (2,5), (3,5), (2,7), (3,7), (5,7)]:
            if D[m]*D[n] != D[m*n]:
                return False
        return True


class FSeries(Series):

    def __init__(self, ring, *items):
        Series.__init__(self, ring)
        self.items = list(items)

    def getitem(self, idx):
        items = self.items
        if idx >= len(items):
            return self.ring.zero
        return items[idx]


class RMul(Series):
    def __init__(self, r, child):
        assert isinstance(child, Series)
        Series.__init__(self, child.ring)
        self.r = r
        self.child = child

    def getitem(self, idx):
        return self.r * self.child[idx]


class Binop(Series):
    def __init__(self, a, b):
        assert isinstance(a, Series)
        assert isinstance(b, Series)
        assert a.ring is b.ring
        Series.__init__(self, a.ring)
        self.a = a
        self.b = b


class Add(Binop):
    def getitem(self, idx):
        return self.a[idx] + self.b[idx]


class Sub(Binop):
    def getitem(self, idx):
        return self.a[idx] - self.b[idx]


class Mul(Binop):
    def getitem(self, idx):
        a, b = self.a, self.b
        value = self.ring.zero
        for j in range(idx+1):
            value += a[j] * b[idx-j]
        return value


class Pow(Series):
    def __init__(self, child, n):
        assert isinstance(child, Series)
        Series.__init__(self, child.ring)
        self.n = n
        assert n>2
        self.child = child

    def getitem(self, idx):
        idxs = list(range(idx+1))
        n = self.n
        child = self.child
        ring = self.ring
        value = ring.zero
        for jdxs in cross([idxs]*(n-1)):
            total = 0
            walue = ring.one
            for j in jdxs:
                walue *= child[j]
                total += j
            if total <= idx:
                walue *= child[idx-total]
                value += walue
        return value


class Compose(Binop):
    def __init__(self, a, b):
        assert b[0] == a.ring.zero, b[0]
        Binop.__init__(self, a, b)

    def getitem(self, idx):
        a, b = self.a, self.b
        value = self.ring.zero
        for i in range(idx+1):
            value += a[i] * ((b**i)[idx])
        return value


class Dirichlet(Binop):
    def __init__(self, a, b):
        assert b[0] == a.ring.zero, b[0] # do we care?
        Binop.__init__(self, a, b)

    def __getitem__(self, idx):
        a = self.ring.zero
        if idx == 0:
            return a
        for d in divisors(idx):
            a += self.a[d] * self.b[idx//d]
        return a


# ----------------------------------------

class Zero(Series):
    pass

class One(Series):
    def getitem(self, idx):
        ring = self.ring
        if idx == 0:
            return ring.one
        return ring.zero

class X(Series):
    def getitem(self, idx):
        ring = self.ring
        if idx == 1:
            return ring.one
        return ring.zero


class Exp(Series):
    def getitem(self, idx):
        ring = self.ring
        return ring.one // factorial(idx)


# ----------------------------------------

#class F(Series):
#    def getitem(self, idx):

def divisors(n):
    assert 1<=n
    i = 1
    while i<=n:
        if n%i == 0:
            yield i
        i += 1
    

class Eisenstein(Series):
    def __init__(self, ring, j, const=0):
        assert j>=0
        Series.__init__(self, ring)
        self.j = j
        self.const = const

    def getitem(self, n):
        assert 0<=n
        if n==0:
            return self.const

        ring = self.ring
        c = 0
        #print("factorize:", n, "=", list(factorize(n)))
        for p in divisors(n):
            #print('\n', p, p**self.j)
            c += p**(self.j - 1)
        return ring.promote(c)
        

class EulerFactor(Series):
    def __init__(self, M, p):
        Series.__init__(self, M.ring)
        assert M.is_multiplicative()
        self.M = M
        self.p = p

    def getitem(self, i):
        return self.M[self.p**i]

        
def eisenstein():
    Series.prec = 6

    ring = Q

    zero = Zero(ring)
    one = One(ring)
    x = X(ring)

    c = FSeries(ring, ring.one // 240)

    const = ring.one // 240
    G4 = Eisenstein(ring, 4, const)
    #print(G4)

    const = -ring.one // 504
    G6 = Eisenstein(ring, 6, const)
    #print(G6)

    D = (240 * G4) ** 3 - (504 * G6) ** 2
    D = (ring.one // ring.promote(1728)) * D
    print("D =", D)

    #return

    # hmm... fail...
    #M = G4
    #M = Eisenstein(ring, 4, 1) * D

    if 1:
        M = Eisenstein(ring, 4, ring.one)*D
        print("M =", M)
        for i,j in [(2,3), (2,5), (3, 5)]:
            lhs = M[i]*M[j]
            rhs = M[i*j]
            print(lhs, rhs, rhs/lhs)
        assert M.is_multiplicative()


    G4 = Eisenstein(ring, 4)
    G6 = Eisenstein(ring, 6)
    assert G4.is_multiplicative()
    assert G6.is_multiplicative()
    assert D.is_multiplicative()

    #return

    G4_2 = EulerFactor(G4, 2)
    #print(G4_2)

    p = 2
    D2 = EulerFactor(D, p)
    #print("D2 =", D2)

    D3 = EulerFactor(D, 3)
    #print("D3 =", D3)

    def check(F, p):
        for n in range(5):
            lhs = F[n+2]
            rhs = -(p**11)*F[n] + F[1]*F[n+1] 
            assert lhs == rhs

    #check(D2, 2)
    #check(D3, 3)

    def find_coeffs(F):
        from gelim import array, solve
        f = [int(str(F[i])) for i in range(4)]
        A = array([[f[0], f[1]], [f[1], f[2]]])
        rhs = [f[2], f[3]]
        v = solve(A, rhs)
        assert v is not None
        a, b = v
        a, b = int(str(a)), int(str(b))
        return a, b

    for p in [2, 3, 5]:
        F = EulerFactor(G4, p)
        #print(F)
        a, b = find_coeffs(F)
        #print(p, a, b)
        assert (a, b) == (-p**3, F[1]) # X_0(3) ?

    for p in [2, 3, 5]:
        F = EulerFactor(G6, p)
        #print(F)
        a, b = find_coeffs(F)
        #print(p, a, b)
        assert (a, b) == (-p**5, F[1]) # X_0(5) ?

    assert find_coeffs(D2) == (-(2**11), D2[1]) # X_0(11) ?
    assert find_coeffs(D3) == (-(3**11), D3[1])


def main():
    ring = Q

    zero = Zero(ring)
    one = One(ring)
    x = X(ring)
    exp = Exp(ring)

    assert ((x+1)**2).eq(x**2 + 2*x + 1)

    a = 2*x+1
    assert (a**3).eq( a*a*a )
    assert (a**4).eq( a*a*a*a )


    #print(exp(zero))
    #print(exp(x+x))
    #print(exp(x**2))
    #print(exp(-x))

    f = exp - 1
    #f = exp - 1 + x**3 # works too...
    print("f =", f)

    # -----------------------------
    # From: http://math.ucr.edu/home/baez/trimble/trimble_lie_operad.pdf
    # Page 8.
    # find compositional inverse

    assert f[0] == 0
    assert f[1] == 1
    
    def delta(h):
        return h(f) - h

    g = zero
    d = x
    for n in range(5):
        s = (-1)**n * d
        d = delta(d)
        g += s

    print("g =", g)

    print("f(g) = ", f(g))
    print("g(f) = ", g(f))

    # ------------------------------------------------------------------------------------------
    # Verifying https://golem.ph.utexas.edu/category/2018/01/more_secrets_of_the_associahed.html

    C = f
    D = g

    c1, c2, c3, c4 = C[2], C[3], C[4], C[5]
    d1, d2, d3, d4 = D[2], D[3], D[4], D[5]
    
    assert(d1 == -c1)
    assert(d2 == -c2 + 2*c1**2)
    assert(d3 == -c3 + 5*c2*c1 - 5*c1**3)
    assert(d4 == -c4 + 6*c3*c1 + 3*c2**2 - 21*c2*c1**2 + 14*c1**4)
    
    




if __name__ == "__main__":

    eisenstein()

    print("OK\n")

        


