#!/usr/bin/env python3

"""
Commutative _polynomials.

See also: poly.py 
"""

# Code converted to python3 from cpoly.py by 2to3


import sys, os
from fractions import Fraction
from random import randint, shuffle, seed


#from bruhat.action import Group
from bruhat.gelim import fstr
from bruhat.util import write


def tpl_add(a, b):
    "tuple add"
    c = tuple(x+y for (x,y) in zip(a, b))
    return c

def tpl_sub(a, b):
    "tuple sub"
    c = tuple(x-y for (x,y) in zip(a, b))
    return c

def tpl_ge(a, b):
    "component-wise tuple comparison"
    for i, j in zip(a, b):
        if i<j:
            return False
    return True

def tpl_union(a, b):
    "component-wise maximum"
    c = tuple(max(ai, bi) for (ai, bi) in zip(a, b))
    return c

def tpl_compare(a, b):
    if sum(a) > sum(b):
        return True
    if sum(a) == sum(b) and a>b:
        return True
    return False


def is_scalar(x):
    return isinstance(x, (int, Fraction))


class Poly(object):
    """
        Commutative _polynomial in <rank> many variables x_1,...,x_<rank> .
    """

    def __init__(self, cs, rank=None, names=None):
        coefs = {}
        degree = 0
        head = None
        for key, value in list(cs.items()):
            assert rank is None or rank == len(key)
            rank = len(key)
            assert is_scalar(value)
            if value == 0:
                continue
            if head is None or tpl_compare(key, head):
                head = key
            coefs[key] = Fraction(value)
            degree = max(degree, sum(key))
        if head is None:
            head = (0,)*rank
        self.cs = coefs
        self.rank = rank
        self.degree = degree
        self.head = head
        self.names = names
        if names is not None:
            assert len(names) == rank
        assert sum(head) == degree, str(self)

    def flatstr(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort(reverse=True)
        ss = [''.join(str(k) for k in key) for key in keys]
        s = ["%s:%s" % (ss[i], cs[k]) for i, k in enumerate(keys)]
        return "{" + ", ".join(s) + "}"

    @classmethod
    def identity(cls, rank, names=None):
        key = (0,)*rank
        return cls({key : 1}, rank)

    @classmethod
    def zero(cls, rank, names=None):
        return cls({}, rank)

    def get_zero(self):
        return Poly({}, self.rank, self.names)

    def get_identity(self):
        key = (0,)*self.rank
        return Poly({key : 1}, self.rank, self.names)

    def promote(self, x):
        if is_scalar(x):
            return Poly({(0,)*self.rank : x}, self.rank, self.names)
        if isinstance(x, Poly):
            return x
        raise TypeError(repr(x))

    def __bool__(self):
        return len(self.cs)>0

    def __eq__(self, other):
        other = self.promote(other)
        assert self.rank == other.rank
        return self.cs == other.cs

    def __ne__(self, other):
        other = self.promote(other)
        assert self.rank == other.rank
        return self.cs != other.cs

    def __hash__(self):
        cs = list(self.cs.items())
        cs.sort(key = lambda k_v:k_v[0])
        cs = tuple(cs)
        return hash(cs)

    def __getitem__(self, key):
        assert len(key) == self.rank
        assert isinstance(key, tuple)
        return self.cs.get(key, 0)

    def __add__(self, other):
        if is_scalar(other):
            return self + Poly({(0,)*self.rank : other}, self.rank, self.names)
        assert self.rank == other.rank
        cs = dict(self.cs)
        for key, value in list(other.cs.items()):
            cs[key] = cs.get(key, 0) + value
        return Poly(cs, self.rank, self.names)

    def __sub__(self, other):
        if is_scalar(other):
            return self - Poly({(0,)*self.rank : other}, self.rank, self.names)
        assert self.rank == other.rank
        cs = dict(self.cs)
        for key, value in list(other.cs.items()):
            cs[key] = cs.get(key, 0) - value
        return Poly(cs, self.rank, self.names)

    def __neg__(self):
        cs = {}
        for key, value in list(self.cs.items()):
            cs[key] = -value
        return Poly(cs, self.rank, self.names)

    def __rmul__(self, r):
        cs = {}
        for key, value in list(self.cs.items()):
            cs[key] = r * value
        return Poly(cs, self.rank, self.names)

    def __mul__(self, other):
        if is_scalar(other):
            return self * Poly({(0,)*self.rank : other}, self.rank, self.names)
        assert self.rank == other.rank
        cs = {}
        for k1, v1 in list(self.cs.items()):
          for k2, v2 in list(other.cs.items()):
            k = tpl_add(k1, k2)
            cs[k] = cs.get(k, 0) + v1*v2
        return Poly(cs, self.rank, self.names)

    def __pow__(self, n):
        if n==0:
            return Poly.identity(self.rank, self.names)
        assert n>0
        p = self
        for i in range(n-1):
            p = self*p
        return p

    def __str__(self):
        rank = self.rank
#        if rank==0:
#            return 
        if self.names is not None:
            names = self.names
        elif rank>1:
            names = ["x_%d"%(i+1) for i in range(rank)]
        else:
            names = ["x"]
        items = []
        cs = list(self.cs.items())
        #cs.sort(key = lambda (k,v):list(reversed(k)))
        cs.sort(key = lambda k_v1:(-sum(k_v1[0], 0), list(reversed(k_v1[0]))))
        for k, val in cs:
            if val==0:
                continue
            item = []
            for i, exp in enumerate(k):
                if exp == 0:
                    pass
                elif exp == 1:
                    item.append("%s" % (names[i],))
                else:
                    s_exp = str(exp)
                    s_exp = s_exp if len(s_exp)==1 else "{"+s_exp+"}"
                    item.append("%s^%s" % (names[i], s_exp))
            item = ' '.join(item)
            if not item:
                item = str(val)
            elif val==1:
                pass
            elif val==-1:
                item = "-"+item
            else:
                item = "%s*%s" % (fstr(val), item)
            items.append(item)
        s = ' + '.join(items)
        s = s.replace("+ -", "- ")
        if not s:
            s = "0"
        return s

    __repr__ = __str__

    def py_str(self):
        s = self.__str__()
        s = s.replace(" + ", "+")
        s = s.replace(" ", "*")
        s = s.replace("{", "")
        s = s.replace("}", "")
        s = s.replace("^", "**")
        return s

    def get_namespace(self):
        s = self.py_str()
        rank = self.rank
        namespace = {}
        for i in range(rank):
            exp = [0]*rank
            exp[i] = 1
            pi = Poly({tuple(exp) : 1}, rank, self.names)
            namespace[pi.py_str()] = pi
        return namespace

    def substitute(self, vs):
        # vs : map name to new Poly
        s = self.py_str()
        ns = self.get_namespace()
        for k, v in vs.items():
            assert type(k) is str
            assert isinstance(v, Poly)
            #assert v.names == self.names # too strict. commented for now.
            #s = s.replace(k, "(%s)"%(v.py_str()))
            assert k in ns, "unknown name %r"%k
            ns[k] = v
        p = eval(s, ns, ns)
        return p
    
    def transform(self, A):
        import numpy
        rank = self.rank
        names = self.names
        A = numpy.array(A)
        assert A.shape == (rank, rank)
        A = A.astype(object)
        vec = []
        for i in range(rank):
            exp = [0]*rank
            exp[i] = 1
            p = Poly({tuple(exp) : 1}, rank, names)
            vec.append(p)
        vec = numpy.array(vec)
        ns = {}
        for i, p in enumerate(vec):
            n = names[i]
            n = n.replace("{", "") # ARGH!
            n = n.replace("}", "") # ARGH!
            ns[n] = numpy.dot(A[i], vec)
            #print(n, "-->", numpy.dot(A[i], vec))
        return self.substitute(ns)

    @classmethod
    def random(cls, rank, degree=3, terms=3):
        cs = {}
        for i in range(terms):
            key = [0]*rank
            remain = degree
            for j in range(rank):
                key[j] = randint(0, remain)
                remain -= key[j]
            shuffle(key)
            key = tuple(key)
            value = randint(-2, 2)
            cs[key] = value
        P = cls(cs, rank)
        return P

    def is_symmetric(self):
        n = self.rank
        ns = list(range(n))
        cs = self.cs

        #keys = list(cs.keys())
        for k0, v0 in list(cs.items()):

            # cycle
            k1 = tuple(k0[(i+1)%n] for i in ns)
            if cs.get(k1) != v0:
                return False

            if n <= 1:
                continue
            
            # swap
            k1 = tuple(k0[i if i>1 else (1-i)] for i in ns)
            if cs.get(k1) != v0:
                return False

        return True

    def reduce1(P, Q):
        "return R, S such that Q=R*P+S "
        m0 = P.head
        for m1 in list(Q.cs.keys()):
            if tpl_ge(m1, m0):
                break
        else:
            return P.get_zero(), Q

        #print "reduce1(%s, %s)" % (P, Q)
        #print "reduce: P.head = m0 = ", m0
        #print "reduce: m1 = ", m1
        m2 = tpl_sub(m1, m0)
        #print "reduce: m1-m0 = ", m2
        m2 = Poly({m2:1}, P.rank, P.names)
        #print "m2:", m2
        R = (Q[m1]/P[m0]) * m2
        S = Q - R * P
        #if S.degree > Q.degree:
        #    print "%s reduce %s" % (P, Q)
        #    print "(%s).head = %s" % (P, P.head,)
        #    print "\t (%s) = (%s)*(%s) + (%s)" % (Q, R, P, S)
        assert S.degree <= Q.degree
        return R, S

    def reduce(P, Q):
        "return R, S such that Q=R*P+S "
        #print "reduce(%s, %s)" % (P, Q)
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
            #print R1, S
        return R, S

    def __div__(self, other):
        if is_scalar(other):
            other = Fraction(1, other)
            other = Poly({(0,)*self.rank : other}, self.rank, self.names)
            return other * self
        R, S = other.reduce(self)
        return R
    __truediv__ = __div__

    def swap(self, i, j=None):
        rank = self.rank
        assert 0<=i<rank
        if j is None:
            j = i+1
        assert 0<=j<rank
        if i==j:
            return self
        cs = {}
        for key, value in list(self.cs.items()):
            key = list(key)
            key[i], key[j] = key[j], key[i]
            key = tuple(key)
            cs[key] = value
        return Poly(cs, rank, self.names)

    def basis(self, i):
        key = [0 for j in range(self.rank)]
        key[i] = 1
        return Poly({tuple(key) : 1}, self.rank, self.names)

    def delta(P, i):
        "divided difference operator"
        assert i+1<P.rank
        Q = P - P.swap(i, i+1)
        #print "delta", P, "=", Q
        R = P.basis(i) - P.basis(i+1)
        #print "divide by", R
        S, T = R.reduce(Q)
        #print "equals:", S
        assert not T
        return S

    def normalized(P):
        "divide by leading coefficient"
        head = P.head
        c = P[head]
        if c != 1:
            P = P / c
        return P

    def critical_pair(P, Q):
        P = P.normalized()
        Q = Q.normalized()
        m = tpl_union(P.head, Q.head)
        mi = tpl_sub(m, P.head)
        mj = tpl_sub(m, Q.head)
        mi = Poly({mi:1}, P.rank, P.names)
        mj = Poly({mj:1}, P.rank, P.names)
        S = mi * P - mj * Q
        return S


#class Ideal(object):
#    def __init__(self, ps):


def reduce_many(polys, P):

    seen = set([P])
    while 1:
        for Q in polys:
            _, P = Q.reduce(P)
        if P in seen:
            break
        seen.add(P)
        #write(str(len(seen)))

    return P


def grobner(polys):
    "Build a grobner basis from polys."

    polys = [p for p in polys if p]
    
    basis = set(polys)
    pairs = set()
    for i in range(len(polys)):
      for j in range(i+1, len(polys)):
        pairs.add((polys[i], polys[j]))

    while pairs:
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




def main():

    seed(0) # deterministic

    I = Poly({(0,) : 1})
    x = Poly({(1,) : 1})

    assert Poly.identity(1) == I
    assert x.get_identity() == I

    assert (I + x*x) * (I+x*x) == (I+x**2)**2
    assert (I + x*x) * (I+x*x) != (I+x)**2
    assert (I+x*x) * (I-x*x) == I-x**4
    assert (I-x**4).is_symmetric()
    assert x*x == x**2

    assert (I+x)**4 == I + 4*x + 6*x**2 + 4*x**3 + x**4

    # ----------------

    I = Poly.identity(3)
    zero = Poly.zero(3)
    x1 = Poly({(1,0,0) : 1})
    x2 = Poly({(0,1,0) : 1})
    x3 = Poly({(0,0,1) : 1})

    p = x1+x2+x3
    q = p.substitute({"x_2": x3})
    assert q == x1 + 2*x3, q

    #print((x1+x2+x3)**3)
    #print(x1*x2*x3)

    assert ((x1+x2+x3)**2 + I).is_symmetric()
    assert (x1*x2*x3).is_symmetric()

    assert not ((x1+x2+x3)**3 + x1).is_symmetric()
    assert not (x1*x2*x3 + x1*x2).is_symmetric()

    half = Fraction(1, 2)
    P = half*x1**3 + x2**2
    assert P.head == (3, 0, 0), P.head
    assert P.normalized() == x1**3 + 2*x2**2

    _, Q = P.reduce(P)
    assert Q == zero

    P = x1**2 * x2 + x1*x2 - 17
    Q = x1**2 * x2**2 + x1*x2**2 - 3
    _, R = P.reduce(Q)
    assert R == 17*x2 - 3

    assert (x1*x2 + x2) / x2 == (x1 + 1)
    assert (x1**2 - x2**2) / (x1+x2) == (x1-x2)

    # -------------------------------------
    # 

    for trial in range(20):
        P = Poly.random(3)
        Q = Poly.random(3)
        #print P, "||", Q
        R, S = P.reduce(Q)

    # -------------------------------------
    # 

    P = x1**2*x2-x1**2
    Q = x1*x2**2-x2**2
    assert P.critical_pair(Q) == x1*x2**2 - x1**2*x2
    assert P.critical_pair(Q) == -Q.critical_pair(P)
    assert P.critical_pair(P) == 0
    assert Q.critical_pair(Q) == 0


    # -------------------------------------
    # 

    P = x1**2*x2-x1**2
    Q = x1*x2**2-x2**2
    polys = [P, Q]
    basis = grobner(polys)
    #print basis

    for i in range(10):
        polys = [Poly.random(3) for _ in range(2)]
        basis = grobner(polys)
        #print len(basis)

    # -------------------------------------
    # divided difference operator:

    P = x1**2 * x2

    assert (x1**2 * x2).delta(0) == x1*x2
    assert (x1*x2).delta(0)      == 0
    assert (x1**2 * x2).delta(1) == x1**2

    assert (x1**2).delta(0)      == x1+x2
    assert (x1**2).delta(1)      == 0
    assert (x1+x2).delta(0)      == 0
    assert (x1+x2).delta(1)      == 1
    assert x1.delta(0)           == 1
    assert x1.delta(1)           == 0

    P = (3*x1**2*x3 + x2**2 + x2*x3)
    Q = (2*x2**2*x3 + x1**2 + x2*x3**2)

    rank = 3
    for trial in range(100):

        P = Poly.random(rank)
        Q = Poly.random(rank)

        for i in range(rank-2):

            lhs = P.delta(i).delta(i+1).delta(i)
            rhs = P.delta(i+1).delta(i).delta(i+1)
            assert lhs == rhs


        for i in range(rank-1):

            assert P.delta(i).delta(i) == 0

            lhs = (P*Q).delta(i)
            rhs = P.delta(i)*Q + P.swap(i, i+1) * Q.delta(i)
            assert lhs == rhs

    print("OK")



    


if __name__ == "__main__":

    main()


