#!/usr/bin/env python

"""

Vector's, Algebra's, Coalgebra's, Hopf Algebra's, 
over the integers.

see also: hopf.py for previous ideas...

TODO: use other monoidal products (categories) ?


reference:
[HaGuKi]
_Algebras,_Rings_and_Modules 
Lie_Algebras_and_Hopf_Algebras (2009)
Michiel Hazewinkel
Nadiya Gubareni
V. V. Kirichenko

"""

import numpy
from random import randint
from functools import reduce
from operator import matmul

from bruhat.util import cross
from bruhat.gset import Perm, Group
#from bruhat.element import Q, Z, Integer, Fraction
from bruhat import element
from bruhat.argv import argv

def mktuple(k):
    if type(k) is tuple:
        return k
    return (k,)


def unify(a_shape, b_shape):
    if a_shape is None:
        return b_shape
    if b_shape is None:
        return a_shape
    if a_shape == b_shape:
        return a_shape
    assert 0, (a_shape, b_shape)

ring = element.Z
ring.tp = element.Integer
scalar = ring.promote

#ring = element.Q
#ring.tp = element.Fraction


class Vector:
    """
        A dict with int values (the scalar coeffs),
        and keys are n-tuples (tensors) of hashable data (stuff).
        <shape> is the length of the tuples, use shape=0 for 
        scalar values.
    """
    def __init__(self, coeffs={}, shape=None):
        coeffs = {mktuple(k):scalar(v)
            for (k,v) in coeffs.items() if scalar(v)}
        keys = list(coeffs.keys())
        keys.sort()
        self.keys = keys
        for k in keys:
            if shape is None:
                shape = len(k)
            assert len(k)==shape
        self._hash = hash(tuple((k,coeffs[k]) for k in keys))
        self.coeffs = coeffs
        self.shape = shape

    @classmethod
    def scalar(cls, item):
        item = scalar(item)
        if type(item) is ring.tp:
            return Vector({ () : item }, 0)
        assert 0, (type(item), scalar)

    @classmethod
    def promote(cls, item):
        if isinstance(item, Vector):
            return item
        item = scalar(item)
        if type(item) is ring.tp:
            return Vector({ () : item }, 0)
        assert 0, (type(item), scalar)

    def __str__(self):
        return "Vector(%s)"%(self.coeffs,)
    __repr__ = __str__

    def __len__(self):
        return len(self.coeffs)

    def __iter__(self):
        return iter(self.coeffs)

    def __getitem__(self, k):
        k = mktuple(k)
        return self.coeffs.get(k, 0)

    def __eq__(self, other):
        #if other==0:
        #    return self.coeffs == {}
        other = Vector.promote(other)
        return self.coeffs==other.coeffs

    def __hash__(self):
        return self._hash

    def __add__(self, other):
        other = Vector.promote(other)
        assert isinstance(other, Vector), (self, other)
        shape = unify(self.shape, other.shape)
        coeffs = dict(self.coeffs)
        for (k,v) in other.coeffs.items():
            coeffs[k] = coeffs.get(k,0) + v
        return Vector(coeffs, shape)

    def __sub__(self, other):
        other = Vector.promote(other)
        assert isinstance(other, Vector), other
        shape = unify(self.shape, other.shape)
        coeffs = dict(self.coeffs)
        for (k,v) in other.coeffs.items():
            coeffs[k] = coeffs.get(k,0) - v
        return Vector(coeffs, shape)

    def __neg__(self):
        return (-1)*self

    def __rmul__(self, r):
        r = scalar(r)
        if r==0:
            return Vector({})
        return Vector({k:r*v for (k,v) in self.coeffs.items()})

    def __matmul__(self, other):
        if type(other) is ring.tp:
            return other*self
        assert isinstance(other, Vector), type(other)
        coeffs = {}
        for u in self:
          for v in other:
            uv = u+v # tuple concat
            coeffs[uv] = coeffs.get(uv, 0) + self[u]*other[v]
        return Vector(coeffs)

    def __rmatmul__(self, other):
        assert isinstance(other, ring.tp)
        return other*self

    def swap(self):
        if not len(self):
            return self
        coeffs = {}
        for k in self:
            assert len(k) == 2
            j = k[1],k[0]
            coeffs[j] = self[k]
        return Vector(coeffs)
        


class Lin:
    """
        Multi-linear (aka tensor) operator.
        Lazy.
        Only defined on "basis tensors", anything
        else we extend linearly...
    """
    def __init__(self, shape):
        assert type(shape) is tuple, type(shape)
        assert len(shape) == 2
        self.shape = shape # (m,n): m <--- n

    def apply(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        assert 0, "abstract base class"

    def __call__(self, v): 
        "linearly extend apply "
        v = Vector.promote(v)
        assert isinstance(v, Vector)
        m, n = self.shape
        unify(v.shape, n)
        u = Vector()
        for k in v:
            assert len(k) == n
            u = u + v[k]*self.apply(*[Vector({ki:1}) for ki in k])
        return u

    def __mul__(self, other):
        assert isinstance(other, Lin)
        return ComposeLin(self, other)

    def __matmul__(self, other):
        return TensorLin(self, other)

    def test(self, vs):
        # assert some multi-linearity...
        m, n = self.shape
        for arg in cross([vs]*n):
            u = self.apply(*arg)
            for i in range(n):
                barg = list(arg)
                r = randint(-5,5)
                barg[i] = r*barg[i]
                assert self.apply(*barg) == r*u
        # XX test additivity.. v+w .. 


class ILin(Lin):
    """ _identity Lin
    """
    def __init__(self, m=1):
        Lin.__init__(self, (m,m))

    def apply(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        #print("__call__", list(vecs))
        if n==1:
            v = vecs[0]
        else:
            v = reduce(matmul, vecs)
        assert isinstance(v, Vector)
        return v


class PermLin(Lin):
    """ _permutation Lin
    """
    def __init__(self, perm):
        m = len(perm)
        items = list(set(perm))
        items.sort()
        assert items == list(range(m)), perm
        Lin.__init__(self, (m,m))
        self.perm = perm

    def apply(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        #print("call", list(vecs))
        perm = self.perm
        vecs = [vecs[i] for i in perm]
        v = reduce(matmul, vecs)
        assert isinstance(v, Vector)
        return v


class OpLin(Lin):
    def __init__(self, op, shape):
        Lin.__init__(self, shape)
        self.op = op

    def apply(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        return self.op(*vecs)


class TensorLin(Lin):
    def __init__(self, aop, bop):
        m = aop.shape[0]+bop.shape[0]
        n = aop.shape[1]+bop.shape[1]
        Lin.__init__(self, (m,n))
        self.ops = (aop, bop)

    def apply(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        aop, bop = self.ops
        a = vecs[:aop.shape[1]]
        b = vecs[aop.shape[1]:]
        a, b = aop.apply(*a) , bop.apply(*b)
        #assert isinstance(a, (int, Vector)), (aop, a)
        #assert isinstance(b, (int, Vector)), (bop, b)
        a = Vector.promote(a)
        b = Vector.promote(b)
        assert isinstance(a, (ring.tp, Vector)), (aop, a)
        assert isinstance(b, (ring.tp, Vector)), (bop, b)
        #print("TensorLin.apply", a, b)
        return a@b


class ComposeLin(Lin):
    def __init__(self, aop, bop):
        assert aop.shape[1] == bop.shape[0]
        Lin.__init__(self, (aop.shape[0], bop.shape[1]))
        self.ops = (aop, bop)

    def apply(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        aop, bop = self.ops
        v = bop.apply(*vecs)
        if type(v) is int:
            v = v*aop(1)
        else:
            v = aop(v)
        return v



class Algebra:
    def __init__(self, unit, mul):
        assert isinstance(unit, Lin), type(unit)
        assert isinstance(mul, Lin), type(mul)
        assert unit.shape == (1,0)
        assert mul.shape == (1,2)
        self.unit = unit
        self.mul = mul

    def test(self, vs, comm=False):
        unit, mul = self.unit, self.mul
        unit.test(vs)
        mul.test(vs)
        # _unital
        ident = unit(1)
        m = mul.apply
        for u in vs:
            assert m(ident, u) == u
            assert m(u, ident) == u
        # assoc
        for u in vs:
          for v in vs:
           for w in vs:
            lhs = m(m(u,v), w)
            rhs = m(u, m(v, w))
            assert lhs == rhs
        if comm:
            for u in vs:
             for v in vs:
                assert m(u,v) == m(v,u)


class Coalgebra:
    def __init__(self, counit, comul):
        assert isinstance(counit, Lin)
        assert isinstance(comul, Lin)
        assert counit.shape == (0,1)
        assert comul.shape == (2,1)
        self.counit = counit
        self.comul = comul

    def test(self, vs, cocomm=False):
        counit, comul = self.counit, self.comul
        counit.test(vs)
        comul.test(vs)
        I = ILin()

        for v in vs:
            m = comul(v)
            # _counital
            u = (I@counit)(m)
            assert u==v
            u = (counit@I)(m)
            assert u==v
            # coassoc
            mmi = (comul@I)(m)
            imm = (I@comul)(m)
            assert mmi==imm

            #if cocomm:

class Hopf:
    def __init__(self, unit, mul, counit, comul, antipode):
        assert isinstance(unit, Lin), type(unit)
        assert isinstance(mul, Lin), type(mul)
        assert unit.shape == (1,0)
        assert mul.shape == (1,2)
        self.unit = unit
        self.mul = mul

        assert isinstance(counit, Lin)
        assert isinstance(comul, Lin)
        assert counit.shape == (0,1)
        assert comul.shape == (2,1)
        self.counit = counit
        self.comul = comul

        assert isinstance(antipode, Lin)
        assert antipode.shape == (1,1)
        self.antipode = antipode

    def test(self, vs):
        unit, mul = self.unit, self.mul
        counit, comul = self.counit, self.comul
        antipode = self.antipode

        # bialgebra
        for u in vs:
          for v in vs:
            assert counit( mul(u@v) ) == counit(u)@counit(v)
        assert counit(unit(1)) == 1
        assert comul(unit(1)) == unit(1)@unit(1)

        I = ILin()
        S = PermLin([0,2,1,3])
        lin = (mul@mul)*S*(comul@comul)
        for u in vs:
          for v in vs:
            lhs = comul(mul(u@v))
            rhs = lin(u@v)
            assert lhs == rhs

        # hopf algebra
        #for v in vs:
        #    print()
        #    print("v =", v)
        #    print("antipode(v) =", antipode(v))
        antipode.test(vs)

        lhs = mul * (antipode @ I) * comul
        mid = unit * counit
        rhs = mul * (I @ antipode) * comul
        for u in vs:
            #print()
            #print("u =", u)
            #print("comul(u) =", comul(u))
            #print("counit(u) =", counit(u))
            m = mid(u)
            l = lhs(u)
            r = rhs(u)
            #print("mid:", m)
            #print("lhs:", l)
            #print("rhs:", r)
            assert l==m
            assert r==m


        #print("Hopf.test: OK")

class Frobenius:
    def __init__(self, unit, mul, counit, comul):
        assert isinstance(unit, Lin), type(unit)
        assert isinstance(mul, Lin), type(mul)
        assert unit.shape == (1,0)
        assert mul.shape == (1,2)
        self.unit = unit
        self.mul = mul

        assert isinstance(counit, Lin)
        assert isinstance(comul, Lin)
        assert counit.shape == (0,1)
        assert comul.shape == (2,1)
        self.counit = counit
        self.comul = comul

    def test(self, vecs):
        unit, mul = self.unit, self.mul
        counit, comul = self.counit, self.comul

        I = ILin()
        lhs = (mul @ I) * (I @ comul)
        mhs = comul * mul
        rhs = (I @ mul) * (comul @ I)
        for u in vecs:
         for v in vecs:
            l = lhs(u@v)
            m = mhs(u@v)
            r = rhs(u@v)
            assert l==m
            assert r==m
        
    @classmethod
    def function(cls, keys):
        # --------------------------------------
        # function Frobenius algebra
    
        zero = ring.zero
        one = ring.one
    
        keys = [mktuple(key) for key in keys]
        basis = [Vector({key:1}) for key in keys]
        unit = OpLin(lambda : 
            Vector({key:1 for key in keys}), (1,0))
        mul = OpLin(lambda a,b: 
            Vector({key:a[key]*b[key] for key in keys}), (1,2))
    
        #algebra = Algebra(unit, mul)
        #algebra.test(basis, True)
        #a, b, c = basis[:3]
        #algebra.test([a+b, a+b+c, c], True)
        
        def counit(v):
            s = zero
            for key in keys:
                s += v[key]
            return Vector.scalar(s)
        counit = OpLin(counit, (0,1))
        comul = OpLin(lambda v: 
            Vector({key+key:v[key] for key in keys}), (2,1))
        #coalgebra = Coalgebra(counit, comul)
        #coalgebra.test(basis, True)
        #coalgebra.test([a, a+2*b, a+b+c, c], True)
        frobenius = Frobenius(unit, mul, counit, comul)
        frobenius.basis = basis
        return frobenius




def all_bits(n, m):
    # XXX this is n+m choose m
    nm = n+m
    for idxs in numpy.ndindex((2,)*nm):
        if sum(idxs) == m:
            yield idxs

    
def str_shuffle(s, t):
    #print("str_shuffle(%r, %r)"%(s, t))

    for bits in all_bits(len(s), len(t)):
        #print(bits)
        idx = 0
        jdx = 0
        word = ""
        for bit in bits:
            if bit==0:
                word += s[idx]
                idx += 1
            else:
                word += t[jdx]
                jdx += 1
        assert len(word) == len(s)+len(t)
        yield word


def test_hopf():
    global ring, scalar
    ring = element.Z
    ring.tp = element.Integer
    scalar = ring.promote

    # ----------------------------------------------

    w = Vector({'abc':1})
    v = Vector({'ac':1})
    #assert str(w) == "Vector({'abc': 1})"
    #assert str(w) == "Vector({('abc',): 1})"
    assert (w+w) == 2*w
    assert 0*w == 0
    assert w-v != 0
    assert w-w == 0

    assert w+v == v+w
    assert -1*(w+v) == -w-v

    vector = lambda s : Vector({s:1})
    a = vector("a")
    b = vector("b")
    c = vector("c")
    assert a@b == Vector( {("a","b"):1} )
    assert a@(b+2*c) == Vector( {("a","b"):1, ("a","c"):2} )

    assert Vector.promote(27) == Vector({():27}, 0)
    assert 3*Vector.promote(2) == 6*Vector.promote(1)

    # ----------------------------------------------
    # the word Algebra
    unit = OpLin(lambda : Vector({'':1}), (1,0))
    assert unit(2) == 2*unit(1)

    def mul(u, v):
        coeffs = {}
        for k, in u:
            for j, in v:
                kj = k+j
                coeffs[kj] = u[k] * v[j] + coeffs.get(kj, 0)
        return Vector(coeffs)
    mul = OpLin(mul, (1,2))

    vs = [vector(s) for s in 'ab a c ac bc abc abbc'.split()]
    ab,a,c,ac,bc,abc,abbc = vs

    assert mul.apply(ab+a, c-bc) == -abbc + ac
    assert mul.apply(unit(1), ab+a) == ab+a
                
    algebra = Algebra(unit, mul)

    for u in [c,ac,abc]:
      for w in [ab,a,c]:
        vs.append(u+w)
    algebra.test(vs)

    for u in vs:
      for v in vs:
        lhs, rhs = mul.apply(u,v) , mul(u@v)
        assert lhs == rhs, (lhs, rhs)

    # ----------------------------------------------
    # copy

    counit = OpLin(lambda v : sum(v[k] for k in v), (0,1))
    counit.test(vs)

    def copy(v):
        coeffs = {}
        for k in v:
            coeffs[k+k] = v[k]
        return Vector(coeffs)
    copy = OpLin(copy, (2,1))
    copy.test(vs)

    coalgebra = Coalgebra(counit, copy)
    coalgebra.test(vs)

    # ----------------------------------------------
    # the shuffle Hopf algebra: [HaGuKi] Ex. 3.4.6

    unit = OpLin(lambda : Vector({'':1}), (1,0))

    assert list(str_shuffle("abc", "ef")) == [
        'abcef', 'abecf', 'abefc', 'aebcf', 'aebfc', 'aefbc',
        'eabcf', 'eabfc', 'eafbc', 'efabc']

    def mul(u, v):
        coeffs = {}
        for k, in u:
          for j, in v:
            for kj in str_shuffle(k,j):
              coeffs[kj] = coeffs.get(kj,0) + u[k]*v[j]
        return Vector(coeffs)
    mul = OpLin(mul, (1,2))

    cd = vector('cd')
    assert mul.apply(ab, cd) == Vector(
        {'abcd': 1, 'acbd': 1, 'acdb': 1, 
        'cabd': 1, 'cadb': 1, 'cdab': 1})

    aaa = vector('aaa')
    assert mul.apply(a,aaa) == 4*vector('aaaa')
    assert mul(a@aaa) == mul.apply(a,aaa) 

    algebra = Algebra(unit, mul)
    algebra.test([a,ab,a-cd], True)

    def comul(u):
        coeffs = {}
        for k, in u:
            for i in range(len(k)+1):
                key = (k[:i],k[i:])
                coeffs[key] = coeffs.get(key,0) + u[k]
        return Vector(coeffs)
    comul = OpLin(comul, (2,1))

    def counit(u):
        r = 0
        for k, in u:
            if len(k)==0:
                r += u[k]
        return r
    counit = OpLin(counit, (0,1))

    assert counit(a) == 0
    assert counit(2*unit(1)+a) == 2

    assert comul(abc) == Vector(
        {('', 'abc'): 1, 
        ('a', 'bc'): 1, 
        ('ab', 'c'): 1, 
        ('abc', ''): 1})

    coalgebra = Coalgebra(counit, comul)
    coalgebra.test(vs)
    

    def antipode(v):
        coeffs = {}
        for k, in v:
            k1 = ''.join(reversed(k))
            r = (-1)**len(k1)
            coeffs[k1,] = v[k]*r
        return Vector(coeffs)



    antipode = OpLin(antipode, (1,1))

    hopf = Hopf(
        unit, mul,
        counit, comul,
        antipode)
    hopf.test(vs)
    

def test_frobenius():
    global ring, scalar
    ring = element.Q
    ring.tp = element.Fraction
    scalar = ring.promote

    # --------------------------------------
    # function Frobenius algebra

    keys = [(i,) for i in range(3)]
    frobenius = Frobenius.function(keys)

    a, b, c = frobenius.basis[:3]
    vecs = [a, a+2*b, a+b-7*c, c]
    algebra = Algebra(frobenius.unit, frobenius.mul)
    algebra.test(vecs, True)
    coalgebra = Coalgebra(frobenius.counit, frobenius.comul)
    coalgebra.test(vecs, True)
    frobenius.test(vecs)

    # --------------------------------------
    # Group Hopf algebra

    G = Group.symmetric(3)
    print(G)
    #for g in G:
    #    print(g)

    i = G.identity
    unit = OpLin(lambda : Vector({(i,):1}), (1,0))
    def mul(a, b):
        coeffs = {}
        for g, in a:
          for h, in b:
            gh = g*h
            coeffs[gh] = coeffs.get(gh, 0) + a[g]*b[h]
        return Vector(coeffs)
    mul = OpLin(mul, (1,2))
    algebra = Algebra(unit, mul)

    frobenius = Frobenius.function(G)
    a, b, c = frobenius.basis[:3]
    vecs = [a, a+2*b, a+b-7*c, c]
    frobenius.test(vecs)

    algebra.test(vecs)



def test():
    #test_hopf()
    test_frobenius()


if __name__ == "__main__":

    from time import time

    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "test"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()



