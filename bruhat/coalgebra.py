#!/usr/bin/env python

"""
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
from bruhat.argv import argv

def mktuple(k):
    if type(k) is tuple:
        return k
    return (k,)


class Vector:
    """
        A dict with int values (the scalar coeffs),
        and keys are n-tuples (tensors) of hashable data (stuff)
    """
    def __init__(self, coeffs={}):
        coeffs = {mktuple(k):int(v)
            for (k,v) in coeffs.items() if int(v)}
        keys = list(coeffs.keys())
        keys.sort()
        self.keys = keys
        shape = None
        for k in keys:
            if shape is None:
                shape = len(k)
            assert len(k)==shape
        self._hash = hash(tuple((k,coeffs[k]) for k in keys))
        self.coeffs = coeffs

    def __str__(self):
        return "Vector(%s)"%(self.coeffs,)
    __repr__ = __str__

    def __len__(self):
        return len(self.coeffs)

    def __iter__(self):
        return iter(self.coeffs)

    def __getitem__(self, k):
        k = mktuple(k)
        return self.coeffs[k]

    def __eq__(self, other):
        if other==0:
            return self.coeffs == {}
        return self.coeffs==other.coeffs

    def __hash__(self):
        return self._hash

    def __add__(self, other):
        coeffs = dict(self.coeffs)
        for (k,v) in other.coeffs.items():
            coeffs[k] = coeffs.get(k,0) + v
        return Vector(coeffs)

    def __sub__(self, other):
        coeffs = dict(self.coeffs)
        for (k,v) in other.coeffs.items():
            coeffs[k] = coeffs.get(k,0) - v
        return Vector(coeffs)

    def __neg__(self):
        return (-1)*self

    def __rmul__(self, r):
        r = int(r)
        if r==0:
            return Vector({})
        return Vector({k:r*v for (k,v) in self.coeffs.items()})

    def __matmul__(self, other):
        if type(other) is int:
            return other*self
        assert isinstance(other, Vector)
        coeffs = {}
        for u in self:
          for v in other:
            uv = u+v # tuple concat
            coeffs[uv] = coeffs.get(uv, 0) + self[u]*other[v]
        return Vector(coeffs)

    def __rmatmul__(self, other):
        assert isinstance(other, int)
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

    def __call__(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        assert 0, "abstract base class"

    def call(self, v): # or use __mul__ ???
        m, n = self.shape
        u = Vector()
        for k in v:
            assert len(k) == n
            u = u + v[k]*self(*[Vector({ki:1}) for ki in k])
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
            u = self(*arg)
            for i in range(n):
                barg = list(arg)
                r = randint(-5,5)
                barg[i] = r*barg[i]
                assert self(*barg) == r*u
        # XX test additivity.. v+w .. 


class ILin(Lin):
    """ _identity Lin
    """
    def __init__(self, m=1):
        Lin.__init__(self, (m,m))

    def __call__(self, *vecs):
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

    def __call__(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        #print("__call__", list(vecs))
        perm = self.perm
        vecs = [vecs[i] for i in perm]
        v = reduce(matmul, vecs)
        assert isinstance(v, Vector)
        return v


class OpLin(Lin):
    def __init__(self, op, shape):
        Lin.__init__(self, shape)
        self.op = op

    def __call__(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        return self.op(*vecs)


class TensorLin(Lin):
    def __init__(self, aop, bop):
        m = aop.shape[0]+bop.shape[0]
        n = aop.shape[1]+bop.shape[1]
        Lin.__init__(self, (m,n))
        self.ops = (aop, bop)

    def __call__(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        aop, bop = self.ops
        a = vecs[:aop.shape[1]]
        b = vecs[aop.shape[1]:]
        a, b = aop(*a) , bop(*b)
        assert isinstance(a, (int, Vector)), (aop, a)
        assert isinstance(b, (int, Vector)), (bop, b)
        #print("TensorLin.__call__", a, b)
        return a@b


class ComposeLin(Lin):
    def __init__(self, aop, bop):
        assert aop.shape[1] == bop.shape[0]
        Lin.__init__(self, (aop.shape[0], bop.shape[1]))
        self.ops = (aop, bop)

    def __call__(self, *vecs):
        m, n = self.shape
        assert len(vecs) == n
        aop, bop = self.ops
        v = bop(*vecs)
        if type(v) is int:
            v = v*aop()
        else:
            v = aop.call(v)
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
        ident = unit()
        for u in vs:
            assert mul(ident, u) == u
            assert mul(u, ident) == u
        # assoc
        for u in vs:
          for v in vs:
           for w in vs:
            lhs = mul(mul(u,v), w)
            rhs = mul(u, mul(v, w))
            assert lhs == rhs
        if comm:
            for u in vs:
             for v in vs:
                assert mul(u,v) == mul(v,u)


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
            u = (I@counit).call(m)
            assert u==v
            u = (counit@I).call(m)
            assert u==v
            # coassoc
            mmi = (comul@I).call(m)
            imm = (I@comul).call(m)
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
            assert counit( mul(u,v) ) == counit(u)*counit(v)
        assert counit(unit()) == 1
        assert comul(unit()) == unit()@unit()

        I = ILin()
        S = PermLin([0,2,1,3])
        lin = (mul@mul)*S*(comul@comul)
        for u in vs:
          for v in vs:
            lhs = comul(mul(u,v))
            rhs = lin(u,v)
            assert lhs == rhs

        # hopf algebra
        antipode.test(vs)

        lhs = mul * (antipode @ I) * comul
        mid = unit * counit
        rhs = mul * (I @ antipode) * comul
        for u in vs:
            print()
            print("u =", u)
            print("comul(u) =", comul(u))
            m = mid.call(u)
            l = lhs.call(u)
            r = rhs.call(u)
            print("mid:", m)
            print("lhs:", l)
            print("rhs:", r)
            assert l==m
            assert r==m


        print("Hopf.test: OK")


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


def test():
    # ----------------------------------------------

    w = Vector({'abc':1})
    v = Vector({'ac':1})
    #assert str(w) == "Vector({'abc': 1})"
    assert str(w) == "Vector({('abc',): 1})"
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

    # ----------------------------------------------
    # the word Algebra
    unit = OpLin(lambda : Vector({'':1}), (1,0))
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

    assert mul(ab+a, c-bc) == -abbc + ac
    assert mul(unit(), ab+a) == ab+a
                
    algebra = Algebra(unit, mul)

    for u in [c,ac,abc]:
      for w in [ab,a,c]:
        vs.append(u+w)
    algebra.test(vs)

    for u in vs:
      for v in vs:
        lhs, rhs = mul(u,v) , mul.call(u@v)
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

    coalebra = Coalgebra(counit, copy)
    coalebra.test(vs)

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
    assert mul(ab, cd) == Vector(
        {'abcd': 1, 'acbd': 1, 'acdb': 1, 
        'cabd': 1, 'cadb': 1, 'cdab': 1})

    aaa = vector('aaa')
    assert mul(a,aaa) == 4*vector('aaaa')

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
    assert counit(2*unit()+a) == 2

    assert comul(abc) == Vector(
        {('', 'abc'): 1, 
        ('a', 'bc'): 1, 
        ('ab', 'c'): 1, 
        ('abc', ''): 1})

    coalebra = Coalgebra(counit, comul)
    coalebra.test(vs)
    

    def antipode(v):
        coeffs = {}
        for k in v:
            k = tuple(reversed(k))
            r = (-1)**len(k)
            coeffs[k] = v[k]*r
        return Vector(coeffs)
    antipode = OpLin(antipode, (1,1))

    hopf = Hopf(
        unit, mul,
        counit, comul,
        antipode)
    hopf.test(vs)
    



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



