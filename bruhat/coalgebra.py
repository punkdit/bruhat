#!/usr/bin/env python

"""
reference:
[HaGuKi]
_Algebras, Rings and _Modules 
Lie _Algebras and Hopf _Algebras (2009)
Michiel Hazewinkel
Nadiya Gubareni
V. V. Kirichenko

"""

import numpy


from bruhat.argv import argv


class Vector:
    "A dict with int values"
    def __init__(self, coeffs={}):
        coeffs = {k:int(v) for (k,v) in coeffs.items() if int(v)}
        keys = list(coeffs.keys())
        keys.sort()
        self.keys = keys
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



class Algebra:
    def __init__(self, unit, mul):
        self.unit = unit
        self.mul = mul

    def test(self, vs, commutative=False):
        unit, mul = self.unit, self.mul
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
        if commutative:
            for u in vs:
             for v in vs:
                assert mul(u,v) == mul(v,u)


class Coalgebra:
    def __init__(self, counit, comul):
        self.counit = counit
        self.comul = comul

    def test(self, vs):
        counit, comul = self.counit, self.comul


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

    w = Vector({'abc':1})
    v = Vector({'ac':1})
    assert str(w) == "Vector({'abc': 1})"
    assert (w+w) == 2*w
    assert 0*w == 0
    assert w-v != 0
    assert w-w == 0

    assert w+v == v+w
    assert -1*(w+v) == -w-v

    # the word Algebra
    unit = lambda : Vector({'':1})
    def mul(u, v):
        coeffs = {}
        for k in u:
            for j in v:
                kj = k+j
                coeffs[kj] = u[k] * v[j] + coeffs.get(kj, 0)
        return Vector(coeffs)

    vector = lambda s : Vector({s:1})
    vs = [vector(s) for s in 'ab a c ac bc abc abbc'.split()]
    ab,a,c,ac,bc,abc,abbc = vs

    assert mul(ab+a, c-bc) == -abbc + ac
    assert mul(unit(), ab+a) == ab+a
                
    algebra = Algebra(unit, mul)

    for u in [c,ac,abc]:
      for w in [ab,a,c]:
        vs.append(u+w)
    algebra.test(vs)


    # the shuffle Hopf algebra: [HaGuKi] Ex. 3.4.6

    unit = lambda : Vector({'':1})

    assert list(str_shuffle("abc", "ef")) == [
        'abcef', 'abecf', 'abefc', 'aebcf', 'aebfc', 'aefbc',
        'eabcf', 'eabfc', 'eafbc', 'efabc']

    def mul(u, v):
        coeffs = {}
        for k in u:
          for j in v:
            for kj in str_shuffle(k,j):
              coeffs[kj] = coeffs.get(kj,0) + u[k]*v[j]
        return Vector(coeffs)

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
        for k in u:
            for i in range(len(k)+1):
                key = (k[:i],k[i:])
                coeffs[key] = coeffs.get(key,0) + u[k]
        return Vector(coeffs)

    print(comul(abc))
    
    



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



