#!/usr/bin/env python

"""
Commutative polynomials.

"""

import sys, os

#from bruhat.action import Group



def tadd(a, b):
    "tuple add"
    c = tuple(x+y for (x,y) in zip(a, b))
    return c


class Poly(object):
    """
        Commutative polynomial in <rank> many variables x_1,...,x_<rank> .
    """

    def __init__(self, cs, rank=None):
        coefs = {}
        for key, value in cs.items():
            assert rank is None or rank == len(key)
            rank = len(key)
            assert long(value) == value
            if value != 0:
                coefs[key] = value
        self.cs = coefs
        self.rank = rank

    @classmethod
    def identity(cls, rank):
        key = (0,)*rank
        return cls({key : 1})

    def __eq__(self, other):
        assert self.rank == other.rank
        return self.cs == other.cs

    def __ne__(self, other):
        assert self.rank == other.rank
        return self.cs != other.cs

    def __hash__(self):
        cs = list(self.cs.items())
        cs.sort(key = lambda (k,v):k)
        cs = tuple(cs)
        return hash(cs)

    def __getitem__(self, key):
        assert len(key) == self.rank
        assert isinstance(key, tuple)
        return self.cs.get(key, 0)

    def __add__(self, other):
        assert self.rank == other.rank
        cs = dict(self.cs)
        for key, value in other.cs.items():
            cs[key] = cs.get(key, 0) + value
        return Poly(cs)

    def __sub__(self, other):
        assert self.rank == other.rank
        cs = dict(self.cs)
        for key, value in other.cs.items():
            cs[key] = cs.get(key, 0) - value
        return Poly(cs)

    def __rmul__(self, r):
        cs = {}
        for key, value in self.cs.items():
            cs[key] = r * value
        return Poly(cs)

    def __mul__(self, other):
        assert self.rank == other.rank
        cs = {}
        for k1, v1 in self.cs.items():
          for k2, v2 in other.cs.items():
            k = tadd(k1, k2)
            cs[k] = cs.get(k, 0) + v1*v2
        return Poly(cs)

    def __pow__(self, n):
        if n==0:
            return Poly.identity(self.rank)
        assert n>0
        p = self
        for i in range(n-1):
            p = self*p
        return p

    def __str__(self):
        rank = self.rank
#        if rank==0:
#            return 
        if rank>1:
            names = ["x_%d"%(i+1) for i in range(rank)]
        else:
            names = ["x"]
        items = []
        cs = list(self.cs.items())
        cs.sort(key = lambda (k,v):list(reversed(k)))
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
                    item.append("%s^%d" % (names[i], exp))
            item = ' '.join(item)
            if not item:
                item = str(val)
            elif val==1:
                pass
            elif val==-1:
                item = "-"+item
            else:
                item = "%s*%s" % (val, item)
            items.append(item)
        s = ' + '.join(items)
        s = s.replace("+ -", "- ")
        return s

    def is_symmetric(self):
        n = self.rank
        ns = range(n)
        cs = self.cs

        #keys = list(cs.keys())
        for k0, v0 in cs.items():

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


def main():


    I = Poly({(0,) : 1})
    x = Poly({(1,) : 1})

    assert Poly.identity(1) == I

    assert (I + x*x) * (I+x*x) == (I+x**2)**2
    assert (I + x*x) * (I+x*x) != (I+x)**2
    assert (I+x*x) * (I-x*x) == I-x**4
    assert (I-x**4).is_symmetric()
    assert x*x == x**2

    assert (I+x)**4 == I + 4*x + 6*x**2 + 4*x**3 + x**4

    # ----------------

    I = Poly.identity(3)
    x1 = Poly({(1,0,0) : 1})
    x2 = Poly({(0,1,0) : 1})
    x3 = Poly({(0,0,1) : 1})

    #print (x1+x2+x3)**3
    #print x1*x2*x3

    assert ((x1+x2+x3)**2 + I).is_symmetric()
    assert (x1*x2*x3).is_symmetric()

    assert not ((x1+x2+x3)**3 + x1).is_symmetric()
    assert not (x1*x2*x3 + x1*x2).is_symmetric()



if __name__ == "__main__":

    main()


