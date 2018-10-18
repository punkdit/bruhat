#!/usr/bin/env python3

"""
Build Reed-Muller codes.
"""


import numpy

from solve import array2, zeros2, dot2, shortstr, rank, find_kernel, span
from argv import argv
from util import choose


class Code(object):
    def __init__(self, G, H=None, d=None, check=True):
        self.G = G.copy()
        self.k, self.n = G.shape
        self.d = d

        if H is None:
            H = list(find_kernel(G))
            H = array2(H)
        if H.shape == (0,):
            H.shape = (0, self.n)
        self.H = H.copy()

        if check:
            self.check()

    def check(self):
        G, H = self.G, self.H
        assert rank(G) == len(G)
        A = dot2(H, G.transpose())
        assert A.sum()==0

    def __str__(self):
        return "Code([[%s, %s, %s]])" % (self.n, self.k, self.d)

    def dump(self):
        G, H = self.G, self.H
        print("G =")
        print(shortstr(G))
    
        print("H =")
        print(shortstr(H))

    def is_selfdual(self):
        G = self.G
        x = dot2(G, G.transpose())
        return x==0

    def get_distance(self):
        G = self.G
        d = None
        for v in span(G):
            w = v.sum()
            if w==0:
                continue
            if d is None or w<d:
                d = w
        return d

    def puncture(self, i):
        assert 0<=i<self.n
        G = self.G
        A = G[:i]
        B = G[i+1:]
        print(A.shape)
        print(B.shape)
        G = numpy.concatenate((A, B), axis=0)
        return Code(G)
    

def build_rm(r, m):

    assert 0<=r<=m

    n = 2**m # length

    one = array2([1]*n)
    basis = [one]

    vs = [[] for i in range(m)]
    for i in range(2**m):
        for j in range(m):
            vs[j].append(i%2)
            i >>= 1
        assert i==0

    vs = [array2(v) for v in vs]

    for k in range(r):
        for items in choose(vs, k+1):
            v = one
            #print(items)
            for u in items:
                v = v*u
            basis.append(v)
        
    G = numpy.array(basis)
    code = Code(G)
    return code


def test():

    for m in range(2, 6):
      for r in range(0, m+1):
        code = build_rm(r, m)
        assert code.n == 2**m
        k = 1
        for i in range(1, r+1):
            k += len(list(choose(list(range(m)), i)))
        assert code.k == k
        if code.k < 12:
            assert code.get_distance() == 2**(m-r)

    print("OK")


def main():
    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)

    code = build_rm(r, m)

    print(code)
    print("d =", code.get_distance())

    code = code.puncture(3)

    print(code)
    print("d =", code.get_distance())




if __name__ == "__main__":

    if argv.main:
        main()
    else:
        test()



