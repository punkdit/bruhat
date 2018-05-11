#!/usr/bin/env python3

import sys, os

from heapq import *

from smap import SMap
from argv import argv


class Point(object):
    def __init__(self, *coords):
        self.coords = tuple(coords)

    def __str__(self):
        return "Point%s"%(self.coords,)
    __repr__ = __str__

    def __hash__(self):
        return hash(self.coords)

    def __add__(self, other):
        cs = tuple(aa+bb for (aa, bb) in zip(self.coords, other.coords))
        return Point(*cs)

    def __sub__(self, other):
        cs = tuple(aa-bb for (aa, bb) in zip(self.coords, other.coords))
        return Point(*cs)

    def __neg__(self):
        cs = tuple(-a for a in self.coords)
        return Point(*cs)

    def __rmul__(self, r):
        assert int(r) == r
        cs = tuple(r*a for a in self.coords)
        return Point(*cs)

    def __eq__(self, other):
        return self.coords == other.coords

    def __ne__(self, other):
        return self.coords != other.coords

    def __le__(self, other):
        return self.norm() <= other.norm()

    def __lt__(self, other):
        return self.norm() < other.norm()

    def __getitem__(self, idx):
        return self.coords[idx]

    def norm(self):
        "norm squared"
        return sum(a**2 for a in self.coords)

    def nbd(self):
        coords = self.coords
        n = len(coords)
        items = []
        for i in range(n):
            cs = list(coords)
            cs[i] = cs[i]-1
            items.append(Point(*cs))
            cs[i] = cs[i]+2
            items.append(Point(*cs))
        return items


def test():

    zero = Point(0, 0, 0)
    a = Point(1, 0, 0)
    b = Point(0, 1, 0)
    c = Point(0, 0, 1)

    assert zero+b==b
    assert zero-b!=b
    assert a-a==zero

    assert zero <= a
    assert a <= b


    nbd = set(c.nbd())
    assert len(nbd)==6
    for p in nbd:
        assert (p-c).norm() == 1


def nsquares(n, ncoeffs):
    zero = (0,)*n
    zero = Point(*zero)

    items = [zero]
    found = set(items)
    counts = {}
    n = 0
    counts[n] = 0
    #for i in range(100):
    while len(counts) <= ncoeffs:
        a = heappop(items)
        assert a in found
        n = a.norm()
        counts[n] = counts.get(n, 0) + 1
        for b in a.nbd():
            if b in found:
                continue
            found.add(b)
            heappush(items, b)

    del counts[n]
    #for p in items:
    #    print(p, p.norm())
    #print

    keys = list(counts.keys())
    keys.sort()
#    items = []
#    for key in keys:
#        items.append((key, counts[key]))

    items = [0]*(keys[-1]+1)
    for key in keys:
        items[key] = counts[key]
    return items



def all_primes(n, ps=None):
    "list of primes < n"

    items = [0]*n
    p = 2 

    while p**2 < n:
        i = 2 
        while p*i < n:
            items[p*i] = 1 
            i += 1

        p += 1
        while p < n and items[p]:
            p += 1

    ps = [i for i in range(2, n) if items[i]==0]
    return ps


def divisors(n):
    divs = [1]
    for i in range(2, n):
        if n%i == 0:
            divs.append(i)
    if n>1:
        divs.append(n)
    return divs


def sigma(k, n):
    "sum of k-th powers of divisors of n"
    i = 0
    for j in divisors(n):
        i += j**k
    return i


def sigma_quo(k, n, quo):
    "sum of k-th powers of divisors of n, whose quotient is divisable only by quo"
#    print("sigma_quo(%s, %s, %s)"%(k, n, quo))
    i = 0
    for j in divisors(n):
        a = n//j
#        print("\tj = %d, n//j = %d"%(j, a), end=' ')
        allow = True
        # a must be divisible by each element of quo
        for b in quo:
            if a%b:
                allow = False
                #print("A", end="")
        # and no other primes
        for b in all_primes(a):
            if a%b==0 and b not in quo:
                allow = False
#                print("(b=%d not in quo)"%b, end="")
#        print("OK" if allow else "")
        if allow:
            i += j**k
#            print("\tj**k", j**k)
    return i


def eisenstein(k, n):
    "E_k eisenstein series: sum of (k-1)-th powers of n"
    assert k>=1
    items = [sigma(k-1, i) for i in range(1, n)]
    return items


def index2(a, b):

    N = 10
    count = 0

    for i in range(-2*N,2*N):
      for j in range(-2*N,2*N):
        p = i*a + j*b
        if 0<=p[0]<N and 0<=p[1]<N:
            count += 1

    d = int(round(N*N / count))
    return d


def lattice_shape(a, b, N):
    items = []
    for i in range(-4*N,4*N):
      for j in range(-4*N,4*N):
        p = i*a + j*b
        if 0<=p[0]<N and 0<=p[1]<N:
            items.append(p)
    items = list(set(items))
    return items


def lattice_draw(points):
    a0 = min(p[0] for p in points)
    assert a0==0
    a1 = max(p[0] for p in points)
    b0 = min(p[1] for p in points)
    assert b0==0
    b1 = max(p[1] for p in points)

    smap = SMap()
    for i in range(a1+1):
      for j in range(b1+1):
        smap[i, j] = '.'
    for p in points:
        smap[p[0], p[1]] = "*"
    print(smap)
    #for p in points:
    #    print(p.coords, end=" ")
    print()


def main():

    # find sublattices of a given index
    n = argv.get("n", 4)

    N = 2*n-1
    uniq = set()
    for i0 in range(0, n+1):
      for i1 in range(0, n+1):
        for j0 in range(0, n+1):
          for j1 in range(0, n+1):
            #j1 = ( n + i1*j0 ) // i0
            area = i0*j1 - i1*j0
            if area != n:
                continue
            items = lattice_shape(Point(i0, i1), Point(j0, j1), N)
            items.sort(key = lambda p:p.coords)
            icon = str(items)
            if icon not in uniq:
                print("%s*A + %s*B, %s*A + %s*B" % (i0, i1, j0, j1))
                print(area)
                lattice_draw(items)
                uniq.add(icon)
    print("count =", len(uniq))


def main():

    assert sigma(0, 4) == 3
    assert sigma(1, 6) == 12

    print(eisenstein(1, 20))
    print(eisenstein(2, 40))
    print(eisenstein(4, 20))

    print()

    print([i//4 for i in nsquares(2, 20)][1:])
    #print([i//4 for i in nsquares(3, 20)][1:])
    print([i//8 for i in nsquares(4, 40)][1:])
    #print([i//16 for i in nsquares(8, 20)][1:])


def main():

    # another experiment...

    N = argv.get("N", 20)
    quo = argv.get("quo", [2])
    #print("quo = %s" % quo)

    E = eisenstein(4, N)
    print(E)

    a = 1
    for b in quo:
        a *= b

    items = []
    for i in range(1, N):
        j = sigma_quo(3, i*a, quo)
        #print("i=%d: %d" % (i*a, j))
        items.append(j)
    print(items)

    glitch = [E[i] - items[i] for i in range(N-1)]
    print(glitch)



if __name__ == "__main__":
    main()





