#!/usr/bin/env python3

from operator import mul
from functools import reduce

import numpy


from bruhat.poly import Poly, Q
from bruhat.util import cross, all_subsets, factorial
from bruhat.argv import argv
from bruhat.elim import solve, shortstr


ring = Q
a = Poly("a", ring)


def interp(vals):
    n = len(vals)

    p = Poly({}, ring)
    one = Poly({():ring.one}, ring)
    
    for i in range(n):
        #print(i, vals[i])
        y0 = p(a=i)
        q = one
        for j in range(i):
            q = q*(a-j)
        #print("q:", q)
        for j in range(i):
            assert q(a=j) == 0
        #print("q(%d)=%s"%(i, q(a=i)))
        r = ring.one / factorial(i)
        #print("r =", r)
        if y0 != vals[i]:
            p = p + (vals[i] - y0)*r*q
        #print("p(%d)=%s"%(i, p(a=i)))
        #print()
    return p


def multi_interp(target):
    shape = target.shape
    n = len(shape)
    #print("multi_interp", shape)
    vs = 'abcde'[:n]
    ms = [Poly(v, ring) for v in vs]
    #print(ms)
    itemss = [list(range(i)) for i in shape]
    coords = []
    for idxs in cross(itemss):
        namespace = dict((vs[i], idxs[i]) for i in range(n))
        #print(namespace)
        coords.append(namespace)
    A = []
    polys = []
    for idxs in cross(itemss):
        p = reduce(mul, [p**i for (p,i) in zip(ms, idxs)])
        polys.append(p)
        #print(idxs, p)
        row = []
        for coord in coords:
            v = p.substitute(coord)
            row.append(ring.promote(v))
        A.append(row)
    A = numpy.array(A, dtype=object)
    A = A.transpose()
    #print(A.shape)
    #print(shortstr(A))
    rhs = target.view()
    rhs.shape = (len(A),1)
    #print(rhs)

    print("solve...")
    v = solve(ring, A, rhs)
    assert v is not None
    print(shortstr(v))
    q = ring.zero
    for i, p in enumerate(polys):
        q = q + v[i, 0]*p
    return q


def multi_factorize(p, N=10, denom=2):
    vs = p.get_vars()
    #print("multi_factorize", vs)
    ring = p.ring
    d = p.degree
    factors = []

    idxss = list(all_subsets(len(vs)))
    idxss.sort(key = len)
    assert idxss[0] == []
    idxss.pop(0)
    for idxs in idxss:
        subvs = [vs[idx] for idx in idxs]
        print("subvs:", subvs)
        coords = [[-ring.promote(x)/denom for x in range(N)] for v in subvs]
        for ii in cross(coords):
            kw = dict((subvs[i], ii[i]) for i in range(len(subvs)))
            y = p(**kw)
            if y!=0:
                continue
            q = ring.zero
            for k,v in kw.items():
                #print("\t", k, v)
                q += Poly(k, ring) - v
            while 1:
                print("factor:", q)
                div, rem = q.reduce(p)
                if rem != 0:
                    break
                factors.append(q)
                p = div
                print("\t", p)
            if p.degree == 1:
                break

    if p != 1:
        factors.append(p)
    return factors


def factorize(p):
    ring = p.ring
    d = p.degree
    factors = []
    for i in range(6*20):
        i = ring.one*i/6
        y = p(a=-i)
        if y!=0:
            continue
        while 1:
            f = (a+i)
            div, rem = f.reduce(p)
            if rem != 0:
                break
            factors.append(a+i)
            p = div
    if p != 1:
        factors.append(p)
    return factors


if argv.vals is not None:
    vals = argv.get("vals", [1, 4, 10, 20, 35, 56])
    p = interp(vals)
    
    print("p =", p)
    print("degree =", p.degree)
    print("factors:", factorize(p))
    
    #print(["%s"%p(a=i) for i in range(n)])
    #print([(a-i).reduce(p)[1] for i in range(n)])

elif 1:

    if 0:
        # B2
        vals = numpy.array(
            [[1, 10, 35, 84, 165], 
            [5, 35, 105, 231, 429], 
            [14, 81, 220, 455, 810], 
            [30, 154, 390, 770, 1326], 
            [55, 260, 625, 1190, 1995]])

    name = argv.next() or "A2"
    N = int(argv.next() or 4)
    import os
    data = os.popen("./sl.sage %s %s"%(name, N)).read()
    vals = eval(data)

    vals = numpy.array(vals)
    vals = vals[:,0,:,0]
    vals = vals.copy()
    print(vals)

    p = multi_interp(vals)
    print("degree:", p.degree)
    print(p)

    #factors = multi_factorize(p)
    #print(factors)
    

elif 0:
    f = lambda a, b, c : (a+1)*(b+1)*(c+1)*(a+b+2)*(b+c+2)*(a+b+c+3)//12
    
    N = 5
    for c in range(3):
        for b in range(N):
            for a in range(N):
                print("%6s"%f(a, b, c), end=" ")
            print()
        print()
    
    
elif 0:
    f = lambda a, b, c, d : (
        (a+1)*(b+1)*(c+1)*(d+1)*
        (a+b+2)*(b+c+2)*(c+d+2)*
        (a+b+c+3)*(b+c+d+3)*
        (a+b+c+d+4)
        //288)
    
    
    N = 5
    for d in range(3):
     for c in range(3):
      for b in range(N):
       for a in range(N):
        print("%6s"%f(a, b, c, d), end=" ")
       print()
      print()
     print()

