#!/usr/bin/env python

import sys
from random import *

from argv import argv



def count_orbits():
    from action import Group, Perm
    from bruhat import incidence as geometry

    X = geometry.projective(3)
    graph = X.get_bag()

    for p in graph:
        print p
    print

    points = [p for p in graph if p.desc=='0']
    lines = [p for p in graph if p.desc=='1']
    assert len(points)==7
    assert len(lines)==7

    p0 = points[0]
    p1 = points[1]
    for line in lines:
        if p0 in line.nbd and p1 in line.nbd:
            break
    else:
        assert 0
    for i in range(2, len(points)):
        p2 = points[i]
        if p2 not in line.nbd:
            break
    else:
        assert 0
    #print p0, p1, p2

    triangle = [p0.idx, p1.idx, p2.idx]

    #for f in X.get_symmetry():
    #    print f
    items = range(14)
    perms = [f for f in X.get_symmetry()]
    perms = [Perm(f, items) for f in perms]
    G = Group(perms, items)

    print "G:", len(G)

    H = G.invariant(*triangle)
    print "H:", len(H)

    cosets = G.left_cosets(H)
    A = G.left_action(cosets)

    B = A.pushout(A)

    print "A*A =", len(B.components())

    print "A*A*A ="
    for B1 in B.components():
        C = A.pushout(B1)
        print '\t',len(C.components())





def main():


    # dimension
    n = argv.get("n", 3)

    def bstr(p):
        s = []
        for i in range(n):
            s.append('1' if p&1 else '0')
            p >>= 1
        s = ''.join(reversed(s))
        return s


    points = []
    for i in range(1, 2**n):
        points.append(i)
    #for p in points:
    #    print bstr(p)

    lines = set()
    for p in points:
      for q in points:
        if p==q:
            continue
        r = p^q
        assert r!=0
        line = [p, q, r]
        line.sort()
        lines.add(tuple(line))
    assert len(lines)==7
    print "lines:", lines

    N = len(points)
    bases = []
    for p in points:
      for q in points:
        if p==q:
            continue
        for r in points:
            if p==r or q==r:
                continue
            base = [p, q, r]
            figure = list(base)
            figure.sort()
            if tuple(figure) in lines:
                continue
            bases.append(tuple(base))
    print "bases:", len(bases)

#    configs = set()
#    for b0 in bases:
#     for b1 in bases:
#      for b2 in bases:
#        configs.add((b0, b1, b2))
#    print "configs:", len(configs)

    b0 = bases[0] # pick one...

#    b0 = choice(bases)
#    b1 = choice(bases)
#    b2 = choice(bases)
#
#    print b0
#    print b1
#    print b2
#    print

    stats = {}
    for b1 in bases:
     for b2 in bases:
      i0 = 0
      j0 = 1
      k0 = 2
      solutions = 0
      for i1 in range(3):
        if b1[i1] == b0[i0]:
            continue
        for i2 in range(3):
            triple = [b0[i0], b1[i1], b2[i2]]
            if len(set(triple))!=3:
                continue
            #print "-----", triple
            figure = list(triple)
            figure.sort()
            if tuple(figure) in lines:
                continue
#            print triple

            for j1 in range(3):
              if j1==i1:
               continue
              for j2 in range(3):
                if j2==i2:
                    continue
                
                triple = [b0[j0], b1[j1], b2[j2]]
                if len(set(triple))!=3:
                    continue
                #print "-----", triple
                figure = list(triple)
                figure.sort()
                if tuple(figure) in lines:
                    continue
#                print "\t", triple

                for k1 in range(3):
                  if k1==i1 or k1==j1:
                   continue
                  for k2 in range(3):
                    if k2==i2 or k2==j2:
                        continue
                    
                    triple = [b0[k0], b1[k1], b2[k2]]
                    if len(set(triple))!=3:
                        continue
                    #print "-----", triple
                    figure = list(triple)
                    figure.sort()
                    if tuple(figure) in lines:
                        continue
#                    print "\t\t", triple
                    solutions += 1

      #print "solutions:", solutions
      stats[solutions] = stats.get(solutions, 0) + 1

    print stats
#    for value in stats.values():
#        print value, 1.*value / 21


#        print b0, b1, b2
#        count = 0
#        for i in range(3):



if __name__ == "__main__":

    _seed = argv.seed
    if _seed is not None:
        seed(_seed)

    fn = argv.next()

    if fn is not None:
        fn = eval(fn)
        fn()

    else:
        main()



