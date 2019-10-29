#!/usr/bin/env python2

from __future__ import print_function

from bruhat.dev.sage_env import *

#from bruhat.argv import argv
#d = argv.get("d", 2) # qudit dimension

from sage.groups.matrix_gps.orthogonal  import SO

def pystr(M):
    s = str(M)
    s = s.replace("\n", ",")
    s = s.replace(" ", ",")
    return "["+s+"]"

SPACE = " "*4
for e in [1, -1]:
 if e==-1:
    continue # skip this 
 for p in [2, 3, 5]:
  for n in [2, 4, 6, 8, 10]:
    G = SO(n, p, e)
    
    ident = SPACE
    print(ident+"@classmethod")
    if e==1:
        print(ident+"def SO_%d_%d(cls, **kw):" % (n, p))
    else:
        print(ident+"def SO_%d_%d_1(cls, **kw):" % (n, p))
    ident = 2*SPACE
    print(ident+"return cls.make(")
    ident = 3*SPACE
    gens = "[%s]"%((",\n"+ident+" ").join(pystr(g) for g in G.gens()))
    print(ident + gens + ",")
    print(ident+pystr(G.invariant_bilinear_form())+",")
    print(ident+pystr(G.invariant_quadratic_form())+",")
    try:
        order = len(G)
    except OverflowError:
        order = None
    print(ident+str(order), ",", p, ")")
    print()





