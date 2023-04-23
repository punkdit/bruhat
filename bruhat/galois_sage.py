#!/usr/bin/env python

"""
abelian galois extensions / cyclotomic number fields 
"""

from time import time
start_time = time()

from bruhat.argv import argv
from bruhat.smap import SMap

from sage.all_cmdline import *
from sage import all_cmdline 


def main():
    #n = argv.get("n", 63)
    #rows = argv.get("rows", 6)
    #cols = argv.get("cols", 6)
    #build(n, rows, cols)

    #build(17, 1, 16)
    #build(20, 2, 4)
    #build(44, 2, 10)
    build(63, 6, 6)


def build(n, rows, cols):
    K = CyclotomicField(n)
    G = K.galois_group()
    Hs = G.subgroups()
    for H in Hs:
        if len(H) != 18:
            continue
        K0 = H.fixed_field()[0]
        #print(K0)
        #print(K0.defining_polynomial())

    items = []
    for H in Hs:
        gens = H.gens()
        s = SMap()
        for row in range(rows):
             for col in range(cols):
                s[row, col] = '.'
        for g in H:
          cs = list(g.cycles()) or [(1,)]
          c = cs[0]
          c = eval(str(c)) # wtf
          for ii, i in enumerate(c):
              i -= 1 # one-based index
              col = i%cols
              row = i//cols
              if ii in [0,1]:
                s[row, col] = "X"
        print(s)
        K0 = H.fixed_field()[0]
        print(K0.defining_polynomial())
        print(H.fixed_field())
        print()
        item = (len(H), s, K0.defining_polynomial())
        items.append(item)

    N = 1
    col = 0
    smap = SMap()
    for item in items:
        if item[0] == N:
            smap.paste(item[1], 0, col)
        else:
            print(smap, rows*cols//N, '\n')
            N = item[0]
            smap = SMap()
            col = 0
            smap.paste(item[1], 0, col)
        col += cols+1

    print(smap, rows*cols//N, '\n')
    





if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


