#!/usr/bin/env python

"""
here we used mobius inversion to find formulas for the trace in number fields
"""

import string

from sage.all_cmdline import *
from sage import all_cmdline

from bruhat.argv import argv
from bruhat.smap import SMap


def main():

    degree = 3
    if degree == 2:
        n = 3
        K = CyclotomicField(n)
        G = K.galois_group()
        print(K)

    elif degree == 3:
        x = polygen(QQ, 'x')
        K = NumberField(x**3 + 5, "x")
    else:
        assert 0

    x = K.gens()[0]
    one = K.one()
    G = K.galois_group()

    print(K)
    print(G)
    A = K.trace_pairing([1,x,x**2])
    #print((x**2-x+1).trace())
    #print(A)
    #print(A[1,2])

    def trace(a, b):
        return K.trace_pairing([a,b])[0,1]
    #t = lambda *args : reduce(trace, args)
    def t(*args):
        n = len(args)
        assert n
        if n==1:
            return (one*args[0]).trace()
        args = list(args)
        while len(args)>1:
            args = [trace(args[0], args[1])] + args[2:]
        return args[0]

    f = lambda a,b,c:(2*t(a*b*c) - t(a*b)*t(c) - t(a*c)*t(b) - t(a)*t(b*c) + t(a)*t(b)*t(c))
    if degree == 2:
        assert t(1) == 2
        assert t(3) == 6
        assert f(1, x, x+1) == 0

    f = lambda a,b,c,d: (
        +t(a)*t(b)*t(c)*t(d)
        -t(a*b)*t(c)*t(d) -t(a*c)*t(b)*t(d) -t(a*d)*t(b)*t(c) 
            -t(a)*t(b*c)*t(d) -t(a)*t(b*d)*t(c) -t(a)*t(b)*t(c*d)
        +t(a*b)*t(c*d) + t(a*c)*t(b*d) +t(a*d)*t(b*c)
        +2*t(a*b*c)*t(d) +2*t(a*b*d)*t(c) + 2*t(a*c*d)*t(b) + 2*t(a)*t(b*c*d)
        -6*t(a*b*c*d)
    )

    els = [2, x, x**2, x+x**2, 1+x+x**2]
    for a in els:
     for b in els:
      for c in els:
       for d in els:
        assert f(a,b,c,d) == 0


if __name__ == "__main__":

    from time import time
    start_time = time()

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



