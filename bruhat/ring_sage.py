#!/usr/bin/env python

"""

"""

from time import time
start_time = time()
import string

from sage.all_cmdline import *
from sage import all_cmdline 

from bruhat.argv import argv
from bruhat.smap import SMap



def main():
    n = argv.get("n")

    K0 = NumberField(x**2 + 1, "x")
    print(K0)

    K1 = K0.extension(x**2+x+1, "w")
    print(K1)

    print(K1.absolute_field("a").galois_group())


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


