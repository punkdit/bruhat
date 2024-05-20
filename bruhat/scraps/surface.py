#!/usr/bin/env python

"""
build some cellulated surfaces... fail
"""

from time import time
from random import choice, randint

from bruhat.todd_coxeter import Schreier

from bruhat.argv import argv


def main():
    ngens = 3
    r, g, b = (0, 1, 2)
    v, e, f = (0, 1, 2)
    rels = [(v, v), (e, e), (f, f)] # self-inverse
    rels.append((f,v,f,v)) # edges

    def make(hgens, maxsize=10000, rev=False):
        hgens = [tuple('vef'.index(item) for item in gen) for gen in hgens]
        if rev:
            hgens = [tuple(reversed(gen)) for gen in hgens]
        graph = Schreier(ngens, rels)
        graph.build(hgens, maxsize)
        return graph

#    # The dessin for the Gaussian elliptic curve
#    hgens = ["rgrb", "rb"*4, "rbgbrb", "grbr", "bgbg"]
#    graph = make(hgens, rev=True)
#    assert len(graph) == 8
#
#    rev = lambda gen : tuple(reversed(gen))
#
#    # this is dessin/code #1 above:
#    hgens = [
#        "bg"*8,  # around red
#        "rb"*4, # around green
#        "grgr", # around blue
#        "grbrbg", # green 4
#        "bgrgrgrgrb", # blue 5
#        rev("bgbrgb"), # green 13
#        rev("rbrgbrgrbr"), # green 5
#        "gbgbrb",  # homology cycle
#        "rbrgbgbgbg", # homology cycle
#        #"gbgrgbgbgbgb", # another homology cycle.. not needed
#    ]
#    graph = make(hgens, rev=True)
#    assert len(graph) == 16

    # ['vefvvev', 'vvvffv', 'evefefvfv', 'vffevefvee', 'efefefvv']
    # ['vevfee', 'efevevevvv', 'efve', 'fefffvvvefef']

    #while 1:
    for _ in range(10000):

        hgens = []
        for i in range(10):
            rel = ''.join(choice('vef') for j in range(randint(4,12)))
            hgens.append(rel)
            graph = make(list(hgens))
            n = len(graph)
            if n < 1000:
                break
        else:
            continue
        if n > 10:
            print()
            print(hgens)
            print(n)
        else:
            print(".", end='', flush=True)

    print()



if __name__ == "__main__":

    start_time = time()


    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        main()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")

