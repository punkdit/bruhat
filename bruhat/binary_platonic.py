#!/usr/bin/env python

"""
see:

An Introduction to the McKay Correspondence
Master Thesis in Physics
Author : Max Lindh
"""

from bruhat.word import Gen, build
from bruhat.repr_sage import dixon_irr
from bruhat.argv import argv


def get_BT():
    a, b, c = [Gen(c) for c in "abc"]
    abc = a*b*c
    graph = build([a, b, c], [a**2==abc, b**3==abc, c**3==abc])
    G = graph.get_gset()
    assert len(G) == 24
    return G
    
def get_BO():
    a, b, c = [Gen(c) for c in "abc"]
    abc = a*b*c
    graph = build([a, b, c], [a**4==abc, b**3==abc, c**2==abc])
    G = graph.get_gset()
    assert len(G) == 48
    return G

def get_BI():
    a, b, c = [Gen(c) for c in "abc"]
    abc = a*b*c
    graph = build([a, b, c], [a**2==abc, b**3==abc, c**5==abc])
    G = graph.get_gset()
    assert len(G) == 120
    return G
    



def test():

    G = get_BT()
    table = dixon_irr(G)
    chi = table[4]

#    G = get_BO()
#    table = dixon_irr(G)
#    chi = table[3]
#
#    G = get_BI()
#    table = dixon_irr(G)
#    chi = table[1]

    print(G)
    print(table)

    print(chi)
    c = chi
    for i in range(1,8):
        print(i, [c.dot(a) for a in table])
        c = chi*c



if __name__ == "__main__":

    from time import time
    start_time = time()

    profile = argv.profile
    name = argv.next() or "test"
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        print("%s()"%name)
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        print("%s()"%name)
        fn = eval(name)
        fn()

    else:
        test()


    t = time() - start_time
    print("OK! finished in %.3f seconds\n"%t)



