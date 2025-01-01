#!/usr/bin/env python

from sage.all_cmdline import GL, matrix, Matrix, GF

from bruhat.argv import argv
from bruhat.equ import Equ
from bruhat.algebraic import Algebraic

def main():
    p = argv.get("p", 2)
    n = argv.get("n", 3)

    G = GL(n,p)
    F = GF(2)

    print("|G| =", len(G))

    if 0:
        found = set()
        for g in G:
            A = Matrix(g)
            char = A.characteristic_polynomial()
            s = str(char)
            found.add(char)
    
        print(len(found))
        for char in found:
            print(char, ":", char.factor())

    G = Algebraic.GL(n, p)
    inv = {g:g.inverse() for g in G}

    remain = set(G)
    itemss = []
    while remain:
        g = iter(remain).__next__()
        remain.remove(g)
        items = [g]
        for h in G:
            k = inv[h]*g*h
            if k in remain:
                remain.remove(k)
                items.append(k)
        itemss.append(items)

    print("char classes:", len(itemss))

    itemss.sort(key = len)
    print([len(items) for items in itemss], len(itemss))

    for items in itemss:
        g = items[0]
        A = matrix(F, n, n, g.A)
        char = A.characteristic_polynomial()
        print(len(items), "-->", char, "\t", char.factor())
        for h in items[1:]:
            B = matrix(F, n, n, h.A)
            assert char == B.characteristic_polynomial()






if __name__ == "__main__":
    from time import time
    start_time = time()
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))



