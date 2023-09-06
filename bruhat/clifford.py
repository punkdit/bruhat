#!/usr/bin/env python


"""
Implement Bruhat decomposition of the Clifford group
as described here:

Hadamard-free circuits expose the structure of the Clifford group
Sergey Bravyi, Dmitri Maslov
https://arxiv.org/abs/2003.09412

maybe ?

"""



from bruhat.action import mulclose
from bruhat.clifford_sage import Clifford
from bruhat.argv import argv


class Gen(object):
    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __hash__(self):
        return hash(self.name)

    def __mul__(self, other):
        return Word([self, other])


class Word(Gen):
    def __init__(self, gens):
        _gens = []
        for gen in gens:
            if isinstance(gen, Word):
                _gens += gen.gens
            else:
                _gens.append(gen)
        self.gens = _gens
        name = "*".join(str(gen) for gen in gens)
        Gen.__init__(self, name)



def test():

    U = Gen("U")
    V = Gen("V")

    assert U==U
    assert U*V != V*U
    assert str(U*V) == "U*V"



if __name__ == "__main__":

    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "test"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()





