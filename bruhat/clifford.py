#!/usr/bin/env python

from bruhat.action import mulclose
from bruhat.argv import argv


class Atom(object):
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
        return Compound([self, other])


class Compound(Atom):
    def __init__(self, atoms):
        _atoms = []
        for atom in atoms:
            if isinstance(atom, Compound):
                _atoms += atom.atoms
            else:
                _atoms.append(atom)
        self.atoms = _atoms
        name = "*".join(str(atom) for atom in atoms)
        Atom.__init__(self, name)


class Cyclotomic(Atom):
    def __init__(self, i, p):
        self.i = i%p
        self.p = p
        #name = "zeta%d^%d"%(p, i)
        name = "w(%d,%d)"%(i,p)
        Atom.__init__(self, name)

    def __mul__(self, other):
        if isinstance(other, Cyclotomic):
            assert self.p == other.p
            atom = Cyclotomic(self.i+other.i, self.p) 
        else:
            atom = Atom.__mul__(self, other)
        return atom


#class Pauli(Atom):
#    def __init__(self, zx):



def test():

    U = Atom("U")
    V = Atom("V")

    assert U==U
    assert U*V != V*U
    assert str(U*V) == "U*V"

    p = 8
    I = Cyclotomic(0, p)
    w = Cyclotomic(1, p)
    assert I*I == I
    assert I*w == w
    assert len(mulclose([w])) == p
    




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





