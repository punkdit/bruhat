#!/usr/bin/env python

from bruhat.argv import argv
from bruhat.todd_coxeter import Schreier


class Word(object):
    def __init__(self, items):
        assert isinstance(items, tuple)
        self.items = items
    def __str__(self):
        return '*'.join(self.items) or "1"
    def __eq__(self, other): # tricky, take care!
        return Relation(self, other)
#    def __eq__(self, other):
#        return self.items == other.items
#    def __hash__(self):
#        return hash(self.items)
    #def __pow__(self, n):
    def __mul__(self, other):
        assert isinstance(other, Word)
        return Word(self.items + other.items)
    def __getitem__(self, idx):
        return self.items[idx]
    def __len__(self):
        return len(self.items)


class Gen(Word):
    def __init__(self, name):
        Word.__init__(self, (name,))
        self.name = name


class Relation(object):
    def __init__(self, lhs, rhs=1):
        if type(rhs) is int and rhs == 1:
            rhs = Word(())
        assert isinstance(lhs, Word), lhs
        assert isinstance(rhs, Word), rhs
        self.lhs = lhs
        self.rhs = rhs
    def __str__(self):
        return "%s == %s"%(self.lhs, self.rhs)
    @classmethod
    def promote(self, item):
        if isinstance(item, Relation):
            return item
        assert isinstance(item, Word)
        return Relation(item)
    def get_idxs(self, get, geti):
        lhs = [get[g] for g in self.lhs]
        rhs = [geti[g] for g in reversed(self.rhs)]
        return tuple(rhs + lhs)


def build(gens, rels, hrels=[]):
    for gen in gens:
        assert isinstance(gen, Gen)
    rels = [Relation.promote(item) for item in rels]
    hrels = [Relation.promote(item) for item in hrels]
    for rel in rels:
        assert isinstance(rel, Relation)
    n = len(gens)
    ngens = 2*n
    get = {gen.name : i for i, gen in enumerate(gens)}
    geti = {gen.name : i+n for i, gen in enumerate(gens)} # inverse
    relis = []
    for gen in gens:
        relis.append( (get[gen.name], geti[gen.name]) )
    relis += [rel.get_idxs(get, geti) for rel in rels]
    hrelis = [rel.get_idxs(get, geti) for rel in hrels]
    print(relis)
    graph = Schreier(ngens, relis)
    #graph.build(maxsize=4000000)
    graph.build(hrelis)
    #print(len(graph))
    return graph
    

def test():
    print("test()")

    gens = [Gen(s) for s in "w HI IH SI IS CZ".split()]
    w, HI, IH, SI, IS, CZ = gens

    rels = [
        w*w*w*w*w*w*w*w,
        w*HI == HI*w,
        w*IH == IH*w,
        w*SI == SI*w,
        w*IS == IS*w,
        w*CZ == CZ*w,
        HI*IH == IH*HI,
        HI*IS == IS*HI,
        SI*IS == IS*SI,
        SI*IH == IH*SI,
        IH*IH,
        IS*IS*IS*IS,
        IS*IH*IS*IH*IS*IH  == w,
        HI*HI,
        SI*SI*SI*SI,
        SI*HI*SI*HI*SI*HI  == w,
        CZ*CZ,
        SI*CZ == CZ*SI,
        IS*CZ == CZ*IS,
        HI*SI*SI*HI*CZ == CZ*HI*IS*SI*IS*SI*HI,
        IH*IS*IS*IH*CZ == CZ*SI*IH*SI*IS*IS*IH,
        CZ*HI*CZ == SI*HI*CZ*SI*IS*HI*SI*w*w*w*w*w*w*w,
        CZ*IH*CZ == IS*IH*CZ*IS*SI*IH*IS*w*w*w*w*w*w*w,
    ]

    ZI = SI*SI
    IZ = IS*IS
    ww = w*w
    hrels = [
        ww, ZI, IZ, HI*ZI*HI, IH*IZ*IH,
    ]

    graph = build(gens, rels) # Cliff(2)
    #graph = build(gens, rels+[w]) # ASp(4)
    #graph = build(gens, rels+hrels) # Mp(4)
    #graph = build(gens, rels+hrels+[w]) # Sp(4)
    print(len(graph))

    perms = graph.get_gens()
    print(len(perms))

    f = open("Cliff2_gens.py", "w")
    for g in perms:
        print(g, file=f)
    f.close()

    
    

if __name__ == "__main__":
    test()

    print("OK\n\n")

