#!/usr/bin/env python3

"""
Here we make some Categories where the objects are
natural numbers: 0, 1, 2, ...

"""

from bruhat.argv import argv

# -------------------------------------------


class Cat(object):
    pass


class Op(object):
    def __init__(self, cat, tgt, src):
        assert isinstance(cat, Cat)
        self.cat = cat
        self.tgt = int(tgt)
        self.src = int(src)

    def __mul__(self, other):
        assert self.cat is other.cat
        assert self.src == other.tgt
        return self.cat.mul(self, other)

    def __str__(self):
        return "(%d <-- %d)"%(self.tgt, self.src)


# -------------------------------------------


class FreeCat(Cat):

    def __init__(self):
        self.cache = {}

    def mul(self, *ops):
        for op in ops:
            assert op.cat is self
        cache = self.cache
        op = cache.get(ops)
        if op is None:
            op = Word(self, ops)
            cache[ops] = op
        return op

    def gen(self, tgt, src):
        return Gen(self, tgt, src)



class Gen(Op):

#    def __init__(self, cat, name):
#        Op.__init__(self, cat)
#        self.name = name
#
#    def __str__(self):
#        return self.name
#
#    def __mul__(self, other):
#        assert self.src == other.tgt
#        return Word(self.cat, [self, other])
    pass
        

class Word(Op):
    def __init__(self, cat, gens):
        items = []
        for g in gens:
            items += g.items if isinstance(g, Word) else [g]
        self.items = items
        Op.__init__(self, cat, items[0].tgt, items[-1].src)

    def __str__(self):
        return "*".join(str(g) for g in self.items)



def test():

    cat = FreeCat()
    gen = cat.gen
    g = gen(1, 1)
    h = gen(1, 1)

    assert g != h
    assert g*h is g*h
    assert g*h == g*h
    assert g*h != h*g
    

# -------------------------------------------



if __name__ == "__main__":

    from time import sleep, time
    start_time = time()
    profile = argv.profile
    name = argv.next() or "test"
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    else:
        fn = eval(name)
        fn()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))


