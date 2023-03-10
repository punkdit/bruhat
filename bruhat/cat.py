#!/usr/bin/env python

from time import time
start_time = time()

from bruhat.argv import argv



class Category(object):

    def __init__(self, name="C", dim=0):
        self.name = name
        self.dim = dim
        self.cells = [set() for i in range(dim+1)]
        self._mul = {} 

    def __str__(self):
        return self.name

    def Cell(self, tgt=None, src=None, codim=0, name="?"):
        cell = Cell(self, tgt, src, codim, name)
        return cell

    def free_mul(self, lhs, rhs):
        key = (lhs, rhs)
        cell = self._mul.get(key)
        if cell is None:
            cell = self.Cell(lhs.tgt, rhs.src, lhs.codim)
            self._mul[key] = cell
        return cell

    def set_mul(self, lhs, rhs, value):
        assert isinstance(lhs, Cell)
        assert isinstance(rhs, Cell)
        assert isinstance(value, Cell)
        assert lhs.tgt == value.tgt
        assert lhs.src == rhs.tgt
        assert rhs.src == value.src
        key = (lhs, rhs)
        assert self._mul.get(key) is None
        self._mul[key] = value

    def check(self):
        cells = self.cells
        if self.dim == 0:
            return
        pairs = []
        for a in cells[0]:
          for b in cells[0]:
            if b.tgt == a.src:
                pairs.append((a, b))
                assert (a, b) in self._mul
        for lhs in pairs:
          for rhs in pairs:
            if lhs[1] == rhs[0]:
                a, b = lhs
                c = rhs[1]
                assert (a*b) * c == a * (b*c)


class Cell(object):
    def __init__(self, cat, tgt=None, src=None, codim=0, name="?"):
        assert isinstance(cat, Category)
        assert tgt is None or isinstance(tgt, Cell)
        assert src is None or isinstance(src, Cell)
        assert 0<=codim<=cat.dim
        assert tgt is None or tgt.codim == codim+1
        assert src is None or src.codim == codim+1
        assert (tgt is None) == (src is None)
        self.cat = cat
        self.tgt = tgt
        self.src = src
        self.codim = codim
        self.name = name
        self.is_object = tgt is None
        cat.cells[codim].add(self)

    def __str__(self):
        if self.is_object:
            return self.name
        else:
            return "%s <-- %s" % (self.tgt, self.src)

    def __mul__(lhs, rhs):
        assert lhs.cat == rhs.cat
        assert rhs.tgt == lhs.src
        return lhs.cat.free_mul(lhs, rhs)


def test():
    S = Category("S", 0) # a set
    X = S.Cell()
    Y = S.Cell()

    C = Category("C", 1) # a Category
    Cell = C.Cell
    A = Cell(codim=1, name="A")
    B = Cell(codim=1, name="B")
    f = Cell(B, A)
    Ai = Cell(A, A)
    Bi = Cell(B, B)
    C.set_mul(f, Ai, f)
    C.set_mul(Bi, f, f)
    C.set_mul(Ai, Ai, Ai)
    C.set_mul(Bi, Bi, Bi)

    assert f*Ai == f*Ai
    assert f*Ai == f

    C.check()


if __name__ == "__main__":

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

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()



