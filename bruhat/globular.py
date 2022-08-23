#!/usr/bin/env python

# https://ncatlab.org/nlab/show/globular+set

bars = ["-", "=", "≡", "≣"] # ...?

class Cell(object):
    def __init__(self, *args):
        """
        call via either:
            Cell(tgt, src)
            Cell(name)
        """
        tgt = None
        src = None
        if len(args)==1:
            name = args[0]
            n = 0
        else:
            assert len(args)==2
            tgt, src = args
            assert isinstance(src, Cell)
            assert isinstance(tgt, Cell)
            n = src.n+1
            assert src.n == tgt.n
            name = "(%s <%s %s)"%(tgt.name, bars[n-1], src.name)
        if n>1:
            # check globular conditions
            assert src.src == tgt.src
            assert src.tgt == tgt.tgt
        self.n = n
        self.src = src
        self.tgt = tgt
        self.name = name

    def __str__(self):
        return self.name

    def __repr__(self):
        if self.tgt is not None:
            return "Cell(%s, %s)"%(self.tgt, self.src)
        else:
            return "Cell(%r)"%(self.name,)

#    @classmethod
#    def whisker(cls, 

    def __mul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        assert n>=2
        return Pair(n-2, lhs, rhs)

    def __lshift__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        assert n>=1
        return Pair(n-1, lhs, rhs)

    def __matmul__(lhs, rhs):
        assert lhs.n == rhs.n 
        n = lhs.n
        assert n>=1
        return Pair(n, lhs, rhs)


class Pair(Cell):

    @classmethod
    def whisker(cls, lhs, rhs):
        while lhs.n < rhs.n:
            lhs = Cell(lhs, lhs)
        while rhs.n < lhs.n:
            rhs = Cell(rhs, rhs)
        assert lhs.n == rhs.n
        return (lhs, rhs)

    cache = {}
    def __new__(cls, codim, lhs, rhs):
        assert isinstance(lhs, Cell)
        assert isinstance(rhs, Cell)
        key = (codim, lhs, rhs)
        if key in cls.cache:
            return cls.cache[key]
        ob = object.__new__(cls)
        cls.cache[key] = ob
        assert codim >= 0
        assert lhs.n == rhs.n, "use whisker first"
        if codim == 0:
            assert rhs.tgt == lhs.src, "not composable"
            tgt = lhs.tgt
            src = rhs.src
        elif codim == 1:
            assert lhs.src.src == rhs.src.tgt
            # lhs.tgt.src == rhs.tgt.tgt
            tgt = Pair(0, lhs.tgt, rhs.tgt)
            src = Pair(0, lhs.src, rhs.src)
        elif codim == 2:
            assert lhs.src.src.src == rhs.src.src.tgt
            # -> lhs.tgt.src.src == rhs.tgt.src.tgt
            # -> lhs.tgt.tgt.src == rhs.tgt.tgt.tgt
            tgt = Pair(1, lhs.tgt, rhs.tgt)
            src = Pair(1, lhs.src, rhs.src)
        else:
            tgt = Pair(codim-1, lhs.tgt, rhs.tgt)
            src = Pair(codim-1, lhs.src, rhs.src)
        Cell.__init__(ob, tgt, src)
        ob.codim = codim
        ob.lhs = lhs
        ob.rhs = rhs
        return ob

    def __init__(self, codim, lhs, rhs):
        pass

    def __str__(self):
        return "(%s *%d %s)" % (self.lhs, self.codim, self.rhs)



def main():

    a, b, c = [Cell(ch) for ch in 'abc']
    assert str(a) == 'a'

    f = Cell(a, b)
    g = Cell(a, b)
    u = Cell(f, g)
    assert str(u) == "((a <- b) <= (a <- b))"

    assert u==u
    assert u != Cell(f, g)

    u = Cell(u, u)
    u = Cell(u, u)

    g = Cell(c, b)
    f = Cell(b, a)
    uu = Pair(0, g, f)

    assert uu is Pair(0, g, f)
    assert uu == Pair(0, g, f)





if __name__ == "__main__":

    main()




