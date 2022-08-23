#!/usr/bin/env python

# https://ncatlab.org/nlab/show/globular+set

bars = ["-", "=", "≡", "≣"] # ...?

class Cell(object):
    def __init__(self, *args):
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
    print(u)



if __name__ == "__main__":

    main()




