#!/usr/bin/env python

from bruhat.element import GenericElement, Ring

class NumberRing(Ring):
    def __init__(self, reduct=(-1, 0)):
        self.one = Number((1, 0), self)
        self.zero = Number((0, 1), self)
        self.i = Number((0, 1), self)
        self.reduct = reduct

    def promote(self, a):
        if isinstance(a, Number):
            assert a.tp == self
            return a
        if isinstance(a, int):
            return Number((a, 0), self)
        try:
            a = int(a)
            return Number((a, 0), self)
        except ValueError:
            pass
        except TypeError:
            pass
        return None

    def add(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]+b[0], a[1]+b[1]), self)

    def sub(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]-b[0], a[1]-b[1]), self)

    def mul(self, a, b):
        a, b = a.value, b.value
        reduct = self.reduct
        c = a[1]*b[1]
        re = a[0]*b[0] + reduct[0]*c
        im = a[0]*b[1] + a[1]*b[0] + reduct[1]*c
        return Number((re, im), self)

    def neg(self, a):
        a = a.value
        return Number((-a[0], -a[1]), self)

    def mod(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]%b[0], a[1]%b[1]), self)



class Number(GenericElement):

    def __str__(self):
        a, b = self.value
        if b==0:
            s = str(a)
        else:
            s = "(%d+%d*i)"%(a, b)
        return s

def main():

    # Gaussian integers
    R = NumberRing((-1, 0))
    one = R.one
    i = R.i

    assert i*i == -one
    assert i**3 == -i

    assert (1+i)**2 == 2*i
    assert (2+i)*(2-i) == 5
    assert (3+2*i)*(3-2*i) == 13

    # Eisenstein integers
    R = NumberRing((-1, -1))
    one = R.one
    i = R.i

    assert i**3 == 1
    assert 1 + i + i**2 == 0

    # Kleinian integers
    R = NumberRing((-2, -1))
    one = R.one
    i = R.i

    assert i * (-1-i) == 2
    assert i * (1+i) == -2




if __name__ == "__main__":

    main()










