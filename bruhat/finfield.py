#!/usr/bin/env python3


"""

"""

import numpy
from numpy import dot, kron

from bruhat.argv import argv
from bruhat.util import cross


scalar = numpy.int64 # ?


class Matrix(object):
    def __init__(self, p, a):
        a = numpy.array(a, dtype=scalar)
        a %= p
        self.a = a
        self.shape = a.shape

    @classmethod
    def zeros(self, p, shape):
        a = numpy.zeros(shape, dtype=scalar)
        return Matrix(p, a)

    def __add__(self, other):
        assert self.p == other.p
        a = self.a + other.a
        return Matrix(self.p, a)

    def __mul__(self, other):
        assert self.p == other.p
        a = dot(self.a, other.a)
        return Matrix(self.p, a)

    def __matmul__(self, other):
        assert self.p == other.p
        a = kron(self.a, other.a)
        return Matrix(self.p, a)




def main():

    dim = argv.get("dim", 3)
    p = argv.get("p", 3)

    n = dim**3
    #for items in cross([tuple(range(p))]*n): # hopeless...
    #    print(items)


    


if __name__ == "__main__":

    main()



