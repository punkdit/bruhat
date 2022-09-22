#!/usr/bin/env python

"""
Build system of equations for finding Belyi projection
from vertex degrees of a dessin.
"""

from sage.all_cmdline import RationalField, PolynomialRing
from sage.all_cmdline import *


def test():
    x = PolynomialRing(RationalField(), 'x').gen()
    f = (x**3 - 1)**2-(x**2-1)**2
    print( f.factor() )


def main():

    pass



if __name__ == "__main__":

    main()




