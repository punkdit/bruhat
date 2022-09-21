#!/usr/bin/env python

from sage.all_cmdline import RationalField, PolynomialRing
from sage.all_cmdline import *


x = PolynomialRing(RationalField(), 'x').gen()
f = (x**3 - 1)**2-(x**2-1)**2
print( f.factor() )




