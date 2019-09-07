#!/usr/bin/env sage

from __future__ import print_function

from sage import *
from sage.rings.number_field.number_field import NumberField

#import sys
#print(sys.argv)


if 1:
    K = NumberField(x^4+1)
    
    r2 = e8+e8^7
    assert r2**2 == 2
    
    G = K.galois_group()
    assert len(G) == 4
