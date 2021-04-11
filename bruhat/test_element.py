#!/usr/bin/env python3

from _element import Fraction


a = Fraction(2,3)
b = Fraction(3,4)

print(a, b)

assert a==a
assert a!=b
assert a!=2
assert Fraction(4,1) == 4
assert a<b<2
assert a<=b<=2
assert not b<a

assert Fraction(4,6) == a

assert 3*a==2
assert a+2 == 2+a == Fraction(8,3)
assert a-2 == Fraction(-4,3)
assert 2-a == Fraction(4,3)

assert a+a+a == 2
assert a**2 == Fraction(4, 9)

assert a//2 == Fraction(1,3)
assert a/2 == Fraction(1,3)
assert 2//a == 3



print("OK")




