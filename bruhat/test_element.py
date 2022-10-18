#!/usr/bin/env python3

from _element import Fraction, FiniteField


a = Fraction(2,3)
b = Fraction(3,4)

assert str(a) == "2/3"
assert str(b) == "3/4"

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


try:
    for _ in range(100):
        a = a + a
    assert 0
except OverflowError:
    pass

a = Fraction(2,1)
try:
    for _ in range(100):
        a = a * a
    assert 0
except OverflowError:
    pass

try:
    b = a ** 100
except OverflowError:
    b = None
assert b is None


ring = FiniteField(3)
zero = ring.zero
one = ring.one
assert zero + zero == zero
assert zero + one == one
#assert one + one == 3*one
assert 1 + one == one+one
assert one + 1 == one+one
assert 2*one == one+one
assert one*2 == one+one
assert one*one == one
two = one+one
assert two*two == one
assert two+two == one

assert two//two == one
assert two/two == one

print("OK")




