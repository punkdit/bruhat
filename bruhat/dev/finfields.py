#!/usr/bin/env python3

from bruhat.argv import argv
from bruhat.util import cross
from bruhat.element import FiniteField, PolynomialRing


def main():

    p = argv.get("p", 3)

    field = FiniteField(p)
    ring = PolynomialRing(field)
    x = ring.x

    found = set()
    items = [list(field.elements)]*p
    for cs in cross(items):
        f = ring.zero
        for c in cs:
            f = x*f
            f = f+c
        vals = tuple(f(j) for j in field.elements)
        found.add(vals)
        print(' '.join([str(v) for v in vals]))
    print(len(found), p**p)


if __name__ == "__main__":

    main()


