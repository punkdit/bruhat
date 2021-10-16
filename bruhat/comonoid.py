#!/usr/bin/env python3
"""
Symbolically find a cocommutative comonoid 
on a 2d vector space.
"""

import numpy

def compose(first, second):
    return numpy.dot(second, first)
tensor = numpy.kron
dot = numpy.dot
array = numpy.array


from sympy.solvers import solve
from sympy import Symbol

I = array([[1, 0], [0, 1]], dtype=object)
SWAP = array([
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
], dtype=object)

#x = Symbol('x')
syms = [Symbol(ch) for ch in 'abcdefghij']
a, b, c, d, e, f, g, h, i, j = syms
i = j = 1
D = array([
    [a, b],
    [c, d],
    [e, f],
    [g, h],
])
E = array([[i, j]])

IE = tensor(I, E)
EI = tensor(E, I)

eqs = []
def make_eqs(lhs, rhs):
    m, n = lhs.shape
    assert rhs.shape == (m, n), rhs.shape
    for i in range(m):
      for j in range(n):
        eq = lhs[i,j] - rhs[i,j]
        print(eq)
        if eq != 0:
            eqs.append(eq)


# unit 
make_eqs(compose(D, IE), I)
make_eqs(compose(D, EI), I)

# assoc
make_eqs(compose(D, tensor(D, I)), compose(D, tensor(I, D)))

# cocommutative
make_eqs(D, compose(D, SWAP))

print(len(eqs))
print(len(set(eqs)))

result = solve(eqs, syms, dict=True)

print(result)
for row in result:
    print(row)
    keys = list(row.keys())
    keys.sort(key=str)
    for k in keys:
        print("\t", k, "=", row[k])

