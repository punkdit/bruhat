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

# https://docs.sympy.org/latest/modules/solvers/solvers.html
from sympy.solvers import solve
from sympy.solvers import nonlinsolve
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
i = j = 1 # specialize ...
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
    print("make_eqs:")
    m, n = lhs.shape
    assert rhs.shape == (m, n), rhs.shape
    for i in range(m):
      for j in range(n):
        eq = lhs[i,j] - rhs[i,j]
        #print(eq)
        if eq != 0:
            print("\t", eq)
            eqs.append(eq)
    print(len(eqs), "eqs")

def test_eq(lhs, rhs):
    m, n = lhs.shape
    assert rhs.shape == (m, n), rhs.shape
    for i in range(m):
      for j in range(n):
        eq = lhs[i,j] - rhs[i,j]
        if hasattr(eq, "simplify"):
            eq = eq.simplify()
        if eq != 0:
            return False
    return True

# unit 
make_eqs(compose(D, IE), I)
make_eqs(compose(D, EI), I)

# _assoc
make_eqs(compose(D, tensor(D, I)), compose(D, tensor(I, D)))

# cocommutative
make_eqs(D, compose(D, SWAP))


result = solve(eqs, a, b, c, d, e, f, g, h, dict=True)

assert len(result) == 1
result = result[0]

print(result)
keys = list(result.keys())
keys.sort(key=str)
for k in keys:
    print("\t", k, "=", result[k])

for eq in eqs:
    val = eq.subs(result)
    val = val.simplify()
    assert val == 0

# -------------------------------------------------------------------------------

if 1:

    # choose something generic 
    g = 7
    h = 22

    # is this equivalent to the standard D0 under basis change? (answer: no)
    D = array([
        [g+1, h-1],
        [-g, 1-h],
        [-g, 1-h],
        [g, h],
    ])
    
    D0 = array([
        [1, 0],
        [0, 0],
        [0, 0],
        [0, 1],
    ])

    # unit 
    assert numpy.alltrue(compose(D, IE) == I)
    assert numpy.alltrue(compose(D, EI) == I)
    
    # _assoc
    assert test_eq(compose(D, tensor(D, I)), compose(D, tensor(I, D)))
    
    # cocommutative
    assert numpy.alltrue(compose(D, SWAP) == D)
    
    eqs = []
    
    A = numpy.array([[a, b], [c, d]])
    
    ai, bi, ci, di = [Symbol(ch) for ch in 'ai bi ci di'.split()]
    Ai = numpy.array([[ai, bi], [ci, di]])

    AA = tensor(A, A)

    lhs = compose(D0, AA)
    rhs = compose(Ai, D)

    make_eqs(compose(A, Ai), I)
    make_eqs(lhs, rhs)

    result = solve(eqs, dict=True)
    
    if len(result) == 0:
        print("no solution")

    assert len(result) == 1, result
    result = result[0]
    
    print(result)
    keys = list(result.keys())
    keys.sort(key=str)
    for k in keys:
        print("\t", k, "=", result[k])
    
    for eq in eqs:
        val = eq.subs(result)
        val = val.simplify()
        assert val == 0


