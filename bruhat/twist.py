#!/usr/bin/env python3


"""
"""

from random import random, randint
from functools import reduce
from operator import mul

import numpy

from bruhat.solve import array2, zeros2, dot2, shortstr, rank, find_kernel, span
from bruhat.solve import linear_independent, parse, pseudo_inverse, eq2, rand2, rank
from bruhat.action import mulclose
from bruhat.util import choose, cross, all_perms


def symplectic(n):
    F = zeros2(2*n, 2*n)
    for i in range(n):
        F[i, n+i] = 1
        F[n+i, i] = 1
    return F


def main_torus():

    n = 8


#   ZZZZZZZZ|XXXXXXXX
#   12345678|12345678
    H = parse("""
    111..1..|........
    1..1....|1.1.....
    ........|11.11...
    .1..1...|.1...1..
    ..1...1.|...1..1.
    ...11.11|........
    .....1.1|....1..1
    ........|..1..111
    """.replace("|",""))

    print()
    print("H=")
    print(shortstr(H))

    F = symplectic(n)
    C = dot2(H, dot2(F, H.transpose()))

    for i in range(n):
      for j in range(i+1, n):
        if C[i, j]:
            print("fail:", i+1, j+1)

    print(rank(H))
    H = linear_independent(H)
    print("H=")
    print(shortstr(H))

    HF = dot2(H, F)
    K = array2(find_kernel(HF))
    print("K=")
    print(shortstr(K))

    HK = numpy.concatenate((H, K))
    L = linear_independent(HK)
    print()
    print(shortstr(L))


def syparse(decl):
    for c in '0123456789':
        decl = decl.replace(c, ' ')
    decl = decl.strip().split()
    n = len(decl[0])
    m = len(decl)
    H = zeros2(m, 2*n)
    for i, row in enumerate(decl):
      for j, c in enumerate(row):
        if c=='X' or c=='Y':
            H[i, j] = 1
        if c=='Z' or c=='Y':
            H[i, n+j] = 1
    return H


def main_fail():

    print()
    print("=="*70)

    H = """
        012345678901234
     0  YXZZ...........
     1  X..X.XX.X......
     2  ZZ..ZZ.......Z.
     3  .X..Y....XZ....
     4  .ZX....Y.Z.....
     5  .....ZZZ.Z....Z
     6  ..Z...ZZ....ZZ.
     7  ...Z....Y.XX...
     8  ..XX.......ZY..
     9  .........ZXX..Y
    10  .....X.....XXXX
    11  ....X.X.X.X..X.
    """
    H = syparse(H)
    print()
    print(shortstr(H))

    m = len(H)
    n = H.shape[1] // 2
    assert n==15
    F = symplectic(n)
    
    C = dot2(H, dot2(F, H.transpose()))

    print(C.shape)

    for i in range(m):
      for j in range(i+1, m):
        if C[i, j]:
            print("fail:", i, j)

    print(rank(H))
    H = linear_independent(H)
    print("H=", H.shape)
    print(shortstr(H))

    HF = dot2(H, F)
    K = array2(find_kernel(HF))
    print("K=", K.shape)
    print(shortstr(K))

    HK = numpy.concatenate((H, K))
    L = linear_independent(HK)
    print()
    print(shortstr(L))



def main():

    # use Z3 solver to replace the I's below with X's or Z's
    H = """
    YIII................
    I..I.....III........
    II.....III..........
    .II.II......I.......
    ..IIIII.............
    .I..I..I.....I...I..
    .....IIII.....I.....
    ........IYI....I....
    ..I...I.....I......Y
    .....I....III.I.....
    ...II......I.I..I...
    ........I....IIII...
    ...............IIIY.
    ......II.........III
    ..........I..III.I..
    ...........II...I.II
    """.strip()

    #H = H.replace("Y", "I") # works !
    H = find_stabilizer_code(H)

    print(H)
    

    H = """
YIII....................................
I..III................I.................
II.............II.....I.................
....IIII........................I.......
......IIIII.............................
..II...I..I......I......................
.I.........I...I..................II....
.II........II............I..............
........I..IIII.........................
............I.III............I..........
.......I.........III............I.......
.........I........IIII..................
.....I.............I.II..........I......
..I..............I.......I..Y...........
................I....III.............I..
.............I......II.II...............
............I..........III...I..........
...II.....I...........................II
........II...I......Y...................
...................I..........I.III.....
...........I.I..........I.....I...I.....
......I.I.....I...........Y.............
..............II..........I........II...
.........II.......I........I...........I
.................II........II..I........
................I............I.......III
........................II..I.II........
.....I..........................IIII....
..............................II.I.II...
...........................I...I....III.
.......................I...I.I.......I.I
....I.I...................I.........I.I.
    """
    H = find_stabilizer_code(H)
    print(H)



def find_stabilizer_code(H):
    H = H.strip().split()
    H = [list(row) for row in H]
    H = numpy.array(H)
    m, n = H.shape
    #print(H)

    import z3
    from z3 import Bool, And, Or, Not, Implies, If, Solver

    vs = {}
    clauses = []
    solver = Solver()
    for i in range(m):
      for j in range(n):
        c = H[i, j]
        if c == '.':
            continue
        X = Bool("X_%d_%d"%(i,j))
        Z = Bool("Z_%d_%d"%(i,j))
        vs[i, j] = (X, Z)
        if c == 'I':
            # X or Z and not both
            solver.add(Or(X, Z))
            solver.add(Not(And(X, Z)))
        elif c == 'X':
            # freeze this as X
            solver.add(X)
            solver.add(Not(Z))
        elif c == 'Z':
            # freeze this Z
            solver.add(Not(X))
            solver.add(Z)
        elif c == 'Y':
            # freeze this as X and Z
            solver.add(X)
            solver.add(Z)
        else:
            assert 0, (i, j, c)

    for i in range(m):
      for j in range(i+1, m):
        bits = []
        for k in range(n):
            if H[i,k] == '.' or H[j,k] == '.':
                continue
            bits.append(k)
        if not bits:
            continue
        #print()
        #print(i, j, bits)
        for lhs in cross([(0,1)]*len(bits)):
          lvars = [vs[i, k][l] for (k,l) in zip(bits, lhs)]
          #print(lhs, lvars)
          for rhs in cross([(0,1)]*len(bits)):
            rvars = [vs[j, k][r] for (k,r) in zip(bits, rhs)]
            syndrome = len([idx for idx in range(len(bits)) if lhs[idx]!=rhs[idx]]) % 2
            #print('\t', rhs, rvars, syndrome)
            if syndrome:
                clause = Not(And(*lvars, *rvars))
                #print('\t', clause)
                solver.add(clause)
            
    result = solver.check()
    if result != z3.sat:
        return None

    H = numpy.empty((m, n), dtype=object)
    H[:] = '.'
    model = solver.model()
    for key, XZ in vs.items():
        #for v in XZ:
            #print(v, model.evaluate(v))
        X, Z = XZ
        X = model.evaluate(X)
        Z = model.evaluate(Z)
        assert X or Z
        assert not (X and Z)
        if X:
            H[key] = 'X'
        if Z:
            H[key] = 'Z'

    H = ("\n".join(''.join(row) for row in H))
    return H


if __name__ == "__main__":

    main()
    print("OK")
    print()


