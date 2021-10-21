#!/usr/bin/env python3
"""
Symbolically find a cocommutative comonoid 
on a 2d vector space.
"""

import numpy
from numpy import dot, array, empty, alltrue

def compose(first, second):
    return numpy.dot(second, first)
tensor = numpy.kron

# https://docs.sympy.org/latest/modules/solvers/solvers.html
from sympy.solvers import solve
from sympy.solvers import nonlinsolve
from sympy import Symbol

from bruhat.argv import argv


class System(object):
    def __init__(self):
        self.eqs = []
        self.idx = 0
        self.vs = []
        self.items = [] # list of array's

    def get_var(self, stem='v'):
        ch = "%s_%d"%(stem, self.idx)
        self.idx += 1
        self.vs.append(ch)
        return Symbol(ch)

    def array(self, rows, cols, name='v'):
        rows = [[self.get_var(name) for j in range(cols)] for i in range(rows)]
        A = numpy.array(rows)
        self.items.append(A)
        return A

    def subs(self, values):
        items = [A.copy() for A in self.items]
        for i, A in enumerate(items):
            for idx in numpy.ndindex(A.shape):
                A[idx] = A[idx].subs(values)
        return items

    def add(self, lhs, rhs):
        #print("add:")
        eqs = self.eqs
        m, n = lhs.shape
        assert rhs.shape == (m, n), rhs.shape
        for i in range(m):
          for j in range(n):
            eq = lhs[i,j] - rhs[i,j]
            #print(eq)
            if eq != 0:
                #print("\tadd:", eq)
                eqs.append(eq)
        #print(len(eqs), "eqs")

    def sage_dump(self):
        vs = ','.join(self.vs)
        print("vs = %s = var('%s')" % (vs, vs))
        print("eqs = [")
        for eq in self.eqs:
            print('    %s == 0,' % eq)
        print("]")
        #print("slns = solve(eqs, vs, solution_dict=True)")
        print("slns = solve(eqs, vs)")
    
    def test_eq(self, lhs, rhs):
        m, n = lhs.shape
        assert rhs.shape == (m, n), rhs.shape
        for i in range(m):
          for j in range(n):
            eq = lhs[i,j] - rhs[i,j]
            eq = eq.simplify()
            if eq != 0:
                return False
        return True

    def solve(self):
    
        eqs = self.eqs
        result = solve(eqs, dict=True)
        
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
        return result

        
dim = argv.get("dim", 2)

dim2 = dim**2

I = empty((dim, dim), dtype=object)
I[:] = 0
for i in range(dim):
    I[i, i] = 1
#print(I)

SWAP = empty((dim, dim, dim, dim), dtype=object)
SWAP[:] = 0
for i in range(dim):
  for j in range(dim):
    SWAP[i, j, j, i] = 1
SWAP.shape = dim2, dim2
#print(SWAP)

if dim==2:
    lhs = array([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
    ], dtype=object)
    assert alltrue(lhs==SWAP)

assert alltrue(dot(SWAP, SWAP)==tensor(I, I))

if 0:
    system = System()
    #M = system.array(2, 4, "M")
    M = array([
        [1, 0, 0, 0],
        [0, 0, 0, 1],
    ], dtype=object)
    G = system.array(2, 2, "G")
    Gi = system.array(2, 2, "Gi")
    
    GG = tensor(G, G)
    
    system.add(compose(G, Gi), I)
    system.eqs.append(G[0,0]*G[1,1] - G[0,1]*G[1,0]-1) # SL(2)
    soln = system.solve()
    
    M1 = compose(compose(GG, M), Gi)
    print(M1)
    for i in range(2):
      for j in range(4):
        M1[i, j] = M1[i, j].subs(soln)
        M1[i, j].simplify()
    print(M1)


if argv.frobenius:
    # Commutative special Frobenius algebra
    
    system = System()
    
    F = system.array(dim, dim**2, "F") # mul
    G = system.array(dim, 1, "G") # unit
    
    D = system.array(dim**2, dim, "D") # comul
    E = system.array(1, dim, "E") # counit
    
    IF = tensor(I, F)
    FI = tensor(F, I)
    
    IG = tensor(I, G)
    GI = tensor(G, I)
    
    ID = tensor(I, D)
    DI = tensor(D, I)
    
    IE = tensor(I, E)
    EI = tensor(E, I)
    
    # unit
    system.add(compose(IG, F), I)
    system.add(compose(GI, F), I)
    
    # _assoc
    system.add(compose(FI, F), compose(IF, F))
    
    # commutative
    #system.add(F, compose(SWAP, F))
    
    # counit 
    system.add(compose(D, IE), I)
    system.add(compose(D, EI), I)
    
    # _coassoc
    system.add(compose(D, DI), compose(D, ID))
    
    # cocommutative
    system.add(D, compose(D, SWAP))
    
    # Frobenius
    system.add(compose(DI, IF), compose(F, D))
    system.add(compose(ID, FI), compose(F, D))
    
    # special
    system.add(compose(D, F), I)
    
    system.sage_dump()

    if 0:
        print(len(system.eqs))
        system.solve()
        print(F)
        print(G)
        print(D)
        print(E)

    exit()

    """
    F_0, F_1, F_2, F_3, F_4, F_5, F_6, F_7, G_8, G_9, 
    D_10, D_11, D_12, D_13, D_14, D_15, D_16, D_17, E_18, E_19

    [[F_0 F_1 F_2 F_3]
     [F_4 F_5 F_6 F_7]]
    [[G_8]
     [G_9]]
    [[D_10 D_11]
     [D_12 D_13]
     [D_14 D_15]
     [D_16 D_17]]
    [[E_18 E_19]]
    """

    vs = "F_0 F_1 F_2 F_3 F_4 F_5 F_6 F_7 G_8 G_9 " 
    vs += "D_10 D_11 D_12 D_13 D_14 D_15 D_16 D_17 E_18 E_19"
    vs = vs.split()
    for v in vs:
        exec("%s = 1"%v)

    # hmmm.... this is all useless....

    # result from sage_dump:
    [
        D_12 - D_14 == 0,
        -D_12 + D_14 == 0,
        D_11*D_12 - D_11*D_14 == 0,
        -D_11*D_12 + D_11*D_14 == 0,
        D_13 - D_15 == 0,
        -D_13 + D_15 == 0,
        -D_11*D_13 + D_11*D_15 == 0,
        D_10*D_12 - D_10*D_14 + D_13*D_14 - D_12*D_15 == 0,
        D_12*D_13 - D_11*D_16 == 0,
        D_14*D_15 - D_11*D_16 == 0,
        -D_12*D_13 + D_11*D_16 == 0,
        -D_14*D_15 + D_11*D_16 == 0,
        D_12*D_16 - D_14*D_16 == 0,
        D_13*D_16 - D_15*D_16 == 0,
        -D_13*D_16 + D_15*D_16 == 0,
        D_11*D_14 - D_10*D_15 + D_15**2 - D_11*D_17 == 0,
        -D_11*D_12 + D_10*D_13 - D_13**2 + D_11*D_17 == 0,
        D_12**2 - D_10*D_16 + D_13*D_16 - D_12*D_17 == 0,
        -D_14**2 + D_10*D_16 - D_15*D_16 + D_14*D_17 == 0,
        D_13*D_14 - D_12*D_15 - D_13*D_17 + D_15*D_17 == 0,
        D_10*E_18 + D_12*E_19 - 1 == 0,
        D_11*E_18 + D_13*E_19 == 0,
        D_10*E_18 + D_14*E_19 - 1 == 0,
        D_11*E_18 + D_15*E_19 == 0,
        D_12*E_18 + D_16*E_19 == 0,
        D_14*E_18 + D_16*E_19 == 0,
        D_13*E_18 + D_17*E_19 - 1 == 0,
        D_15*E_18 + D_17*E_19 - 1 == 0,
        D_10*F_0 + D_12*F_1 + D_14*F_2 + D_16*F_3 - 1 == 0,
        D_11*F_0 + D_13*F_1 + D_15*F_2 + D_17*F_3 == 0,
        F_1*F_3 - F_2*F_3 == 0,
        D_14*F_1 - D_11*F_4 == 0,
        D_12*F_2 - D_11*F_4 == 0,
        -D_14*F_1 + D_11*F_4 == 0,
        -D_12*F_2 + D_11*F_4 == 0,
        D_16*F_1 - D_13*F_4 == 0,
        -D_16*F_1 + D_13*F_4 == 0,
        D_16*F_2 - D_15*F_4 == 0,
        -D_16*F_2 + D_15*F_4 == 0,
        F_1*F_4 - F_2*F_4 == 0,
        -F_1*F_4 + F_2*F_4 == 0,
        D_11*F_0 - D_10*F_1 + D_15*F_1 - D_11*F_5 == 0,
        D_12*F_3 - D_11*F_5 == 0,
        -D_12*F_3 + D_11*F_5 == 0,
        D_13*F_0 - D_12*F_1 + D_17*F_1 - D_13*F_5 == 0,
        -D_14*F_0 + D_10*F_4 - D_15*F_4 + D_14*F_5 == 0,
        D_16*F_3 - D_15*F_5 == 0,
        -D_16*F_3 + D_15*F_5 == 0,
        -D_16*F_0 + D_12*F_4 - D_17*F_4 + D_16*F_5 == 0,
        F_3*F_4 - F_1*F_5 == 0,
        -F_3*F_4 + F_1*F_5 == 0,
        D_11*F_0 - D_10*F_2 + D_13*F_2 - D_11*F_6 == 0,
        D_14*F_3 - D_11*F_6 == 0,
        -D_14*F_3 + D_11*F_6 == 0,
        -D_12*F_0 + D_10*F_4 - D_13*F_4 + D_12*F_6 == 0,
        D_16*F_3 - D_13*F_6 == 0,
        -D_16*F_3 + D_13*F_6 == 0,
        D_15*F_0 - D_14*F_2 + D_17*F_2 - D_15*F_6 == 0,
        -D_16*F_0 + D_14*F_4 - D_17*F_4 + D_16*F_6 == 0,
        F_0*F_1 - F_0*F_2 + F_2*F_5 - F_1*F_6 == 0,
        F_3*F_4 - F_2*F_6 == 0,
        -F_3*F_4 + F_2*F_6 == 0,
        F_3*F_5 - F_3*F_6 == 0,
        -F_3*F_5 + F_3*F_6 == 0,
        -F_4*F_5 + F_4*F_6 == 0,
        D_11*F_1 - D_10*F_3 + D_13*F_3 - D_11*F_7 == 0,
        D_11*F_2 - D_10*F_3 + D_15*F_3 - D_11*F_7 == 0,
        -D_12*F_1 + D_10*F_5 - D_13*F_5 + D_12*F_7 == 0,
        D_13*F_2 - D_12*F_3 + D_17*F_3 - D_13*F_7 == 0,
        -D_14*F_2 + D_10*F_6 - D_15*F_6 + D_14*F_7 == 0,
        D_15*F_1 - D_14*F_3 + D_17*F_3 - D_15*F_7 == 0,
        -D_16*F_1 + D_14*F_5 - D_17*F_5 + D_16*F_7 == 0,
        D_10*F_4 + D_12*F_5 + D_14*F_6 + D_16*F_7 == 0,
        -D_16*F_2 + D_12*F_6 - D_17*F_6 + D_16*F_7 == 0,
        D_11*F_4 + D_13*F_5 + D_15*F_6 + D_17*F_7 - 1 == 0,
        F_1**2 - F_0*F_3 + F_3*F_5 - F_1*F_7 == 0,
        -F_2**2 + F_0*F_3 - F_3*F_6 + F_2*F_7 == 0,
        F_2*F_4 - F_0*F_6 + F_6**2 - F_4*F_7 == 0,
        -F_1*F_4 + F_0*F_5 - F_5**2 + F_4*F_7 == 0,
        F_2*F_5 - F_1*F_6 - F_5*F_7 + F_6*F_7 == 0,
        F_0*G_8 + F_1*G_9 - 1 == 0,
        F_0*G_8 + F_2*G_9 - 1 == 0,
        F_1*G_8 + F_3*G_9 == 0,
        F_2*G_8 + F_3*G_9 == 0,
        F_4*G_8 + F_5*G_9 == 0,
        F_4*G_8 + F_6*G_9 == 0,
        F_5*G_8 + F_7*G_9 - 1 == 0,
        F_6*G_8 + F_7*G_9 - 1 == 0,
    ]



    # result from sage_dump:
    soln = [
        D_12 - D_14, -D_12 + D_14,
        D_11*D_12 - D_11*D_14, -D_11*D_12 + D_11*D_14,
        D_13 - D_15, -D_13 + D_15,
        -D_11*D_13 + D_11*D_15, D_10*D_12 - D_10*D_14 + D_13*D_14 - D_12*D_15,
        D_12*D_13 - D_11*D_16, D_14*D_15 - D_11*D_16,
        -D_12*D_13 + D_11*D_16, -D_14*D_15 + D_11*D_16,
        D_12*D_16 - D_14*D_16, D_13*D_16 - D_15*D_16,
        -D_13*D_16 + D_15*D_16, D_11*D_14 - D_10*D_15 + D_15**2 - D_11*D_17,
        -D_11*D_12 + D_10*D_13 - D_13**2 + D_11*D_17,
        D_12**2 - D_10*D_16 + D_13*D_16 - D_12*D_17,
        -D_14**2 + D_10*D_16 - D_15*D_16 + D_14*D_17,
        D_13*D_14 - D_12*D_15 - D_13*D_17 + D_15*D_17,
        D_10*E_18 + D_12*E_19 - 1, D_11*E_18 + D_13*E_19,
        D_10*E_18 + D_14*E_19 - 1, D_11*E_18 + D_15*E_19,
        D_12*E_18 + D_16*E_19, D_14*E_18 + D_16*E_19,
        D_13*E_18 + D_17*E_19 - 1, D_15*E_18 + D_17*E_19 - 1,
        D_10*F_0 + D_12*F_1 + D_14*F_2 + D_16*F_3 - 1,
        D_11*F_0 + D_13*F_1 + D_15*F_2 + D_17*F_3,
        F_1*F_3 - F_2*F_3, D_14*F_1 - D_11*F_4,
        D_12*F_2 - D_11*F_4, -D_14*F_1 + D_11*F_4,
        -D_12*F_2 + D_11*F_4, D_16*F_1 - D_13*F_4,
        -D_16*F_1 + D_13*F_4, D_16*F_2 - D_15*F_4,
        -D_16*F_2 + D_15*F_4, F_1*F_4 - F_2*F_4,
        -F_1*F_4 + F_2*F_4, D_11*F_0 - D_10*F_1 + D_15*F_1 - D_11*F_5,
        D_12*F_3 - D_11*F_5, -D_12*F_3 + D_11*F_5,
        D_13*F_0 - D_12*F_1 + D_17*F_1 - D_13*F_5,
        -D_14*F_0 + D_10*F_4 - D_15*F_4 + D_14*F_5,
        D_16*F_3 - D_15*F_5, -D_16*F_3 + D_15*F_5,
        -D_16*F_0 + D_12*F_4 - D_17*F_4 + D_16*F_5,
        F_3*F_4 - F_1*F_5, -F_3*F_4 + F_1*F_5,
        D_11*F_0 - D_10*F_2 + D_13*F_2 - D_11*F_6,
        D_14*F_3 - D_11*F_6, -D_14*F_3 + D_11*F_6,
        -D_12*F_0 + D_10*F_4 - D_13*F_4 + D_12*F_6,
        D_16*F_3 - D_13*F_6, -D_16*F_3 + D_13*F_6,
        D_15*F_0 - D_14*F_2 + D_17*F_2 - D_15*F_6,
        -D_16*F_0 + D_14*F_4 - D_17*F_4 + D_16*F_6,
        F_0*F_1 - F_0*F_2 + F_2*F_5 - F_1*F_6, F_3*F_4 - F_2*F_6,
        -F_3*F_4 + F_2*F_6, F_3*F_5 - F_3*F_6,
        -F_3*F_5 + F_3*F_6, -F_4*F_5 + F_4*F_6,
        D_11*F_1 - D_10*F_3 + D_13*F_3 - D_11*F_7,
        D_11*F_2 - D_10*F_3 + D_15*F_3 - D_11*F_7,
        -D_12*F_1 + D_10*F_5 - D_13*F_5 + D_12*F_7,
        D_13*F_2 - D_12*F_3 + D_17*F_3 - D_13*F_7,
        -D_14*F_2 + D_10*F_6 - D_15*F_6 + D_14*F_7,
        D_15*F_1 - D_14*F_3 + D_17*F_3 - D_15*F_7,
        -D_16*F_1 + D_14*F_5 - D_17*F_5 + D_16*F_7,
        D_10*F_4 + D_12*F_5 + D_14*F_6 + D_16*F_7,
        -D_16*F_2 + D_12*F_6 - D_17*F_6 + D_16*F_7,
        D_11*F_4 + D_13*F_5 + D_15*F_6 + D_17*F_7 - 1,
        F_1**2 - F_0*F_3 + F_3*F_5 - F_1*F_7,
        -F_2**2 + F_0*F_3 - F_3*F_6 + F_2*F_7,
        F_2*F_4 - F_0*F_6 + F_6**2 - F_4*F_7,
        -F_1*F_4 + F_0*F_5 - F_5**2 + F_4*F_7,
        F_2*F_5 - F_1*F_6 - F_5*F_7 + F_6*F_7,
        F_0*G_8 + F_1*G_9 - 1, F_0*G_8 + F_2*G_9 - 1,
        F_1*G_8 + F_3*G_9, F_2*G_8 + F_3*G_9,
        F_4*G_8 + F_5*G_9, F_4*G_8 + F_6*G_9,
        F_5*G_8 + F_7*G_9 - 1, F_6*G_8 + F_7*G_9 - 1]


    print(len(soln))



