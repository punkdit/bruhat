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

    def show(self, x):
        assert len(x) == len(self.vs)
        ns = dict((str(v), xi) for (v, xi) in zip(self.vs, x))
        for item in self.subs(ns):
            print(item)

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

    def __call__(self, x):
        assert len(x) == len(self.vs)
        eqs = self.eqs
        ns = dict((str(v), xi) for (v, xi) in zip(self.vs, x))
        value = []
        for eq in eqs:
            value.append(eq.subs(ns))
        return value

    def sage_dump(self):
        vs = ','.join(self.vs)
        print("vs = %s = var('%s')" % (vs, vs))
        print("eqs = [")
        for eq in self.eqs:
            print('    %s == 0,' % eq)
        print("]")
        #print("slns = solve(eqs, vs, solution_dict=True)")
        print("slns = solve(eqs, vs)")

    def py_func(self):
        arg = ",".join(str(v) for v in self.vs)
        #lines = ["def f(%s):"%arg]
        lines = ["def f(x):"]
        lines.append("  %s = x" % (arg,))
        lines.append("  value = [")
        for eq in self.eqs:
            lines.append("    %s,"%(eq,))
        lines.append("  ]")
        lines.append("  return value")
        code = '\n'.join(lines)
        #print(code)
        ns = {}
        exec(code, ns, ns)
        return ns['f']
    
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


def main():

    # Commutative special Frobenius algebra
    
    system = System()
    
    F = system.array(dim, dim**2, "F") # mul
    G = system.array(dim, 1, "G") # unit
    
    D = system.array(dim**2, dim, "D") # comul
    E = system.array(1, dim, "E") # counit

    for item in system.items:
        print(item)
    
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
    
    #system.sage_dump()
    f = system.py_func()

    n = len(system.vs)
    x0 = numpy.random.rand(n)
    #print(x0)
    #v = system(x0)
    #print(v)

    from scipy.optimize import root
    print("solving...")
    solution = root(f, x0, method="lm", tol=1e-4) # 43 secs
    assert solution.success
    x = solution.x
    print("solution:")
    print(solution.x)
    print(f(x))

    system.show(x)


if __name__ == "__main__":

    main()



