#!/usr/bin/env python3
"""
Symbolically find a cocommutative comonoid 
on a 2d vector space.
"""

import numpy
from numpy import dot, array, empty, alltrue, exp
from numpy.linalg import norm

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

    def subs(self, values=None, x=None):
        if values is None:
            values = dict((str(v), xi) for (v, xi) in zip(self.vs, x))
        items = [A.copy() for A in self.items]
        for i, A in enumerate(items):
            for idx in numpy.ndindex(A.shape):
                A[idx] = A[idx].subs(values)
        items = [A.astype(float) for A in items]
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

    def py_func(self, verbose=False):
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
        if verbose:
            print(code)
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

    def root(self, trials=1, scale=1., method="lm", tol=1e-6):
        from scipy.optimize import root
        n = len(self.vs)
        f = self.py_func()
        for trial in range(trials):
            x0 = numpy.random.normal(size=n)*scale
            solution = root(f, x0, method=method, tol=tol, options={"maxiter":100000})
            if solution.success:
                break
            #print(solution)
        else:
            return None
        x = solution.x
        values = self.subs(x=x)
        return values

    def minimize(self, trials=10, scale=1, method=None, exclude=[], sigma=1000):
        from scipy.optimize import minimize
        n = len(self.vs)
        f = self.py_func()
        def f1(x):
            values = f(x)
            #result = sum([abs(v) for v in values])
            result = sum([v**2 for v in values])
            for v in exclude:
                r = sigma*(norm(v-x)**2)
                #print(x, exp(-r))
                result += exp(-r)
            return result
        for trial in range(trials):
            x0 = numpy.random.normal(size=n)*scale
            solution = minimize(f1, x0, tol=1e-7, method=method)
            if solution.success and solution.fun < 1e-5:
                break
        else:
            return None
        #print(solution)
        x = solution.x
        values = self.subs(x=x)
        return values

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

        
def test():
    dim = 2
    I = empty((dim, dim), dtype=object)
    I[:] = 0
    for i in range(dim):
        I[i, i] = 1
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

    dim = argv.get("dim", 3)
    dim2 = dim**2
    
    I = empty((dim, dim), dtype=float)
    I[:] = 0
    for i in range(dim):
        I[i, i] = 1
    #print(I)
    
    SWAP = empty((dim, dim, dim, dim), dtype=float)
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
    
    # -----------------------------------------------------------
    # Build some Frobenius algebra's ... find the copyable states

    system = System()
    
    F = system.array(dim, dim**2, "F") # mul
    G = system.array(dim, 1, "G") # unit
    
    D = system.array(dim**2, dim, "D") # comul
    E = system.array(1, dim, "E") # counit

    #for item in system.items:
    #    print(item)
    
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
    commutative = argv.get("commutative", True)
    if commutative:
        system.add(F, compose(SWAP, F))
    
    # counit 
    system.add(compose(D, IE), I)
    system.add(compose(D, EI), I)
    
    # _coassoc
    system.add(compose(D, DI), compose(D, ID))
    
    # cocommutative
    #system.add(D, compose(D, SWAP))
    
    # Frobenius
    system.add(compose(DI, IF), compose(F, D))
    system.add(compose(ID, FI), compose(F, D))
    
    # special
    special = argv.get("special", True)
    if special:
        system.add(compose(D, F), I)

    f = system.py_func()

    trials = argv.get("trials", 100)
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed:", _seed)
        numpy.random.seed(_seed)

    #for trial in range(trials):
    while 1:

        print(".", end="", flush=True)
        values = system.root()
        if values is None:
            continue
        F, G, D, E = values[:4]

        lhs = compose(SWAP, F)
        commutative = numpy.allclose(lhs, F, atol=1e-8)
        lhs = compose(D, SWAP)
        cocommutative = numpy.allclose(lhs, D, atol=1e-8)

        if commutative != cocommutative:
            lhs = compose(SWAP, F)
            print("\ncommutative:", commutative)
            print(lhs)
            print(F)
            print("\ncocommutative:", cocommutative)
            lhs = compose(D, SWAP)
            print(lhs)
            print(D)
        
        assert commutative == cocommutative

        ##break
        #if numpy.min(D) < 1e-4:
        #    continue

        break

    #print(F)
    #print(D)

    #A = dot(tensor(dot(E, F), I), tensor(I, dot(D, G)))
    #A = tensor(I, dot(D, G))
    A = dot(E, F)
    A.shape = (dim, dim)
    print("A =")
    print(A)
    vals = numpy.linalg.eig(A)[0]
    print("vals:", vals)

    def find_copyable(D):
        #print("\nminimize...")
        system = System()
        v = system.array(dim, 1)
        system.add(compose(v, D), tensor(v, v))
        #system.py_func(True)
        found = []
        trials = 0
        exclude = [numpy.zeros(dim)]
        exclude = []
        #while trials < 1000:
        while len(found) < dim and trials<1000:
            trials += 1
            values = system.root(scale=10)
            if values is None:
                values = system.minimize(scale=10, method="TNC", exclude=exclude, sigma=3000)
            if values is None:
                continue
            v = values[0]
            v = v.transpose()[0]
            #print(v)
            if norm(v) < 1e-2:
                continue
            for w in found:
                if norm(v-w) < 1e-2:
                    break
            else:
                found.append(v)
                print("found:", len(found))
                print(v)
                exclude.append(v)
        print()
        g = numpy.array(found).transpose()
        return g

    g = find_copyable(D)

    if g.shape == (dim,dim) and dim==3:
        # Check the Frobenius structure is invariant if
        # we swap basis elements around...
        print(D)
        print(g)
        
        gi = numpy.linalg.inv(g)
        P = numpy.array([[1,0,0],[0,0,1],[0,1,0]])
        Q = dot(g, dot(P, gi))
        Qi = numpy.linalg.inv(Q)
        QQ = tensor(Q, Q)
        QQi = tensor(Qi, Qi)
    
        D1 = dot(QQ, dot(D, Qi))
        print(D1)
    
        g1 = find_copyable(D1)
        print(g1)

        assert numpy.allclose(D, D1, rtol=1e-6)

        F1 = dot(Q, dot(F, QQi))
        assert numpy.allclose(F, F1, rtol=1e-6)


if __name__ == "__main__":

    if argv.test:
        test()
    else:
        main()



