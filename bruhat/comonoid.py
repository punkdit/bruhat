#!/usr/bin/env python3
"""
Symbolically find a cocommutative comonoid 
on a 2d vector space.
"""

from functools import reduce

import numpy
from numpy import dot, array, empty, alltrue, allclose, exp, zeros, kron
from numpy.linalg import norm

#def compose(first, second):
    #return numpy.dot(second, first)
def compose(*items):
    items = reversed(items)
    return reduce(numpy.dot, items)
#tensor = numpy.kron
tensor = lambda *items : reduce(numpy.kron, items)

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

    def subs(self, values=None, x=None, dtype=float):
        if values is None:
            values = dict((str(v), xi) for (v, xi) in zip(self.vs, x))
        items = [A.copy() for A in self.items]
        for i, A in enumerate(items):
            for idx in numpy.ndindex(A.shape):
                A[idx] = A[idx].subs(values)
        items = [A.astype(dtype) for A in items]
        return items

    def doshow(self, x):
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
    
    def py_jac(self, verbose=False):
        arg = ",".join(str(v) for v in self.vs)
        #lines = ["def f(%s):"%arg]
        lines = ["def f(x):"]
        #lines.append("  print('jac')")
        lines.append("  %s = x" % (arg,))
        lines.append("  value = [")
        vs = self.vs
        for eq in self.eqs:
          row = []
          for v in vs:
            d = eq.diff(v)
            row.append(str(d))
          lines.append("    [%s]," % (', '.join(row)))
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

    def root(self, trials=1, scale=1., method="lm", tol=1e-6, maxiter=1000000):
        from scipy.optimize import root
        n = len(self.vs)
        f = self.py_func()
        jac = self.py_jac()
        for trial in range(trials):
            x0 = numpy.random.normal(size=n)*scale
            solution = root(f, x0, jac=jac, method=method, 
                tol=tol, options={"maxiter":maxiter})
            #print(solution)
            if solution.success:
                break
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

        
def test_sympy():
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


def test():
    dim = 2
    I = empty((dim, dim), dtype=float)
    I[:] = 0
    for i in range(dim):
        I[i, i] = 1
    #M = system.array(2, 4, "M")
    M = array([
        [1, 0, 0, 0],
        [0, 0, 0, 1],
    ], dtype=float)

    system = System()
    G = system.array(2, 2, "G")
    Gi = system.array(2, 2, "Gi")
    GG = tensor(G, G)
    
    system.add(compose(G, Gi), I)
    #system.eqs.append(G[0,0]*G[1,1] - G[0,1]*G[1,0]-1) # SL(2)
    system.py_func(verbose=True)
    values = system.minimize()
    
    assert values is not None
    G, Gi = values
    print(compose(G, Gi))



def shortstr(F):
    s = str(F)
    s = s.replace(". ", "  ")
    s = s.replace(".]", " ]")
    s = s.replace(" 0 ", " . ")
    return s


def get_swap(dim):
    dim2 = dim**2
    SWAP = empty((dim, dim, dim, dim), dtype=float)
    SWAP[:] = 0
    for i in range(dim):
      for j in range(dim):
        SWAP[i, j, j, i] = 1
    SWAP.shape = dim2, dim2
    #print(SWAP)
    return SWAP


def get_mul(dim):
    F = zeros((dim, dim**2))
    for i in range(dim):
        F[i, i*dim + i] = 1
    return F

def get_comul(dim):
    F = zeros((dim**2, dim))
    for i in range(dim):
        F[i*dim + i, i] = 1
    return F

def got_comul(dim):
    F = zeros((dim**2, dim))
    for i in range(dim):
        F[i*dim + i, i] = 1
    return F

def get_unit(dim):
    F = zeros((dim, 1))
    F[:] = 1.
    return F

def get_counit(dim):
    F = zeros((1, dim))
    F[:] = 1.
    return F

def get_identity(dim):
    I = empty((dim, dim), dtype=float)
    I[:] = 0
    for i in range(dim):
        I[i, i] = 1
    return I


def test_root():
    dim = 4

    system = System()
    I = get_identity(dim)
    A = system.array(dim, dim)
    B = system.array(dim, dim)

    system.add(compose(A, A), I)
    system.add(compose(B, B), I)
    system.add(compose(A, B), compose(B, A))
    
    f = system.py_func(verbose=False)
    jac = system.py_jac(verbose=False)

    from scipy.optimize import root
    scale = 1.0
    n = len(system.vs)
    x0 = numpy.random.normal(size=n)*scale
    #jac = False
    solution = root(f, x0, jac=jac, method='lm', tol=1e-6, options={"maxiter":1000000})
    #print("solution.nfev:", solution.nfev)
    assert solution.success, solution
    
    x = solution.x
    values = system.subs(x=x)
    A, B = values
    #print(A)
    #print(B)
    assert numpy.allclose(compose(A, A), I)
    assert numpy.allclose(compose(B, B), I)
    assert numpy.allclose(compose(A, B), compose(B, A))


def test_symplectic():

    n = argv.get("n", 2)

    F = zeros((2*n, (2*n)**2), dtype=float)
    for i in range(n):
        F[i, i*2*n + i+n] = 1
        F[i+n, (i+n)*2*n + i ] = -1
#    print(F)

    for i in range(2*n):
      for j in range(2*n):
        a = zeros((2*n,))
        b = zeros((2*n,))
        a[i] = 1
        b[j] = 1
        u = kron(a, b)
        v = dot(F, u)
        #print(i, j, v)

    swap = get_swap(2*n)
    unit = get_unit(2*n)
    unit[:n] = -1
    mul = F
    #counit = get_counit(2*n)
    counit = unit.transpose()
    copy = get_comul(2*n)
    comul = mul.transpose()
    I = get_identity(2*n)

    print(shortstr(dot(mul, tensor(I, unit))))
    print(shortstr(dot(mul, tensor(unit, I))))

    lhs = dot(copy, mul)
    rhs = dot(tensor(mul, I), tensor(I, copy))
    print(allclose(lhs, rhs))
    print(shortstr(lhs))
    print(shortstr(rhs))

    return

    cup = dot(comul, unit)
    cap = dot(counit, F)
    print(shortstr(cap))
    print(shortstr(dot(tensor(I, cap), tensor(cup, I))))

    print(dot(F, comul))

#    lhs = dot(F, tensor(F, I))
#    rhs = dot(F, tensor(I, F))
#    print(shortstr(lhs))
#    print(shortstr(rhs))

#    # jacobi identity ? No...
#    IS = tensor(I, swap)
#    SI = tensor(swap, I)
#    FF = dot(F, tensor(I, F))
#    A = FF
#    B = dot(FF, dot(IS, SI))
#    C = dot(FF, dot(SI, IS))
#    print((A + B + C))


def main():

    dim = argv.get("dim", 3)
    
    I = empty((dim, dim), dtype=float)
    I[:] = 0
    for i in range(dim):
        I[i, i] = 1
    #print(I)
    
    SWAP = get_swap(dim)
    
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
    anticomm = argv.anticomm
    commutative = argv.get("commutative", not anticomm)
    if anticomm:
        #system.add(F, compose(-SWAP, F))
        FE = compose(F, E)
        system.add(FE, compose(-SWAP, FE))
    elif commutative:
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
        values = system.root(maxiter=10000000)
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


def main_hopf(dim=2):

    dim = argv.get("dim", dim)
    
    I = empty((dim, dim), dtype=float)
    I[:] = 0
    for i in range(dim):
        I[i, i] = 1
    #print(I)
    
    swap = get_swap(dim)
    
    # -----------------------------------------------------------
    # Build some Frobenius algebra's & Hopf algebra's

    system = System()
    
    # green spiders
    g_gg = D = system.array(dim**2, dim, "D") # comul
    g_   = E = system.array(1, dim, "E") # counit
    gg_g = F = system.array(dim, dim**2, "F") # mul
    _g   = G = system.array(dim, 1, "G") # unit

    # red spiders
    r_rr = J = system.array(dim**2, dim, "J") # comul
    r_   = K = system.array(1, dim, "K") # counit
    rr_r = L = system.array(dim, dim**2, "L") # mul
    _r   = M = system.array(dim, 1, "M") # unit

    r_cup = compose(_r, r_rr)
    g_cup = compose(_g, g_gg)
    g_cap = compose(gg_g, g_)
    r_cap = compose(rr_r, r_)

    def mk_frobenius(system, D, E, F, G):
        ID, DI = tensor(I, D), tensor(D, I)
        IE, EI = tensor(I, E), tensor(E, I)
        IF, FI = tensor(I, F), tensor(F, I)
        IG, GI = tensor(I, G), tensor(G, I)
    
        # unit
        system.add(compose(IG, F), I)
        system.add(compose(GI, F), I)
        
        # _assoc
        system.add(compose(FI, F), compose(IF, F))
        
        # counit 
        system.add(compose(D, IE), I)
        system.add(compose(D, EI), I)
        
        # _coassoc
        system.add(compose(D, DI), compose(D, ID))
        
        # Frobenius
        system.add(compose(DI, IF), compose(F, D))
        system.add(compose(ID, FI), compose(F, D))

    mk_frobenius(system, D, E, F, G)
    mk_frobenius(system, J, K, L, M)

    #print(tensor(I, swap, I).shape)
    #print(tensor(r_rr, r_rr).shape)
    lhs = compose(tensor(r_rr, r_rr), tensor(I, swap, I))
    lhs = compose(tensor(r_rr, r_rr), tensor(I, swap, I), tensor(gg_g, gg_g))
    rhs = compose(gg_g, r_rr)
    system.add(lhs, rhs)

    system.add(compose(gg_g, r_), tensor(r_, r_))
    system.add(compose(_g, r_rr), tensor(_g, _g))
    system.add(compose(_g, r_), numpy.array([[1.]]))

    system.add(compose(_r, r_rr, tensor(I, g_)), _g)
    system.add(compose(_r, r_rr, tensor(g_, I)), _g)

    system.add(compose(tensor(I, _r), gg_g, g_), r_)
    system.add(compose(tensor(_r, I), gg_g, g_), r_)

    inv = compose(tensor(I, r_cup), tensor(swap, I), tensor(I, g_cap))
    if argv.inv:
        print("inv")
        system.add(compose(g_gg, tensor(I, inv), rr_r), compose(g_, _r))
        system.add(compose(g_gg, tensor(inv, I), rr_r), compose(g_, _r))

    print("py_func")

    f = system.py_func()

    trials = argv.get("trials", 100)
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed:", _seed)
        numpy.random.seed(_seed)


    values = system.root(maxiter=10000000)
    if values is None:
        assert 0

    print()

    D, E, F, G, J, K, L, M = values

    g_gg = D
    g_   = E
    gg_g = F
    _g   = G

    r_rr = J
    r_   = K
    rr_r = L
    _r   = M

    print(r_rr)
    print(rr_r)
    print(g_gg)
    print(gg_g)

    print("commutative:")
    print("r_rr", allclose(compose(r_rr, swap), r_rr))
    print("g_gg", allclose(compose(g_gg, swap), g_gg))
    print("rr_r", allclose(compose(swap, rr_r), rr_r))
    print("gg_g", allclose(compose(swap, gg_g), gg_g))

    r_cup = compose(_r, r_rr)
    g_cup = compose(_g, g_gg)
    g_cap = compose(gg_g, g_)
    r_cap = compose(rr_r, r_)

    print("g_cap", allclose(compose(swap, g_cap), g_cap))
    print("r_cap", allclose(compose(swap, r_cap), r_cap))
    print("r_cup", allclose(compose(r_cup, swap), r_cup))
    print("g_cup", allclose(compose(g_cup, swap), g_cup))

    inv = compose(tensor(I, r_cup), tensor(swap, I), tensor(I, g_cap))

    return locals()

    linv = compose(tensor(I, r_cup), tensor(swap, I), tensor(I, g_cap))
    rinv = compose(tensor(r_cup, I), tensor(I, swap), tensor(g_cap, I))
    print(allclose(linv, rinv, atol=1e-4))
    inv = linv
    print(inv)

    lhs = compose(g_gg, tensor(I, inv), rr_r)
    rhs = compose(g_, _r)
    print(rhs)
    print(lhs)
    print(allclose(lhs, rhs, atol=1e-4))
    lhs = compose(g_gg, tensor(inv, I), rr_r)
    print(lhs)
    print(allclose(lhs, rhs, atol=1e-4))
    print()

    lhs = compose(r_rr, tensor(I, inv), gg_g)
    rhs = compose(r_, _g)
    print(rhs)
    print(lhs)
    print(allclose(lhs, rhs, atol=1e-4))
    lhs = compose(r_rr, tensor(inv, I), gg_g)
    print(lhs)
    print(allclose(lhs, rhs, atol=1e-4))
    print()


# -----------------------------------------------------------------------------
#
#

from bruhat.util import cross

scalar = numpy.int64 # ?


class Matrix(object):
    def __init__(self, p, a):
        a = numpy.array(a, dtype=scalar)
        a %= p
        self.a = a
        self.shape = a.shape

    @classmethod
    def zeros(self, p, shape):
        a = numpy.zeros(shape, dtype=scalar)
        return Matrix(p, a)

    def __add__(self, other):
        assert self.p == other.p
        a = self.a + other.a
        return Matrix(self.p, a)

    def __mul__(self, other):
        assert self.p == other.p
        a = dot(self.a, other.a)
        return Matrix(self.p, a)

    def __matmul__(self, other):
        assert self.p == other.p
        a = kron(self.a, other.a)
        return Matrix(self.p, a)



def ffield():
    "look for finite field solutions"

    p = argv.get("p", 3)
    dim = argv.get("dim", 3)
    dim2 = dim**2
    
    I = empty((dim, dim), dtype=int)
    I[:] = 0
    for i in range(dim):
        I[i, i] = 1
    #print(I)
    
    SWAP = empty((dim, dim, dim, dim), dtype=int)
    SWAP[:] = 0
    for i in range(dim):
      for j in range(dim):
        SWAP[i, j, j, i] = 1
    SWAP.shape = dim2, dim2
    #print(SWAP)
    
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

    # https://theory.stanford.edu/~nikolaj/programmingz3.html
    import z3

    def all_smt(s, initial_terms):
        def block_term(s, m, t):
            s.add(t != m.eval(t, model_completion=True))
        def fix_term(s, m, t):
            s.add(t == m.eval(t, model_completion=True))
        def all_smt_rec(terms):
            if z3.sat == s.check():
               m = s.model()
               yield m
               for i in range(len(terms)):
                   s.push()
                   block_term(s, m, terms[i])
                   for j in range(i):
                       fix_term(s, m, terms[j])
                   yield from all_smt_rec(terms[i:])
                   s.pop()
        yield from all_smt_rec(list(initial_terms))

    def block_model(s):
        m = s.model()
        s.add(z3.Or([f() != m[f] for f in m.decls() if f.arity() == 0]))

    ns = {}
    clauses = []
    for v in system.vs:
        v = str(v)
        iv = z3.Int(v)
        ns[v] = iv
        clauses.append(iv < p)
        clauses.append(0 <= iv)
    print(ns)

    for eq in system.eqs:
        eq = eval(str(eq), ns) % p == 0
        clauses.append(eq)

    solver = z3.Solver()

    for c in clauses:
        solver.add(c)

#    print("solve...")
#    result = solver.check()
#    assert result == z3.sat, result
#    print("done.")

    #for model in all_smt(solver, []):

    while solver.check() == z3.sat:
        model = solver.model()

        values = {}
        for d in model:
            #print(d, model[d])
            values[str(d)] = int(str(model[d]))
        print(values)
    
        items = []
        for A in system.items:
            B = numpy.zeros(A.shape, int)
            for idx in numpy.ndindex(A.shape):
                B[idx] = values[str(A[idx])]
            items.append(B)
    
        F, G, D, E = items
        print(F)
    
        for vals in cross([tuple(range(dim))]*dim):
            v = numpy.array(vals)
            v.shape = (dim,1)
            if v.sum() == 0:
                continue
            lhs = dot(D, v)
            rhs = tensor(v, v)
            if numpy.allclose(lhs, rhs):
                print(v.transpose())
    
        block_model(solver)
        print()

    print("done")


if __name__ == "__main__":

    name = argv.next() or "main"
    fn = eval(name)

    print("%s()"%name)
    fn()

    print("OK")
    print()





