#!/usr/bin/env python3

"""
https://golem.ph.utexas.edu/category/2019/03/how_much_work_can_it_be_to_add.html#c055688
"""

import sys
from math import log
from random import shuffle, choice

import numpy
from matplotlib import pyplot

#from bruhat.gelim import row_reduce, shortstr, kernel
from qupy.dev import linalg
from bruhat.argv import argv

EPSILON = 1e-8

def W(*items):

    sitems = sum(items)
    r = 0.
    for n in items:
        r += n * log(n)
    return -1*(r - sitems*log(sitems))


if 0:
    def getname(*items):
        items = list(items)
        items.sort(reverse=True)
        return "W(%s)"%(','.join(str(i) for i in items))
    
    def namerank(key):
        assert key[0] == "W"
        i, j = eval(key[1:])
        return (i+j, -i, -j)

else:
    def getname(*items):
        items = list(items)
        items.sort(reverse=True)
        return tuple(items)
    
    def namerank(key):
        i, j = key
        return (i+j, -i, -j)


def fstr(x):
    s = "%.4f"%x 
    s = s.rjust(8)
    return s


def astr(A):
    m, n = A.shape
    lines = []
    for i in range(m):
        row = ' '.join(fstr(x) for x in A[i, :])
        lines.append(row)
    return '\n'.join(lines)


def zeros(n):
    return numpy.zeros((n,), dtype=int)


def support(A):
    # non-zero indexes
    return numpy.where(A)[0]


def propagate(A, x):

    x = x.copy()
    m, n = A.shape

    changed = True
    while changed:
        changed = False

        #print("propagate", x)

        for i in range(m):
            row = A[i]
            freevar = None
            value = 0
            for j in support(row):
                if x[j] != 0:
                    value += row[j] * x[j]
                elif freevar is None:
                    freevar = j
                else:
                    # too many free variables
                    break
            else:
                if freevar is None and numpy.dot(row, x) != 0:
                    return None # no freevar and not a solution. Fail.

                if freevar is not None:
                    #print("row:", row, "freevar =", freevar, "value =", value)
                    value = -row[freevar] * value
                    if value > 0:
                        x[freevar] = value
                        changed = True
                    else:
                        return None # no positive solution

    return x


_maxdepth = 0
def solve(A, x=None, depth=0, maxval=5):

    global _maxdepth
    _maxdepth = max(depth, _maxdepth)

    m, n = A.shape

    x = x.copy()

    freevars = numpy.where(x==0)[0]
    #print('  '*depth, "solve")
    #print('  '*depth, "freevars:", freevars)

    if len(freevars) == 0:
        y = numpy.dot(A, x)
        if numpy.allclose(y, 0):
            yield x

        else:
            # Fail, backtrack... ?
            assert 0, "how did we get here... do we care?"
            return

    else:

        vals = list(range(1, maxval+1))
        if argv.random:
            shuffle(vals)

        for val in vals:

            if argv.random and 0:
                idx = choice(freevars)
            else:
                idx = freevars[0]

            x[idx] = val
            x1 = propagate(A, x)
            if x1 is None:
                #print("X", end="", flush=True)
                continue

            for x2 in solve(A, x1, depth+1, maxval):
                yield x2


def latex_nosep(rows, desc):
    n = len(rows[0])
    lines = []
    #lines.append("$$")
    lines.append(r"\begin{array}{%s}"%(desc))
    for i, row in enumerate(rows):
        row = list(row)
        if type(row)==list:
            line = " & ".join(str(fld) for fld in row) + r" \\"
        else:
            line = str(row)
        lines.append(line)
    lines.append(r"\end{array}")
    #lines.append("$$")
    s = "\n".join(lines)
    return s



def do_pulp(keys, rows):
    # Try to refute entropy inequalities using linear programming.
    import pulp
    system = pulp.LpProblem("system", pulp.LpMinimize)

    names = {}
    vs = {}
    for (i, j) in keys:
        name = "W_%d_%d"%(i, j)
        names[i, j] = name
        var = pulp.LpVariable(name)
        vs[i, j] = var
        system += var >= 1.0

    for i, row in enumerate(rows):
        expr = 0
        for j, key in enumerate(keys):
            value = row.get(key, 0)
            expr += value * vs[key]
        #print(expr)
        system += expr == 0.0

    #system += vs[3, 1] >= vs[2, 2] # Infeasible at N=9
    #system += vs[5, 1] >= vs[4, 2] # Infeasible at N=25
    system += vs[4, 2] >= vs[3, 3]  # Infeasible at N=64

    status = system.solve()
    s = pulp.LpStatus[status]
    print("pulp:", s)
    if s == "Infeasible":
        return

    if 0:
        for key in keys:
            v = vs[key]
            print(v, "=", pulp.value(v))
        


def main():

    N = argv.get("N", 9)
    modulo = argv.get("modulo")

    # build a linear system of equations
    rows = []
    keys = set()
    
    #keys.add((1, 1))
    #rows.append({(1, 1) : 1})
    
    for n in range(4, N+1):
    
        # a + b + c == n
        triples = set()
        for a in range(1, n):
         for b in range(1, n):
          for c in range(1, n):
            if a+b+c != n:
                continue
            items = [a, b, c]
            items.sort(reverse=True)
            triples.add(tuple(items))
        triples = list(triples)
        triples.sort()
    #    print("n =", n)
    #    print(triples)
    
        for a, b, c in triples:
    
            row = {}
            lhs = [getname(a, b), getname(a+b, c)]
            lhs.sort()
            rhs = [getname(a, b+c), getname(b, c)]
            rhs.sort()
    
            if lhs == rhs:
                continue
    
            #print(lhs, "=", rhs)
    
            for key in lhs:
                row[key] = 1
                keys.add(key)
    
            for key in rhs:
                row[key] = -1
                keys.add(key)
        
            #print(row)
            rows.append(row)
    
    for key in keys:
        i, j = key
        m = 2
        while 1:
            keym = (m*i, m*j)
            if keym not in keys:
                break
            row = {key : m, keym : -1}
            #print(key, keym, row)
            rows.append(row)
            m += 1
    
    keys = list(keys)
    keys.sort(key = namerank)
    #print("keys:", keys)

    #for k in keys:
    #    print("W(%d,%d)"%k)
    
    m = len(rows)
    n = len(keys)
    #A = numpy.zeros((m, n), dtype=object)
    A = numpy.zeros((m, n), dtype=int)
    for i, row in enumerate(rows):
        for j, key in enumerate(keys):
            value = row.get(key, 0)
            A[i, j] = value

    if argv.pulp:
        do_pulp(keys, rows)
        return
    
    #print("A.shape:", A.shape)
    #print(A)
    #print(keys)

    if 0:
        desc = ''.join(['r' for key in keys])
        rows = [list(row) for row in A]
        header = ['W(%d,%d)'%key for key in keys]
        s = latex_nosep([header]+rows, desc)
    
        print(s)
        return

    if argv.kernel:
        A = A.astype(float)
        B = linalg.row_reduce(A, truncate=True)
        m, n = B.shape
        print("B.shape:", B.shape)
        print("kernel:", n-m)
        return

    maxval = argv.get("maxval", 2*N)
    print("maxval:", maxval)

    x = zeros(n)
    x[0] = x0 = argv.get("x0", 1)
    x = propagate(A, x)

    maxcount = argv.maxcount
    N0 = argv.get("N0", N)
    ys = range(1, N0)

    solutions = []
    if argv.random:

        assert argv.maxcount is not None
        while len(solutions) < argv.maxcount:
            try:
                v = solve(A, x, maxval=maxval).__next__()
            except StopIteration:
                continue
            solutions.append(list(v))

            if argv.verbose:
                s = ' '.join(str(x).ljust(2) for x in v)
                print("solution:", s)
    
    else:
        for v in solve(A, x, maxval=maxval):
            solutions.append(list(v))
    
            if argv.verbose:
                s = ' '.join(str(x).ljust(2) for x in v)
                print("solution:", s)
    
            if maxcount is not None and len(solutions)>maxcount:
                break


    if 0:
        for v in solutions:
            for idx, k in enumerate(keys):
                value = v[idx]
                if argv.modulo is not None:
                    value %= argv.modulo
                i, j = k
                print("W(%d,%d)=%s"%(i, j, value), end=" ")
            print()

    #s = latex_nosep(solutions)

    print("maxdepth:", _maxdepth)
    print("solutions:", len(solutions))

    #lookup

    if argv.plot and not argv.random:
        vecs = set()
        for v in solutions:
            vec = []
            for i in range(1, N0):
                if i < N0-i:
                    idx = keys.index((N0-i, i))
                else:
                    idx = keys.index((i, N0-i))
                vec.append(v[idx])
            #print("vec:", vec)
            vec = tuple(vec)
            vecs.add(vec)
        print("vecs:", len(vecs))

        v0 = sum(vec[0] for vec in vecs) / len(vecs)
        for vec in vecs:
            pyplot.plot(ys, vec)

        vec = []
        r = v0 / W(1, N0-1)
        for i in range(1, N0):
            vec.append(r*W(i, N0-i))
        pyplot.plot(ys, vec, 'x-')
        pyplot.show()

    elif argv.plot:

        idxs = []
        for i in range(1, N0):
            if i < N0-i:
                idx = keys.index((N0-i, i))
            else:
                idx = keys.index((i, N0-i))
            idxs.append(idx)

        if argv.average:
            vec = []
            for idx in idxs:
                val = sum(v[idx] for v in solutions) / len(solutions)
                vec.append(val)
            pyplot.plot(ys, vec)
        else:
            for v in solutions:
                vec = [v[idx] for idx in idxs]
                pyplot.plot(ys, vec)

        #r = vec[0] / W(1, N0-1)
        r = x0 / W(1, 1)
        vec = []
        for i in range(1, N0):
            vec.append(r*W(i, N0-i))
        pyplot.plot(ys, vec, 'x-')

        pyplot.show()

    return

    A = A.astype(float)
    B = linalg.row_reduce(A, truncate=True)
    print("B.shape:", B.shape)
    print(B)

    v = numpy.array([W(i, j) for (i, j) in keys])
    #print(numpy.dot(A, v))
    assert numpy.abs(numpy.dot(B, v)).sum() < EPSILON

    return

    K = linalg.kernel(B)
    BK = numpy.dot(B, K)
    assert numpy.abs(BK).sum() < EPSILON
    K = K.transpose()

    K = linalg.row_reduce(K)

    print(keys)
    print("kernel:")
    print(repr(K))
    print(astr(K))

    if N==6:
        v = [1, 1, 2, 2, 1, 2, 2, 2, 3]
        v = [1, 2, 1, 2, 1, 1, 4, 4, 3]
        v = [1, 2, 1, 2, 2, 2, 3, 4, 3]
        print(numpy.dot(B, v))



if __name__ == "__main__":

    main()



