#!/usr/bin/env python3

import sys
from math import log

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
            return

    else:
        idx = freevars[0]

        for val in range(1, maxval+1):

            x[idx] = val
            x1 = propagate(A, x)
            if x1 is None:
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




def main():

    N = argv.get("N", 9)

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
        
    #        print(row)
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
    
    m = len(rows)
    n = len(keys)
    #A = numpy.zeros((m, n), dtype=object)
    A = numpy.zeros((m, n), dtype=int)
    for i, row in enumerate(rows):
        for j, key in enumerate(keys):
            value = row.get(key, 0)
            A[i, j] = value
    
    #print("A.shape:", A.shape)
    #print(A)
    print(keys)

    if 0:
        desc = ''.join(['r' for key in keys])
        rows = [list(row) for row in A]
        header = ['W(%d,%d)'%key for key in keys]
        s = latex_nosep([header]+rows, desc)
    
        print(s)
        return

    maxval = argv.get("maxval", 2*N)

    x = zeros(n)
    x[0] = argv.get("x0", 1)
    x = propagate(A, x)

    count = 0
    solutions = []
    for v in solve(A, x, maxval=maxval):
        solutions.append(list(v))
        s = ' '.join(str(x).ljust(2) for x in v)
        print("solution:", s)
        count += 1

        vec = []
        for i in range(1, N):
            if i < N-i:
                idx = keys.index((N-i, i))
            else:
                idx = keys.index((i, N-i))
            vec.append(v[idx])
        #print("vec:", vec)
        pyplot.plot(vec)

    #s = latex_nosep(solutions)

    vec = []
    for i in range(1, N):
        vec.append(W(i, N-i))
    pyplot.plot(vec, 'x-')
    

    print("_maxdepth:", _maxdepth)
    print("solutions:", count)

    if argv.plot:
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



