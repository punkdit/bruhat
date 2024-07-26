#!/usr/bin/env python

"""

"""

import string

import numpy

from sage.all_cmdline import GF, NumberField
from sage import all_cmdline 

from bruhat.util import all_primes
from bruhat.argv import argv
from bruhat.smap import SMap



def slow_count_points(p, l):
    q = p**l

    F = GF(q)

    e = F.gen()
    if l==1:
        # e is additive gen
        items = [i*e for i in range(q)]
    else:
        # e is multiplicative gen
        items = [e**i for i in range(1, q)] + [0*e]
    #print(p, l, len(set(items)), len(items))
    assert len(set(items)) == len(items) == q, [e, e**2]

    f = lambda x,y : y**2 == (x-1)*x*(x+1)
    sols = [f(x,y) for x in items for y in items]
    n = sols.count(True)

    lhs = lambda y : y**2
    rhs = lambda x : (x-1)*x*(x+1)

    lsend = {a : 0 for a in items}
    for y in items:
        lsend[lhs(y)] += 1

    rsend = {a : 0 for a in items}
    for x in items:
        rsend[rhs(x)] += 1

    count = 0
    for a in items:
        count += lsend[a]*rsend[a]
    assert count == n

    return n+1 # add projective point


def count_points(p, l, lhs, rhs):
    q = p**l

    F = GF(q)

    e = F.gen()
    if l==1:
        # e is additive gen
        items = [i*e for i in range(q)]
    else:
        # e is multiplicative gen
        items = [e**i for i in range(1, q)] + [0*e]
    #print(p, l, len(set(items)), len(items))
    #assert len(set(items)) == len(items) == q, [e, e**2]
    assert len(items) == q, [e, e**2]

    #f = lambda x,y : y**2 == (x-1)*x*(x+1)
    #sols = [f(x,y) for x in items for y in items]
    #n = sols.count(True)
    #print('.', end='', flush=True)

    lsend = {a : 0 for a in items}
    rsend = {a : 0 for a in items}
    for a in items:
        lsend[lhs(a)] += 1
        rsend[rhs(a)] += 1
        # does not help much:
        #lsend[a**2] += 1
        #rsend[(a-1)*a*(a+1)] += 1
    #print('. ', end='', flush=True)

    count = 0
    for a in items:
        count += lsend[a]*rsend[a]
    #assert count == n

    return count+1 # add projective point


def test():

    # good reduction at p=2
    lhs = lambda y : y**2 + y
    rhs = lambda x : x**3 + x
    p = 2

    # additive reduction at p=2
    lhs = lambda y : y**2
    rhs = lambda x : x**3
    p = 2

    # split multiplicative reduction at p=5
    lhs = lambda y : y**2
    rhs = lambda x : x**3 - x**2 + 5
    p = 5
    
    # nonsplit multiplicative reduction at p=3
    lhs = lambda y : y**2
    rhs = lambda x : x**3 - x**2
    p = 3
    
    ls = list(range(1, 9))
    for l in ls:
        print(p, l, p**l, end=" ", flush=True)
        n = count_points(p, l, lhs, rhs)
        print(n, n-(p**l+1))

def test_bolza():

    lhs = lambda y : y**2
    rhs = lambda x : x**5 - x
    p = argv.get("p", 5)
    
    l = argv.get("l", 7)
    ls = list(range(1, l))
    for l in ls:
        print(p, l, p**l, end=" ", flush=True)
        n = count_points(p, l, lhs, rhs)
        print(n, n-(p**l+1))


def find_coeffs(F, degree=2):
    from bruhat.gelim import array, solve
    f = [int(str(F[i])) for i in range(2*degree)]
    A = [[f[i+j] for i in range(degree)] for j in range(degree)]
    #A = array([[f[0], f[1]], [f[1], f[2]]])
    A = array(A)
    #rhs = [f[2], f[3]]
    rhs = [f[i] for i in range(degree, 2*degree)]
    v = solve(A, rhs)
    assert v is not None
    coeffs = [int(str(x)) for x in v]
    return coeffs


def test_jacobian_bolza():

    lhs = lambda y : y**2
    rhs = lambda x : x**3 +x**2 -3*x + 1
    p = argv.get("p", 3)
    
    ls = list(range(1, 7))
    diffs = []
    for l in ls:
        print(p, l, p**l, end=" ", flush=True)
        n = count_points(p, l, lhs, rhs)
        diff = n-(p**l+1)
        print(n, diff)
        diffs.append(diff)
    a, b = find_coeffs(diffs)
    print(a, b)


def main():

    # Ref.
    # Advanced Topics in the Arithmetic of Elliptic Curves
    # Joseph H. Silverman
    # Appendix A
    d = argv.get("d", 4)

    lhs = lambda y : y**2
    if d == 3:
        # eisenstein
        rhs = lambda x : x**3+1
    elif d == 4:
        # gauss
        rhs = lambda x : (x-1)*x*(x+1)
        #rhs = lambda x : x**3 - 11*x + 14
        #rhs = lambda x : x**3 + x
    elif d == 7:
        # y**2 + x*y == x**3 - x**2 - 2*x - 1 # conductor == 1
        rhs = lambda x : x**3 - 595*x + 5586 # conductor == 2 (??)
    elif d == 8:
        rhs = lambda x : x**3 + 4*x**2 + 2*x
    elif d == 11:
        lhs = lambda y : y**2 + y
        rhs = lambda x : x**3 - x**2 -7*x + 10
    elif d == 19:
        lhs = lambda y : y**2 + y
        rhs = lambda x : x**3 - 38*x + 90
    elif d == 43:
        lhs = lambda y : y**2 + y
        rhs = lambda x : x**3 - 860*x + 9707
    elif d == 67:
        lhs = lambda y : y**2 + y
        rhs = lambda x : x**3 - 7370*x + 243528
    elif d == 163:
        lhs = lambda y : y**2 + y
        rhs = lambda x : x**3 - 2174420*x + 1234136692
    else:
        print("d not found")
        return

    p = argv.get("p")
    l = argv.get("l", 1)
    if argv.all_primes:
        ps = list(all_primes(60))
        ls = [1,2,3,4]
    elif p is None:
        ps = [2, 3, 5, 7, 11]
        ls = list(range(1,5))
    else:
        ps = [p]
        if p < 5:
            ls = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        elif p < 10:
            ls = [1, 2, 3, 4, 5, 6]
        elif p < 50:
            ls = [1, 2, 3, 4]
        else:
            ls = [1, 2, 3]

    print(ps)

    for p in ps:
        if argv.mod is not None and p%d != argv.mod:
            continue
        print("p=%d, p mod %d=%d"%(p, d, p%d) )
        sols = []
        for l in ls:
            #print(p, l, p**l, end=" ", flush=True)
            n = count_points(p, l, lhs, rhs)
            diff=n-(p**l+1)
            sols.append(diff)
            #print(n, diff)
        if len(sols)>3:
            a, b = find_coeffs(sols)
            a = -a
            b = -b
            c = 1
            disc = b**2 - 4*a*c
            print("1 + %dx + %dx^2, disc=%s"%(b, a, disc), disc%d==0)


        print()


def test_number():
    n = argv.get("n")

    K0 = NumberField(x**2 + 1, "x")
    print(K0)

    K1 = K0.extension(x**2+x+1, "w")
    print(K1)

    print(K1.absolute_field("a").galois_group())


def count_frobenius():
    def order(A, n):
        i = 1
        B = A%n
        while not numpy.alltrue(B==I):
            i += 1
            B = (B@A)%n
        return i
    I = numpy.array([[1,0],[0,1]])

    A = numpy.array([[2,-1],[1,2]])
    for j in range(1,13):
        print("j=%d"%j, "3^j=%d"%(3**j), "order =",order(A, 3**j))


if __name__ == "__main__":

    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


