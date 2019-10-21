#!/usr/bin/env python3

"""
higher genus weight enumerators.
Check _they satisfy various identities.
"""

from random import random, randint, shuffle
from functools import reduce
from operator import mul

import numpy

from bruhat.solve import parse, shortstr, shortstrx, span, int_scalar, find_kernel
from bruhat.solve import row_reduce, dot2
from bruhat.poly import Poly
from bruhat import element
from bruhat.gset import Group
from bruhat.argv import argv
from bruhat.util import choose, cross, all_perms
from bruhat.dev.geometry import all_codes


ring = element.Z

def named_poly(cs, names):
    items = {}
    for k, v in cs.items():
        tpl = tuple((names[i], exp) for (i, exp) in enumerate(k) if exp)
        items[tpl] = v
    return Poly(items, ring)


xpoly1 = lambda cs : named_poly(cs, "x_0 x_1".split())
ypoly1 = lambda cs : named_poly(cs, "y_0 y_1".split())
xpoly2 = lambda cs : named_poly(cs, "x_{00} x_{01} x_{10} x_{11}".split())
xpoly3 = lambda cs : named_poly(cs,
    "x_{000} x_{001} x_{010} x_{011} x_{100} x_{101} x_{110} x_{111}".split())
xpoly4 = lambda cs : named_poly(cs,
    """x_0000 x_0001 x_0010 x_0011 x_0100 x_0101 x_0110 x_0111
    x_1000 x_1001 x_1010 x_1011 x_1100 x_1101 x_1110 x_1111""".split())

x_0 = xpoly1({(1,0):1})
x_1 = xpoly1({(0,1):1})

y_0 = ypoly1({(1,0):1})
y_1 = ypoly1({(0,1):1})

x_00 = xpoly2({(1,0,0,0):1})
x_01 = xpoly2({(0,1,0,0):1})
x_10 = xpoly2({(0,0,1,0):1})
x_11 = xpoly2({(0,0,0,1):1})

x_000 = xpoly3({(1,0,0,0,0,0,0,0):1})
x_001 = xpoly3({(0,1,0,0,0,0,0,0):1})
x_010 = xpoly3({(0,0,1,0,0,0,0,0):1})
x_011 = xpoly3({(0,0,0,1,0,0,0,0):1})
x_100 = xpoly3({(0,0,0,0,1,0,0,0):1})
x_101 = xpoly3({(0,0,0,0,0,1,0,0):1})
x_110 = xpoly3({(0,0,0,0,0,0,1,0):1})
x_111 = xpoly3({(0,0,0,0,0,0,0,1):1})

x_0000 = xpoly4({(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0):1})
x_0001 = xpoly4({(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0):1})
x_0010 = xpoly4({(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0):1})
x_0011 = xpoly4({(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0):1})
x_0100 = xpoly4({(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0):1})
x_0101 = xpoly4({(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0):1})
x_0110 = xpoly4({(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0):1})
x_0111 = xpoly4({(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0):1})
x_1000 = xpoly4({(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0):1})
x_1001 = xpoly4({(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0):1})
x_1010 = xpoly4({(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0):1})
x_1011 = xpoly4({(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0):1})
x_1100 = xpoly4({(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0):1})
x_1101 = xpoly4({(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0):1})
x_1110 = xpoly4({(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0):1})
x_1111 = xpoly4({(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1):1})




def genus_enum1(G, verbose=False):
    m, n = G.shape
    cs = {}
    for v0 in span(G):
        exp = [0, 0]
        for i in range(n):
            exp[v0[i]] += 1
        exp = tuple(exp)
        cs[exp] = cs.get(exp, 0) + 1
    p = xpoly1(cs)

    return p


def genus_enum2(G, verbose=False):
    m, n = G.shape
    cs = {}
    items = list(span(G))
    for v0 in items:
        #print(".",end='',flush=True)
        for v1 in items:
            exp = [0, 0, 0, 0]
            #vv = numpy.array([2*v0, v1])
            for i in range(n):
                exp[2*v0[i] + v1[i]] += 1
            exp = tuple(exp)
            cs[exp] = cs.get(exp, 0) + 1
        #break
    #print()
    p = xpoly2(cs)

    return p


def genus_enum3(G, verbose=False):
    #print(shortstr(G))
    m, n = G.shape
    cs = {}
    items = list(span(G))
    for v0 in items:
        #print(".",end='',flush=True)
        for v1 in items:
          for v2 in items:
            exp = [0, 0, 0, 0, 0, 0, 0, 0]
            #vv = numpy.array([2*v0, v1])
            for i in range(n):
                exp[4*v0[i] + 2*v1[i] + v2[i]] += 1
            exp = tuple(exp)
            cs[exp] = cs.get(exp, 0) + 1
        #break
    #print()
    p = xpoly3(cs)

    return p


def genus_enum4(G, verbose=False):
    #print(shortstr(G))
    m, n = G.shape
    cs = {}
    exp = numpy.array([0]*16, dtype=int)
    items = list(span(G))
    for v0 in items:
        #print(".",end='',flush=True)
        for v1 in items:
          for v2 in items:
           for v3 in items:
            exp[:] = 0
            for i in range(n):
                exp[8*v0[i] + 4*v1[i] + 2*v2[i] + v3[i]] += 1
            key = tuple(exp)
            cs[key] = cs.get(key, 0) + 1
    p = xpoly4(cs)

    return p


def genus_enum5(G, verbose=False):
    #print(shortstr(G))
    m, n = G.shape
    cs = {}
    exp = numpy.array([0]*32, dtype=int)
    items = list(span(G))
    for v0 in items:
        #print(".",end='',flush=True)
        for v1 in items:
         for v2 in items:
          for v3 in items:
           for v4 in items:
            exp[:] = 0
            for i in range(n):
                exp[16*v0[i] + 8*v1[i] + 4*v2[i] + 2*v3[i] + v4[i]] += 1
            key = tuple(exp)
            cs[key] = cs.get(key, 0) + 1
    p = xpoly5(cs) # todo...

    return p


def test_enum(G, verbose=False):

    w1 = genus_enum1(G)
    w2 = genus_enum2(G)
    w3 = genus_enum3(G)

    if verbose:
        print(w1.flatstr())
        print(w2.flatstr())
        print(w3.flatstr())

    zero1 = xpoly1({})
    zero2 = xpoly2({})

    w2_10 = w2.substitute({ "x_00" : x_0, "x_10" : x_1, "x_01":zero1, "x_11":zero1 })
    assert w2_10 == w1

    w2_01 = w2.substitute({ "x_00" : x_0, "x_10" : zero1, "x_01":x_1, "x_11":zero1 })
    assert w2_01 == w1

    w2_11 = w2.substitute({ "x_00" : x_0, "x_10" : zero1, "x_01":zero1, "x_11":x_1 })
    assert w2_11 == w1

    lhs = w2.substitute({ "x_00" : zero1, "x_10" : x_0, "x_01":zero1, "x_11":x_1 })
    assert lhs == w1 or lhs == 0

    lhs = w2.substitute({ 
        "x_00" : x_0*y_0, 
        "x_10" : x_1*y_0, 
        "x_01" : x_0*y_1, 
        "x_11" : x_1*y_1,
    })
    w1_y = w1.substitute({ "x_0" : y_0, "x_1" : y_1, })
    rhs = w1 * w1_y
    assert lhs == rhs

    lhs = w3.substitute({ 
        "x_000" : x_00*y_0, 
        "x_010" : x_01*y_0, 
        "x_001" : x_00*y_1, 
        "x_011" : x_01*y_1,
        "x_100" : x_10*y_0, 
        "x_110" : x_11*y_0, 
        "x_101" : x_10*y_1, 
        "x_111" : x_11*y_1,
    })
    rhs = w2 * w1_y
    assert lhs == rhs

    lhs = w3.substitute({ 
        "x_000" : y_0*x_00, 
        "x_010" : y_0*x_10, 
        "x_001" : y_0*x_01, 
        "x_011" : y_0*x_11,
        "x_100" : y_1*x_00, 
        "x_110" : y_1*x_10, 
        "x_101" : y_1*x_01, 
        "x_111" : y_1*x_11,
    })
    rhs = w2 * w1_y
    assert lhs == rhs, "%s != %s"%(lhs, rhs)

    lhs = w3.substitute({ 
        "x_000" : y_0*x_00, 
        "x_010" : y_1*x_00, 
        "x_001" : y_0*x_01, 
        "x_011" : y_1*x_01,
        "x_100" : y_0*x_10, 
        "x_110" : y_1*x_10, 
        "x_101" : y_0*x_11, 
        "x_111" : y_1*x_11,
    })
    rhs = w2 * w1_y
    assert lhs == rhs, "%s != %s"%(lhs, rhs)

    if argv.verbose:
        print("G =")
        print(shortstr(G))
        print("w1 =", w1)
        print("w2 =", w2)
        print("w3 =", w3)
        print()


def test():

    G = parse("11 ..")
    test_enum(G)

    G = parse("11. .11")
    #test_enum(G, verbose=True)

    #print(p.flatstr())
    #return

    #G = parse("11.. .11. 1..1")
    #test_enum(G)

    for m in range(1, 5):
      for n in range(m, 6):
        for G in all_codes(m, n):
            test_enum(G)

    print("OK")


def hecke():

    n = 5 # cols

    codes = []
    lookup = {}
    for m in range(0, n): # rows
        for G in all_codes(m, n):
            codes.append(G)
            lookup[G.tostring()] = G
    N = len(codes)
    pairs = []
    for i in range(N):
      for j in range(i+1, N):
        GH = intersect(codes[i], codes[j])
        assert GH.tostring() in lookup, GH


def freeze(G):
    return (G.tostring(), G.shape)

def thaw(fG):
    s, shape = fG
    G = numpy.fromstring(s, dtype=int_scalar)
    G.shape = shape
    return G


def mk_latex(G):
    head = r"\left[\begin{array}{cc}"
    tail = r"\end{array}\right]"
    lines = ["&".join(str(i) for i in row)+"\\\\" for row in G]
    lines = [head] + lines + [tail]
    return "\n".join(lines)


def find_equ():

    m = argv.get("m", 3)
    n = argv.get("n", 6)

    print("m=%d n=%d" % (m, n))
    show = argv.show

    ps_1 = {}
    ps_2 = {}
    ps_3 = {}
    count = 0
    for G in all_codes(m, n):
        w1 = genus_enum1(G)
        w2 = genus_enum2(G)
        #w3 = genus_enum3(G)
        value = freeze(G)
        if show and w2 not in ps_2 and G[:, 0].sum():
            if argv.latex:
                print("$$" + mk_latex(G)+ "$$")
                print()
                print(w1.flatstr())
                print()
                print(w2.flatstr())
                print()
                print()
            else:
                print(G)
                print(w1.flatstr())
                print(w2.flatstr())
                print()
        if w1 not in ps_1:
            ps_1[w1] = value
        if w2 not in ps_2:
            ps_2[w2] = value
        #ps_3[w3] = value
        count += 1
    print(count, len(ps_1), len(ps_2), len(ps_3))
    v1 = set(ps_1.values())
    v2 = set(ps_2.values())
    v3 = set(ps_3.values())


    for fG in v2:
      if fG not in v1:
        G = thaw(fG)
        w1 = genus_enum1(G)
        fG1 = ps_1[w1]
        G1 = thaw(fG1)
        print(G)
        print(G1)
        print(w1.flatstr())
        print(genus_enum2(G).flatstr())
        print(genus_enum2(G1).flatstr())


def CX_12(w2):
    #w2 = w2.substitute({"x_00":x_00, "x_01":x_01, "x_10":x_10, "x_11":x_11})
    w2 = w2.substitute({"x_00":x_00, "x_01":x_01, "x_10":x_11, "x_11":x_10})
    return w2

def CX_21(w2):
    #w2 = w2.substitute({"x_00":x_00, "x_01":x_01, "x_10":x_10, "x_11":x_11})
    w2 = w2.substitute({"x_00":x_00, "x_01":x_11, "x_10":x_10, "x_11":x_01})
    return w2


def CCX(w3):
    w3 = w3.substitute({ 
        "x_000" : x_000, 
        "x_010" : x_010, 
        "x_001" : x_001, 
        "x_011" : x_011,
        "x_100" : x_100, 
        "x_110" : x_111, 
        "x_101" : x_101, 
        "x_111" : x_110,
    })
    return w3

def CCZ(w3):
    w3 = w3.substitute({ 
        "x_000" : x_000, 
        "x_010" : x_010, 
        "x_001" : x_001, 
        "x_011" : x_011,
        "x_100" : x_100, 
        "x_110" : x_110, 
        "x_101" : x_101, 
        "x_111" : -x_111,
    })
    return w3



def find_CX():
    m = argv.get("m", 3)
    n = argv.get("n", 6)

    assert CX_12(x_10) == x_11

    print("m=%d n=%d" % (m, n))
    show = argv.show

    ps_1 = {}
    ps_2 = {}
    ps_3 = {}
    count = 0
    found = set()
    for G in all_codes(m, n):
        #w1 = genus_enum1(G)
#        w2 = genus_enum2(G)
#        a = CX_12(w2)
#        b = CX_21(w2)
#        if w2 == a and w2 == b:
#            print("", end="", flush=True)
#        else:
#            print(".", end="", flush=True)
        w3 = genus_enum3(G)
        a = CCX(w3)
        b = CCZ(w3)
        if w3 == a and w3 == b:
            #print("3", end="", flush=True)
            #print()
            #print(G)
            pass
        elif w3 == a:
            #print("1", end="", flush=True)
            if w3 not in found:
                print(G)
        elif w3 == b:
            print("2", end="", flush=True)
        #else:
            #print(".", end="", flush=True)
        found.add(w3)


    print()


def get_phi_2(w2):
    ns = dict((k, 0) for k in globals().keys())
    items = []
    for ns1 in [
        {"x_00":x_00, "x_10":x_10},
        {"x_00":x_00, "x_01":x_01},
        {"x_00":x_00, "x_11":x_11},
        {"x_10":x_10, "x_01":x_01},
        {"x_10":x_10, "x_11":x_11},
        {"x_01":x_01, "x_11":x_11}]:
        ns2 = dict(ns)
        ns2.update(ns1)
        items.append(w2.substitute(ns2))
    return items


def get_phi_3(w3, w2=None):
    vs = dict(globals())
    ns = dict((k, 0) for k in globals().keys())
    items = []
    names = "x_000 x_001 x_010 x_011 x_100 x_101 x_110 x_111".split()
    count = 0
    for _ns1 in choose(names, 4):
        ns1 = dict((k, vs[k]) for k in _ns1)
        ns2 = dict(ns)
        ns2.update(ns1)
        w = w3.substitute(ns2)
        items.append(w)
        #print(' '.join(_ns1), w.otherstr())
#        if w2 is not None and w.otherstr() == w2.otherstr():
#            x = sum(eval("0b"+n[2:]) for n in _ns1) % 2
#            print(' '.join(_ns1), w.otherstr())
#            assert x==0, x
#            #print(w.otherstr(), w2.otherstr())
#            count += 1
#    print("count:", count)
    return items

def get_phi_4(w4, w3=None):
    vs = dict(globals())
    ns = dict((k, 0) for k in globals().keys())
    items = []
    names = """
        x_0000 x_0001 x_0010 x_0011 x_0100 x_0101 x_0110 x_0111
        x_1000 x_1001 x_1010 x_1011 x_1100 x_1101 x_1110 x_1111
    """.strip().split()
    count = 0
    for _ns1 in choose(names, 8):
        ns1 = dict((k, vs[k]) for k in _ns1)
        ns2 = dict(ns)
        ns2.update(ns1)
        w = w4.substitute(ns2)
        items.append(w)
#        #print(' '.join(_ns1), w.otherstr())
#        if w3 is not None and w.otherstr() == w3.otherstr():
#            x = sum(eval("0b"+n[2:]) for n in _ns1) % 2
#            print(' '.join(_ns1), "==", x)
#            count += 1
#    print("count:", count)
    return items



def test_phi_3():
    names = """
        x_000 x_001 x_010 x_011 x_100 x_101 x_110 x_111
    """.strip().split()
    count = 0
    for ns in choose(names, 4):
        x = sum(eval("0b"+n[2:]) for n in ns) % 2
        if x==0:  # <------ subscripts sum to zero
            print(' '.join(ns), "==", x)
            count += 1
    print("count:", count)



def test_phi_4():
    names = """
        x_0000 x_0001 x_0010 x_0011 x_0100 x_0101 x_0110 x_0111
        x_1000 x_1001 x_1010 x_1011 x_1100 x_1101 x_1110 x_1111
    """.strip().split()
    count = 0
    for ns in choose(names, 8):
        x = sum(eval("0b"+n[2:]) for n in ns) % 2
        if x==0:  # <------ subscripts sum to zero
            print(' '.join(ns), "==", x)
            count += 1
    print("count:", count)



def main():

    if argv.code == "RM14":
        s = """
        1111111111111111
        .1.1.1.1.1.1.1.1
        ..11..11..11..11
        ....1111....1111
        ........11111111
        """
    else:
        s = argv.get("G")
    if s is None:
        G = parse("111")
        #G = parse("1. .1")
        #G = parse("11. .11")
        #G = parse("1.. .11")
    else:
        G = parse(s)

    print(G)
    w1 = genus_enum1(G)
    w2 = genus_enum2(G)
    #print(w1)
    print("w1 =")
    print(w1.flatstr())
    #print(w2)
    print("w2 =")
    print(w2.flatstr())

    items = get_phi_2(w2)
    #for w in items:
    #    print("\t", w.flatstr())
    w = sum(items)
    print("sum(items) =")
    print(w.flatstr())
    print("diff:")
    print((w2 - w).flatstr())
    
    w3 = genus_enum3(G)
    print("w3 =")
    print(w3.flatstr())
    #for w in items:
    #    print("\t", w.flatstr(), w.flatstr()==w2.flatstr())
    items = get_phi_3(w3, w2)
    #w = sum(items)
    #print("diff:")
    #print((w3 - w).flatstr())
    

    if 0:
        w4 = genus_enum4(G)
        print("w4 =")
        print(w4.flatstr())
        #for w in items:
        #    print("\t", w.flatstr(), w.flatstr()==w2.flatstr())
        items = get_phi_4(w4, w3)
        #w = sum(items)
        #print("diff:")
        #print((w4 - w).flatstr())
    

if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:
        name = argv.next() or "main"

        eval(name)()


