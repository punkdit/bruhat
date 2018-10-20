#!/usr/bin/env python3

"""
Build Reed-Muller codes.
"""


import numpy

from solve import array2, zeros2, dot2, shortstr, rank, find_kernel, span
from solve import linear_independent
from argv import argv
from util import choose


class Code(object):
    """
        binary linear code, as defined by a generator matrix.
    """
    def __init__(self, G, H=None, d=None, desc="", check=True):
        assert len(G.shape)==2
        self.G = G.copy()
        self.k, self.n = G.shape
        self.d = d
        self.desc = desc

        if H is None:
            H = list(find_kernel(G))
            H = array2(H)
        if H.shape == (0,):
            H.shape = (0, self.n)
        self.H = H.copy()

        if check:
            self.check()

    def check(self):
        G, H = self.G, self.H
        assert rank(G) == len(G)
        A = dot2(H, G.transpose())
        assert A.sum()==0

    def __str__(self):
        desc = ', "%s"'%self.desc if self.desc else ""
        return "Code([[%s, %s, %s]]%s)" % (self.n, self.k, self.d, desc)

    def dump(self):
        G, H = self.G, self.H
        print("G =")
        print(shortstr(G))
    
        print("H =")
        print(shortstr(H))

    def is_selfdual(self):
        #G = self.G
        #x = dot2(G, G.transpose())
        #return x.sum()==0
        return self.eq(self.get_dual())

    def get_dual(self):
        return Code(self.H, self.G)

    def eq(self, other):
        "Two codes are equal if their generating matrices have the same span."
        G1, G2 = self.G, other.G
        if len(G1) != len(G2):
            return False
        A = dot2(self.H, other.G.transpose())
        B = dot2(other.H, self.G.transpose())
        assert (A.sum()==0) == (B.sum()==0)
        return A.sum() == 0

    def get_distance(self):
        G = self.G
        d = None
        for v in span(G):
            w = v.sum()
            if w==0:
                continue
            if d is None or w<d:
                d = w
        if self.d is None:
            self.d = d
        return d

    def puncture(self, i):
        assert 0<=i<self.n
        G = self.G
        A = G[:, :i]
        B = G[:, i+1:]
        G = numpy.concatenate((A, B), axis=1)
        G = linear_independent(G, check=True)
        return Code(G)

    def is_morthogonal(self, m):
        G = self.G
        k = self.k
        assert m>=2
        if m>2 and not self.is_morthogonal(m-1):
            return False
        items = list(range(k))
        for idxs in choose(items, m):
            v = G[idxs[0]]
            for idx in idxs[1:]:
                v = v * G[idx]
            if v.sum()%2 != 0:
                return False
        return True

    def is_triorthogonal(self):
        return self.is_morthogonal(3)

    def weight_enum(self, A, B):
        G = self.G
        the_op = None
        for v in span(G):
            op = A if v[0]==0 else B
            for vi in v[1:]:
                op = op * (A if vi==0 else B)
            the_op = op if the_op is None else the_op + op
        return the_op


def reed_muller(r, m, puncture=False):
    "Build Reed-Muller code"

    assert 0<=r<=m, "r=%s, m=%d"%(r, m)

    n = 2**m # length

    one = array2([1]*n)
    basis = [one]

    vs = [[] for i in range(m)]
    for i in range(2**m):
        for j in range(m):
            vs[j].append(i%2)
            i >>= 1
        assert i==0

    vs = [array2(v) for v in vs]

    for k in range(r):
        for items in choose(vs, k+1):
            v = one
            #print(items)
            for u in items:
                v = v*u
            basis.append(v)
        
    G = numpy.array(basis)

    code = Code(G, d=2**(m-r), desc="reed_muller(%d, %d)"%(r, m))

    if puncture:
        code = code.puncture(0)

    return code


class Tensor(object):

    """ Some kind of graded ring element... I*I*I + X*X*X etc.
    """

    zero = 0
    one = 1
    def __init__(self, items, grade=None):
        # map key -> coeff, key is ("A", "B") etc.
        assert items or (grade is not None)
        keys = list(items.keys())
        keys.sort()
        self.items = {}
        nz = []
        for key in keys:
            assert grade is None or grade==len(key)
            grade = len(key)
            v = int(items[key])
            if v != 0:
                self.items[key] = v # uniquify
                nz.append(key)
        self.keys = nz
        self.grade = grade

    def __add__(self, other):
        assert self.grade == other.grade
        items = dict(self.items)
        for (k, v) in other.items.items():
            items[k] = items.get(k, self.zero) + v
        return Tensor(items, self.grade)

    def __sub__(self, other):
        assert self.grade == other.grade
        items = dict(self.items)
        for (k, v) in other.items.items():
            items[k] = items.get(k, self.zero) - v
        return Tensor(items, self.grade)

    def __mul__(self, other):
        items = {}
        for (k1, v1) in self.items.items():
          for (k2, v2) in other.items.items():
            k = k1+k2
            assert k not in items
            items[k] = v1*v2
        return Tensor(items, self.grade+other.grade)
    tensor = __mul__

    def __rmul__(self, r):
        items = {}
        for (k, v) in self.items.items():
            items[k] = r*v
        return Tensor(items, self.grade)

    def subs(self, rename):
        the_op = Tensor({}, self.grade) # zero
        one = self.one
        for (k, v) in self.items.items():
            final = None
            for ki in k:
                op = rename.get(ki, Tensor({ki : one}))
                if final is None:
                    final = op
                else:
                    final = final * op # tensor
            the_op = the_op + v*final
        return the_op

    def __str__(self):
        ss = []
        for k in self.keys:
            v = self.items[k]
            s = ''.join(str(ki) for ki in k)
            if v == 1:
                pass
            elif v == -1:
                s = "-"+s
            else:
                s = str(v)+"*"+s
            ss.append(s)
        ss = '+'.join(ss) or "0"
        ss = ss.replace("+-", "-")
        return ss

    def __repr__(self):
        return "Tensor(%s)"%(self.items)

    def __eq__(self, other):
        return self.items == other.items

    def __ne__(self, other):
        return self.items != other.items

    def __hash__(self):
        return hash((str(self), self.grade))



def test():

    for m in range(2, 7):
      for r in range(0, m+1):
        code = reed_muller(r, m)
        assert code.n == 2**m
        k = 1
        for i in range(1, r+1):
            k += len(list(choose(list(range(m)), i)))
        assert code.k == k
        if code.k < 12:
            assert code.get_distance() == 2**(m-r)

        if 0<=r<=m-1:
            dual = code.get_dual()
            code1 = reed_muller(m-r-1, m)
            #print(m-r-1 == r, dual.eq(code), code)
            assert code1.eq(dual)

    print("OK")


def gen():
    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)

    code = reed_muller(r, m)

    print(code)
    #print("d =", code.get_distance())
    code.dump()

    code = code.puncture(3)

    print(code)
    code.dump()
    print("d =", code.get_distance())


    for m in range(2, 10):
      for r in range(0, m+1):
        code = reed_muller(r, m)
        print(code, end=" ")
        if code.is_selfdual():
            print("is_selfdual", end=" ")
        if code.is_morthogonal(3):
            print("is_triorthogonal", end=" ")
        p = code.puncture(0)
        if p.is_morthogonal(3):
            print("puncture.is_triorthogonal", end=" ")
        print()


def main():
    I = Tensor({"I" : 1})
    X = Tensor({"X" : 1})
    Y = Tensor({"Y" : 1})
    Z = Tensor({"Z" : 1})

    II = I*I
    XI = X*I
    IX = I*X
    XX = X*X
    assert II+II == 2*II

    assert X*(XI + IX) == X*X*I + X*I*X

    assert ((I-Y)*I + I*(I-Y)) == 2*I*I - I*Y - Y*I
    assert (XI + IX).subs({"X": I-Y}) == ((I-Y)*I + I*(I-Y))

    A = Tensor({"A":1})
    B = Tensor({"B":1})
    p = A*A*A + B*B*A + B*A*B + A*B*B
    q = A*A*A + B*B*B
    p1 = p.subs({"A": A+B, "B": A-B})
    assert p1 == 4*A*A*A + 4*B*B*B

    verbose = argv.verbose

    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)

    code = reed_muller(r, m)
    code.dump()

    p = code.weight_enum(A, B)
    if verbose:
        print("code:")
        print(p)

    dual = code.get_dual()
    q = dual.weight_enum(A, B)
    if verbose:
        print("dual:")
        print(q)
    print("p==q:", p==q)
    print("code.is_selfdual:", code.is_selfdual())

    #r = p.subs({"A": A+B, "B": A-B})
    r = code.weight_enum(A+B, A-B)
    if verbose:
        print("P(A+B, A-B)")
        print(r)
    coeff = 2**len(code.G)
    print("MacWilliams:", r == coeff*q)


    print("OK")


if __name__ == "__main__":

    if argv.gen:
        gen()
    elif argv.main:
        main()
    else:
        test()



