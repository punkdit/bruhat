#!/usr/bin/env python

"""
Applying distributivity to computing the determinant of a matrix.
What is the optimal solution ?
How close is Dodgson condensation to being optimal?

"""

import sys, os
from operator import add, mul
from random import shuffle, random

#import numpy

from bruhat.util import cross, allsignedperms


class Op(object):
    "commutative mult-arity operation"

    def __init__(self, items, check=True):
        items = list(items)
        items.sort()
        self.items = items
        if check:
            self.check()

    def check(self):
        pass

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self.items == other.items

    def __ne__(self, other):
        if self.__class__ != other.__class__:
            return True
        return self.items != other.items

#    def __le__(self, other):
#        assert 0
#        return str(self) <= str(other)

    def __cmp__(self, other):
        #return cmp(str(self), str(other))
        i = cmp(self.__class__.order, other.__class__.order)
        if i!=0:
            return i
        i = cmp(len(self), len(other))
        if i!=0:
            return i
        i = cmp(self.items, other.items) # lexicographic
        return i

    def __hash__(self):
        key = (self.__class__,) + tuple(self.items)
        return hash(key)

    def __repr__(self):
        s = "%s(%s)"%(
            self.__class__.__name__, 
            ', '.join(repr(item) for item in self.items))
        return s

    def __str__(self):
        s = self.str()
        s = s.replace("+ -", "- ")
        s = s.replace(" 1*", " ")
        s = s.replace("(1*", "(")
        return s

    def allvars(self):
        if isinstance(self, Var):
            return [self[0]]
        elif isinstance(self, Integer):
            return []
        vs = reduce(add, (op.allvars() for op in self))
        vs = list(set(vs))
        vs.sort()
        return vs

    def getlambda(self):
        vs = self.allvars()
        s = "lambda %s : %s" % (', '.join(vs), self)
        f = eval(s)
        return f

    def str(self):
        return self.__repr__()

    def __getitem__(self, i):
        return self.items[i]

    def __len__(self):
        return len(self.items)

    def __add__(self, other):
        if self.__class__==Add and other.__class__==Add:
            op = Add(self.items + other.items)
        elif self.__class__==Add:
            op = Add(self.items + [other])
        elif other.__class__==Add:
            op = Add([self] + other.items)
        else:
            op = Add([self] + [other])
        return op

    def __radd__(self, other):
        if isinstance(other, (int, long)):
            return Integer(other)+self
        raise TypeError

    def __mul__(self, other):
        if self.__class__==Mul and other.__class__==Mul:
            items = self.items + other.items
        elif self.__class__==Mul:
            items = self.items + [other]
        elif other.__class__==Mul:
            items = [self] + other.items
        else:
            items = [self] + [other]
        counts = {}
        for item in items:
            n = 1
            while isinstance(item, Reciprocal):
                n *= -1
                item = item[0]
            counts[item] = counts.get(item, 0) + n
        ops = []
        for op, n in counts.items():
            if n == 0:
                continue
            if n < 0:
                n *= -1
                op = Reciprocal(op)
            for i in range(n):
                ops.append(op)
        op = Mul(ops)
        return op

    def __rmul__(self, other):
        if isinstance(other, (int, long)):
            return Integer(other)*self
        raise TypeError

    def __sub__(self, other):
        return self + Integer(-1)*other

    def __div__(self, other):
        return self * Reciprocal(other)

    def __pow__(self, n):
        if n==0:
            return Mul([])
        if n<0:
            return Reciprocal(self.__pow__(-n))
        op = self
        while n>1:
            op = op*self
            n -= 1
        return op

    def dist(self):
        return self

    def simplify(self):
        return self

    def equal(lhs, rhs):
        # distribute * over + to make normal form
        lhs = lhs.dist()
        rhs = rhs.dist()
        return lhs==rhs

    def opcount(self):
        return 0

    def factorize_1(self):
        return self

    def factorize(self):
        op = self
        while 1:
            op1 = op.factorize_1()
            if op1 == op:
                break
            op = op1
        return op


class Var(Op):
    order = 1 # unique per-class
    def __init__(self, value):
        Op.__init__(self, [value])

    def str(self):
        return str(self.items[0])


class Reciprocal(Op):
    order = 4

    def __init__(self, value):
        #if isinstance(value, Reciprocal): XX
        Op.__init__(self, [value])

    def str(self):
        return "(1./(%s))"%(self[0])


class Integer(Op):
    order = 0 # unique per-class
    def __init__(self, value):
        assert long(value) == value
        Op.__init__(self, [value])

    def str(self):
        return "%s"%(self.items[0],)

    def __cmp__(self, other):
        if isinstance(other, Integer):
            return -cmp(self[0], other[0])
        else:
            return -cmp(other, self)

    def __add__(self, other):
        if isinstance(other, Integer):
            return Integer(self[0] + other[0])
        else:
            return other + self

    def __mul__(self, other):
        if isinstance(other, Integer):
            return Integer(self[0] * other[0])
        else:
            return other * self

    def __div__(self, other):
        if isinstance(other, Integer):
            assert self[0] % other[0] == 0, (self, other)
            return Integer(self[0] // other[0])
        else:
            return self * Reciprocal(other)
        


class Add(Op):
    order = 2 # unique per-class
    op = "+"

    def str(self):
        return "(%s)"%(' + '.join(item.str() for item in self.items))

    def check(self):
        for item in self:
            if item.__class__ is Add:
                assert 0, repr(self)

    def dist(self):
        "depth-first _apply dist"
        if not self:
            return self
        op = self[0].dist()
        i = 1
        while i < len(self):
            op = op + self[i].dist()
            i += 1
        return op

    def simplify(self):
        if not self:
            return self
        ops = [op.simplify() for op in self]
        counts = {}
        for op in ops:
            counts[op] = counts.get(op, 0) + 1
        ops = []
        for key, val in counts.items():
            if val==0:
                continue
            elif val==1:
                ops.append(key)
            else:
                ops.append(Integer(val)*key)
        if len(ops)>1:
            op = reduce(add, ops)
        elif len(ops)==0:
            op = Integer(0)
        else:
            op = ops[0]
        return op

    def opcount(self):
        return reduce(add, [op.opcount() for op in self])

    def __XX__factorize(self):
        opss = []
        remain = []
        for op in self:
            if not isinstance(op, Mul):
                remain.append(op)
                continue
            ops = set(op)
            opss.append(ops)
        n = len(opss)
        pairs = {}
        for i in range(n):
          for j in range(i+1, n):
            ops = opss[i].intersection(opss[j])
            ops = [op for op in ops if not isinstance(Integer)]
            if not ops:
                continue
            pairs[i, j] = ops

    def factorize_1(self):
        ops = [op.factorize() for op in self]
        #shuffle(ops)
        n = len(ops)
        for i in range(n):
            if ops[i].__class__ != Mul:
                continue
            for j in range(i+1, n):
                if ops[j].__class__ != Mul:
                    continue
                op = ops[i].gdc(ops[j])
                if op is None:
                    continue
                opj = ops.pop(j) # pop j first
                opi = ops.pop(i)
                op = op*(opi.div(op) + opj.div(op))
                if len(ops):
                    return reduce(add, ops+[op])
                else:
                    return op
        return self



class Mul(Op):
    order = 3 # unique per-class
    op = "*"

    def str(self):
        return ('*'.join(item.str() for item in self.items))

    def check(self):
        for item in self:
            if item.__class__ is Mul:
                assert 0, repr(self)

    def dist(self):
        "depth-first _apply dist"
        if not self:
            return self
        op = self[0].dist()
        i = 1
        while i < len(self):
            op = op * self[i].dist()
            i += 1
        itemss = []
        for item in self:
            if isinstance(item, (Var, Integer, Reciprocal)):
                itemss.append([item])
            elif isinstance(item, Add):
                itemss.append(item.items)
            else:
                assert 0, repr(item)
        ops = []
        for items in cross(itemss):
            ops.append(reduce(mul, items))
        assert len(ops)
        if len(ops)==1:
            op = ops[0]
        else:
            op = Add(ops)
        return op

    def factorize_1(self):
        ops = [op.factorize() for op in self]
        op = reduce(mul, ops)
        return op

    def opcount(self):
        ops = [op for op in self if not isinstance(op, (Integer, Reciprocal))]
        count = max(0, len(ops)-1)
        for op in self:
            count += op.opcount()
        return count

    def gdc(self, other):
        if not isinstance(other, Mul):
            other = Mul([other])
        a_counts = {}
        for op in self:
            if isinstance(op, Integer):
                continue
            a_counts[op] = a_counts.get(op, 0) + 1
        b_counts = {}
        for op in other:
            b_counts[op] = b_counts.get(op, 0) + 1
        ops = []
        for op in a_counts.keys():
            n = min(a_counts[op], b_counts.get(op, 0))
            if n:
                ops += [op]*n
        if not ops:
            return
        if len(ops)==1:
            return ops[0]
        return reduce(mul, ops)

    def div(self, other):
        if not isinstance(other, Mul):
            other = Mul([other])
        a_counts = {}
        for op in self:
            a_counts[op] = a_counts.get(op, 0) + 1
        b_counts = {}
        for op in other:
            b_counts[op] = b_counts.get(op, 0) + 1
            assert op in a_counts
        ops = []
        for op in a_counts.keys():
            n = a_counts[op] - b_counts.get(op, 0)
            assert n>=0
            if n:
                ops += [op]*n
        if not ops:
            return
        if len(ops)==1:
            return ops[0]
        return reduce(mul, ops)


class Matrix(object):
    def __init__(self, rows):
        self.rows = list(list(row) for row in rows)
        self.m = len(rows)
        self.n = len(rows[0])
        for row in rows:
            assert len(row) == self.n

    def __getitem__(self, key):
        i, j = key
        row = self.rows[i]
        val = row[j]
        return val

    def __setitem__(self, key, val):
        i, j = key
        row = self.rows[i]
        row[j] = val

    def condense(self):
        assert self.m == self.n
        assert self.n > 1
        n = self.n
        rows = []
        for i in range(n-1):
            row = []
            for j in range(n-1):
                x = self[i, j]*self[i+1, j+1] - self[i+1, j] * self[i, j+1]
                row.append(x)
            rows.append(row)
        return Matrix(rows)

    def dodgson(self):
        # https://en.wikipedia.org/wiki/Dodgson_condensation
        #print "dodgson"
        A = self
        B = A.condense()
        while B.n > 1:
            #print "while:"
            #print "\tA =", A.rows
            #print "\tB =", B.rows
            assert B.n == A.n-1
            C = B.condense()
            #print "\tC =", C.rows
            for i in range(C.n):
              for j in range(C.n):
                C[i, j] = C[i, j] / A[i+1, j+1]
            #print "\tC =", C.rows
            A, B = B, C
        return B[0, 0]

    def determinant(self):
        assert self.m == self.n
        assert self.m > 0
        m = self.m

        vals = []
        idxs = range(m)
        for sign, idxs in allsignedperms(range(m)):
            ops = [Integer(sign)]
            for row, col in enumerate(idxs):
                ops.append(self[row, col])
            #op = Mul(ops)
            op = reduce(mul, ops)
            vals.append(op)
        #return Add(vals)
        return reduce(add, vals)


def matrix(n):
    rows = [[Var("x%d%d"%(i+1, j+1)) for j in range(n)] for i in range(n)]
    return Matrix(rows)


def evaleq(p, q):
    vs = p.allvars()
    assert q.allvars() == vs
    f = p.getlambda()
    g = q.getlambda()
    vs = [0.1 + random() for v in vs]
    error = abs(f(*vs) - g(*vs))
    return error < 1e-8


def main():

    x, y, z = [Var(_) for _ in 'xyz']
    a, b, c, d = [Var(_) for _ in 'abcd']

    #assert str(x+y) == "(x+y)"
    assert x==x
    assert x!=z
    assert x+y == x+y
    assert x+y == y+x
    assert x+y != x+z

    #print (a*c + a*b*c + b*d)
    A = (a*c + a*b*c + b*d)
    assert A==A
    B = A.dist()

    assert A.__class__ is B.__class__
    assert len(A) == len(B)
    for i in range(len(A)):
        assert A[i] == B[i], (A[i], B[i])
    assert A.items == B.items
    assert A == B

    assert (a*(b+c)).dist() == a*b + a*c
    assert ((a+b)*(c+d)).dist() == a*c + a*d + b*c + b*d

    assert (a*b*c).gdc(d) == None
    assert (a*b*c).gdc(b) == b
    assert (a*b*c).gdc(b*c*c) == b*c
    assert (a*b*c).div(b*c) == a

    p = a*b + a*c
    assert p.factorize() == a*(b+c)

    p = a*(c + d) + b*(c + d)
    assert p.factorize() == (a+b) * (c+d)

    p = a*c+a*d+b*c+b*d
    assert p.factorize() == (a+b) * (c+d)

    A = A*A
    A.dist()
    A = ((a+b)**4).dist()
    
    rows = [
        [-2 , -1 , -1 , -4 ],
        [-1 , -2 , -1 , -6 ],
        [-1 , -1 , 2 , 4 ],
        [2 , 1 , -3 , -8]]
    rows = [[Integer(i) for i in row] for row in rows]
    A = Matrix(rows)
    assert A.determinant() == Integer(-8)
    assert A.dodgson() == Integer(-8)

    a, b, c, d, e, f, g, h, i = [Var(_) for _ in 'abcdefghi']
    A = Matrix([[a,b,c],[d,e,f],[g,h,i]])

    for n in range(3, 6):
        print "n =", n
        A = matrix(n)
        p = A.determinant()
        #print p
        print p.opcount()
        q = p.factorize()
        #print q
        print q.opcount()
    
        r = A.dodgson()
        print r.opcount()

        assert evaleq(p, p)
        assert evaleq(p, q)
        assert evaleq(p, r)

    print "n = 6"
    A = matrix(6)
    p = A.determinant()
    print p.opcount()
    r = A.dodgson()
    print r.opcount()






if __name__ == "__main__":

    main()

