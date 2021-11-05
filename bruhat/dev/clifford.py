#!/usr/bin/env python3

import numpy
from numpy import dot, alltrue, zeros, array, identity


from bruhat.action import mulclose_fast
from bruhat.smap import SMap


#int_scalar = numpy.int8 # dangerous...
int_scalar = numpy.int32 # less dangerous ?


class ___Clifford(object):

    def __init__(self, A, dyadic=0):
        assert len(A.shape)==3
        m, n, k = A.shape
        assert k==2
        assert m==n
        self.A = A
        self.dyadic = dyadic

    def __str__(self):
        smap = SMap()
        A = self.A
        N = A.shape[0]
        smap[0, 0] = '['
        col = lambda j : 2+3*j
        for i in range(N):
          smap[i, 1] = '['
          for j in range(N):
            re, im = A[i, j, :]
            if re == 0 and im == 0:
                s = '.'
            elif im == 0:
                s = str(re)
            elif re == 0 and im==1:
                s = "i"
            elif re == 0 and im==-1:
                s = "-i"
            elif re == 0:
                assert 0, im
            else:
                s = "%d+%di"%(re, im)
            smap[i, col(j)] = s
          smap[i, col(N)] = ']'
        smap[i, 1+col(N)] = ']'
        return str(smap)

    @classmethod
    def identity(cls, n):
        N = 2**n
        A = numpy.zeros((N, N, 2), dtype=int_scalar)
        for i in range(N):
            A[i, i, 0] = 1
        return cls(A)


def normalize(re, im, dyadic):
    while dyadic < 0:
        if not alltrue(re & 1 == 0) or not alltrue(im & 1 == 0):
            break
        re = re >> 1
        im = im >> 1
        dyadic += 1
    return re, im, dyadic


class Clifford(object):

    def __init__(self, re, im, dyadic=0, inv=None):
        assert dyadic <= 0, self
        if dyadic < 0:
            re, im, dyadic = normalize(re, im, dyadic)
        self.re = re
        self.im = im
        self.dyadic = dyadic
        self.inv = inv

    def __str__(self):
        dyadic = self.dyadic
        s = "%s + \ni%s" % (self.re, self.im)
        if dyadic != 0:
            i = 2**dyadic
            s = "%s(%s)"%(i, s)
        return s
    __repr__ = __str__

    def __eq__(self, other):
        return (self.dyadic == other.dyadic and
            alltrue(self.re == other.re) and 
            alltrue(self.im == other.im))

    def __hash__(A):
        re, im, dyadic = A.re, A.im, A.dyadic
        key = re.tobytes(), im.tobytes(), dyadic
        return hash(key)

    def mul(A, B, inv=None):
        re = dot(A.re, B.re) - dot(A.im, B.im)
        im = dot(A.re, B.im) + dot(A.im, B.re)
        dyadic = A.dyadic + B.dyadic
        op = Clifford(re, im, dyadic)
        if inv is None and B.inv is not None and A.inv is not None:
            inv = B.inv.mul(A.inv, op) # recurse
        op.inv = inv
        return op
    __mul__ = mul

    def __neg__(A):
        re, im, dyadic = A.re, A.im, A.dyadic
        op = Clifford(-re, -im, dyadic)
        op.inv = -A.inv if A.inv is not A else op
        return op

    @classmethod
    def identity(cls, n):
        N = 2**n
        re = identity(N, dtype=int_scalar)
        im = zeros((N, N), dtype=int_scalar)
        I = cls(re, im)
        I.inv = I
        return I

    @classmethod
    def phase(cls, n):
        N = 2**n
        re = zeros((N, N), dtype=int_scalar)
        im = identity(N, dtype=int_scalar)
        op = cls(re, im)
        op.inv = cls(-re, -im, 0, op)
        return op

    @classmethod
    def zgate(cls):
        N = 2
        re = identity(N, dtype=int_scalar)
        re[1, 1] = -1
        im = zeros((N, N), dtype=int_scalar)
        op = cls(re, im)
        op.inv = op
        return op

    @classmethod
    def sgate(cls):
        N = 2
        re = zeros((N, N), dtype=int_scalar)
        im = zeros((N, N), dtype=int_scalar)
        re[0, 0] = 1
        im[1, 1] = 1
        op = cls(re, im)
        re = zeros((N, N), dtype=int_scalar)
        im = zeros((N, N), dtype=int_scalar)
        re[0, 0] = 1
        im[1, 1] = -1
        op.inv = cls(re, im, 0, op)
        return op

    @classmethod
    def xgate(cls):
        N = 2
        re = zeros((N, N), dtype=int_scalar)
        im = zeros((N, N), dtype=int_scalar)
        re[1, 0] = 1
        re[0, 1] = 1
        op = cls(re, im)
        op.inv = op
        return op

    @classmethod
    def ygate(cls):
        N = 2
        re = zeros((N, N), dtype=int_scalar)
        im = zeros((N, N), dtype=int_scalar)
        im[1, 0] = 1
        im[0, 1] = -1
        op = cls(re, im)
        op.inv = op
        return op

    @classmethod
    def hgate(cls):
        N = 2
        re = array([[1, -1], [1, 1]], dtype=int_scalar)
        im = array([[-1, +1], [-1, -1]], dtype=int_scalar)
        op = cls(re, im, -1)
        re = array([[1, 1], [-1, 1]], dtype=int_scalar)
        im = array([[+1, +1], [-1, +1]], dtype=int_scalar)
        op.inv = cls(re, im, -1, op)
        return op




def main():

    I = Clifford.identity(1)
    iI = Clifford.phase(1)
    nI = iI * iI
    X = Clifford.xgate()
    Z = Clifford.zgate()
    Y = Clifford.ygate()
    S = Clifford.sgate()

    assert I != X != Z != S
    assert X*X == I
    assert Z*Z == I
    assert Y*Y == I
    assert X*Z != Z*X
    assert X*Z == nI*Z*X
    assert S*S == Z
    assert S*S*S*S == I
    assert S*X*S*X == iI

    H = Clifford.hgate()
    assert H*H == -Y
    assert H*H*H*H == I

    gen = [H, S, X]
    G = mulclose_fast(gen)
    assert len(G) == 96

    for g in G:
        assert g * g.inv == I
    



if __name__ == "__main__":

    main()


