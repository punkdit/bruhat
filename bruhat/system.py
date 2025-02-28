#!/usr/bin/env python


import numpy
int_scalar = object

from sage.all_cmdline import (
    FiniteField, CyclotomicField, latex, block_diagonal_matrix,
    PolynomialRing, ZZ, QQ)
from sage import all_cmdline


from bruhat.matrix_sage import Matrix
from bruhat.argv import argv

def dot(A, B):
    return numpy.dot(A, B)

def identity(m):
    return numpy.identity(m, dtype=int_scalar)

def array(A):
    return numpy.array(A, dtype=int_scalar)

def zeros(m, n):
    return numpy.zeros((m,n), dtype=int_scalar)



class LinearCombination(object):
    def __init__(self, data={}):
        self.data = dict(data) # map ob to coefficient
        for key in data:
            assert type(key) is tuple or key==1

    @classmethod
    def promote(cls, item):
        if type(item) is cls:
            lin = item
        elif type(item) in (int, float): #, complex, numpy.complex128): # wah...
            if item==0:
                lin = cls()
            else:
                lin = cls({1 : item})
        elif type(item) is tuple:
            lin = cls({item : 1})
        else:
            raise TypeError(item, type(item))
        return lin

    def __str__(self):
        if not self.data:
            return "0"
        #s = ["%d*(%d,%d)" % (value, key[0], key[1]) 
        #    for (key, value) in self.data.items()]
        s = ["%s*%s" % (value, key)
            for (key, value) in list(self.data.items())]
        s = "+".join(s)
        s = s.replace("1*(", "(")
        return s
        #return "lin(%s)"%(self.data)
    __repr__ = __str__

    def __add__(self, other):
        other = self.promote(other)
        keys = set(list(self.data.keys()) + list(other.data.keys()))
        data = {}
        for key in keys:
            value = self.data.get(key, 0) + other.data.get(key, 0)
            if value != 0:
                data[key] = value
        return self.__class__(data)
    __radd__ = __add__

    def __neg__(self):
        return 0-self

    def __pos__(self):
        return self

    def __sub__(self, other):
        other = self.promote(other)
        keys = set(list(self.data.keys()) + list(other.data.keys()))
        data = {}
        for key in keys:
            value = self.data.get(key, 0) - other.data.get(key, 0)
            if value != 0:
                data[key] = value
        return self.__class__(data)

    def __rsub__(self, other):
        other = self.promote(other)
        keys = set(list(self.data.keys()) + list(other.data.keys()))
        data = {}
        for key in keys:
            value = other.data.get(key, 0) - self.data.get(key, 0)
            if value != 0:
                data[key] = value
        return self.__class__(data)

    def __mul__(self, other):
        #other = self.promote(other)
        #if type(other) in (int, long):
        if other==0:
            return self.__class__()
        if other==1:
            return self
        data = {}
        for key, value in list(self.data.items()):
            data[key] = other*value
        return self.__class__(data)
    __rmul__ = __mul__

    def __mod__(self, i):
        data = dict((key, value%i) for (key, value) in list(self.data.items()))
        return self.__class__(data)

    def __eq__(self, other):
        other = self.promote(other)
        return self.data == other.data

    def __ne__(self, other):
        return not self==other



def MatrixEntry(i, j):
    return LinearCombination({(i, j) : 1})


def Unknown(*args):
    if len(args)==2:
        m, n = args
    elif len(args)==1:
        m, n = args[0]
    else:
        raise TypeError
    data = [[MatrixEntry(i, j) for j in range(n)] for i in range(m)]
    A = numpy.array(data)
    return A


class System(object):
    "system of linear equations"
    def __init__(self, ring, *Ts):
        assert Ts
        self.ring = ring
        m = sum(T.shape[0] for T in Ts)
        n = sum(T.shape[1] for T in Ts)
        T = Unknown(m, n)
        self.lhss = []
        self.rhss = []
        i = 0
        blocks = []
        for T0 in Ts:
            assert T0.dtype == object # the Unknown
            j = 0
            row = []
            for T1 in Ts:
                row.append( [[i, i+T1.shape[0]], [j, j+T1.shape[1]]] )
                j += T1.shape[1]
            blocks.append(row)
            i += T0.shape[0]
        self.blocks = blocks = numpy.array(blocks, dtype=int_scalar)
        n = self.n = len(Ts)
        for i in range(n):
          for j in range(n):
            blk = blocks[i, j]
            T1 = T[blk[0, 0]:blk[0, 1], blk[1, 0]:blk[1, 1]]
            if i==j:
                Ts[i][:] = T1
            else:
                self.append(T1, 0)
        self.Ts = Ts
        self.T = T

    def append(self, lhs, rhs):
        if type(rhs) in (int, float) and rhs==0:
            rhs = zeros(*lhs.shape)
        if type(rhs) in (int, float) and rhs==1:
            n = min(*lhs.shape)
            rhs = identity(n)
#        if lhs.dtype==int_scalar:
#            lhs, rhs = rhs, lhs # swap
#        elif lhs.dtype==object and rhs.dtype==object:
#            lhs = lhs-rhs
#            rhs = zeros(*rhs.shape)
#        assert lhs.dtype == object, lhs.dtype
#        assert rhs.dtype == int_scalar, rhs.dtype
#        if len(lhs.shape)==1:
#            lhs = lhs.view()
#            lhs.shape = (lhs.shape[0], 1)
#        if len(rhs.shape)==1:
#            rhs = rhs.view()
#            rhs.shape = (rhs.shape[0], 1)
        # must have the unknowns in the lhs, rhs is constant
        assert lhs.shape == rhs.shape, (lhs.shape, rhs.shape)
        self.lhss.append(lhs)
        self.rhss.append(rhs)

    def build(self):
        T, lhss, rhss = self.T, self.lhss, self.rhss
        m = sum(lhs.shape[0]*lhs.shape[1] for lhs in lhss)
        n = T.shape[0]*T.shape[1]
        H = zeros(
            m, # this many constraints
            n) # this many unknowns
        v = zeros(m, 1)
        row = 0
        for lhs, rhs in zip(lhss, rhss):
          #print(lhs)
          #print(rhs)
          for i0 in range(lhs.shape[0]):
            for j0 in range(lhs.shape[1]):
                lin = lhs[i0, j0]
                #lin = LinearCombination.promote(lin)
                assert isinstance(lin, LinearCombination), repr(lhs)
                #print(lin.data)
                v[row] = rhs[i0, j0]
                #for (i1, j1), w in list(lin.data.items()):
                for key,val in lin.data.items(): 
                    if type(key) is tuple:
                        i1,j1 = key
                        #print(key, val)
                        H[row, T.shape[1]*i1+j1] = val
                    else:
                        assert 0
                        assert key == 1
                        v[row] -= val
                row += 1
        return H, v

    def kernel(self):
        H, v = self.build()
        #basis = find_kernel(H)
        H = Matrix(self.ring, H)
        basis = H.kernel().t
        return basis

    def solve(self, unpack=True, check=False):
        H, v = self.build()
        #print("solve:")
        #print(H)
        H = Matrix(self.ring, H)
        v = Matrix(self.ring, v)
        #Hinv = ~H
        Hinv = H.pseudoinverse()
        #print("Hinv:")
        #print(Hinv)
        #u = dot(Hinv, v)
        u = Hinv * v
        #print shortstrx(H, u, v)
        if H*u != v:
            # No solution
            return None
        #if not eq(dot(H, u), v):
            #return None
        #u.shape = self.T.shape
        u = u.reshape(*self.T.shape)
        #print shortstrx(T)
        if unpack:
            u = self.unpack(u)
        return u

    def unpack(self, u):
        if len(self.Ts)>1:
            us = []
            i, j = 0, 0
            for T in self.Ts:
                u0 = u[i:i+T.shape[0], j:j+T.shape[1]]
                i += T.shape[0]
                j += T.shape[1]
                us.append(u0)
            return us
        else:
            return u

    def all_solutions(self):
        u = self.solve(unpack=False)
        kern = self.kernel()
        if not kern:
            return
        #print("kern:")
        #print(kern)
        #for v in span(kern):
        for v in kern:
            #v.shape = u.shape
            v = v.reshape(*u.shape)
            yield self.unpack((u+v))

    def solve_homogeneous(self):
        kern = self.kernel()
        if not kern:
            return
        for v in kern:
            #v = v.reshape(*u.shape)
            yield self.unpack(v)

    def random_solution(self):
        u = self.solve(unpack=False)
        kern = self.kernel()
        if not kern:
            return
        for v in kern:
            v.shape = u.shape
            if random()<=0.5:
                u += v
        return self.unpack(u)


def test():
    #H = Matrix(QQ, [[1,0,1],[1,0,0],[0,1,1]])
    #print(H * ~H)

    one = QQ.gen()
    r = one/2
    H = array([[r,0,1],[1,0,0],[0,1,1]])

    m, n = H.shape

    T = Unknown(n, m)
    print(T)
    print(dot(H,T))

    system = System(QQ, T)
    system.append(dot(H, T), identity(m))

    T = system.solve()
    print("T")
    print(T)

    H = Matrix(QQ, H)
    assert H*T == Matrix.identity(QQ, m)




if __name__ == "__main__":

    from time import time
    start_time = time()

    profile = argv.profile
    name = argv.next() or "test"
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        print("\n%s()"%name)
        fn = eval(name)
        fn()

    else:
        test()


    t = time() - start_time
    print("OK! finished in %.3f seconds\n"%t)




