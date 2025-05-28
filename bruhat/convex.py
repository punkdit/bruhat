#!/usr/bin/env python

from operator import mul, matmul, add
from functools import reduce

from sage.all_cmdline import ZZ, QQ

#from bruhat.matrix_sage import Matrix
from bruhat import matrix_sage
from bruhat.argv import argv


ring = QQ
one = ring.one()

def Matrix(v):
    return matrix_sage.Matrix(ring, v)

def identity(n):
    return matrix_sage.Matrix.identity(ring, n)

def make_free(n):
    if not n:
        return []
    r = one / n
    vs = []
    for i in range(n):
        v = []
        for j in range(n):
            if j==i:
                v.append((n-1)*one/n)
            else:
                v.append((-1)*one/n)
        assert sum(v) == 0
        v = Matrix(v).t
        vs.append(v)
    assert reduce(add, vs).sum() == 0
    return vs


def conv(a, b, r=one/2):
    return r*a + (1-r)*b


#def stack(vs):
#    u = None
#    for v in vs:


class Conv:
    def __init__(self, gens):
        self.gens = list(gens)
        # put gens in the columns of u:
        u = None
        for v in gens:
            u = v if u is None else u.augment(v)
        self.u = u
        shape = None
        for v in gens:
            assert v.sum() == 0
            assert shape is None or v.shape == shape
            shape = v.shape

    def __getitem__(self, idx):
        return self.gens[idx]

    def __len__(self):
        return len(self.gens)

    def __str__(self):
        #s = str(self.gens)
        s = str(self.u)
        #s = s.replace("\n", "")
        #s = s.replace(" ", "")
        return "Conv(\n%s, dim=%d)"%(s, self.dim)

    @property
    def dim(self):
        return self.u.t.rank()

    @classmethod
    def free(cls, n):
        vs = make_free(n)
        return Conv(vs)

    def _quotient1(self, l, r, rels):
        u = l-r
        uu = u@u.d       
        K = uu.cokernel()
        #print("_quotient1:")
        #print(K)
        #print(self)
        gens = [K*v for v in self.gens]
        rels = [(K*l, K*r) for (l,r) in rels]
        return Conv(gens), rels

    def quotient(self, *rels):
        #print("\nquotient:")
        other = self
        rels = list(rels)
        while rels:
            rel = rels.pop(0)
            other, rels = other._quotient1(*rel, rels)
            #print(other, len(rels))
        return other

    def _fail_quotient(self, *rels):
        # argh, how to do this in one go?
        u = None
        for (l,r) in rels:
            lr = l-r
            op = lr@lr.d 
            u = op if u is None else u.stack(op)
        K = u.cokernel()
        print("quotient:")
        print(u)
        print(K)
        gens = [K*v for v in self.gens]
        return Conv(gens)


def test():

    v = Matrix([1,0,0]).t

    a, b, c = make_free(3)
    assert (a+b+c).sum() == 0
    
    assert conv(a, b) == (one/2)*(a+b)
    assert conv(a, b, one/3) == (one/3)*(a+2*b)


    vs = make_free(4)
    a, b, c, d = vs
    u = conv(a,b,one/3) - conv(c, d)

    uu = ( u@u.d )
    assert uu.rank() == 1

    K = uu.cokernel()

    ws = [K*v for v in vs]
    a1, b1, c1, d1 = ws

    assert conv(a1,b1,one/3) == conv(c1, d1)

    for w in ws:
        assert w.sum() == 0

    # -------------------------------------------------------

    space = Conv.free(4)
    assert space.dim == 3

    # now construct the square by imposing one relator
    a, b, c, d = space
    other = space.quotient( (conv(a,b), conv(c,d) ) )
    a, b, c, d = other
    assert conv(a,b) == conv(c,d)
    assert other.dim == 2

    a, b, c, d = space
    other = space.quotient( (conv(a,b,one/3), conv(c,d) ) )
    a, b, c, d = other
    assert conv(a,b,one/3) == conv(c,d)
    assert other.dim == 2

    a, b, c, d = space
    other = space.quotient( (conv(a,b), conv(c,d)), (a, d) )
    assert other.dim == 1
    a, b, c, d = other
    assert a==d

    #return

    # now we construct the tensor product of the square with
    # itself, giving an 8 dimensional conv space
    space = Conv.free(16)
    lookup = {(i,j):space[i+4*j] for i in range(4) for j in range(4)}

    rels = []
    for i in range(4):
        a, b, c, d = [lookup[i,j] for j in range(4)]
        rels.append( (conv(a,b), conv(c,d)) )

        a, b, c, d = [lookup[j,i] for j in range(4)]
        rels.append( (conv(a,b), conv(c,d)) )

    other = space.quotient(*rels)
    assert other.dim == 8

    # -------------------------------------------------------

    space = Conv.free(6)
    a, b, c, d, e, f = space

    # construct an Octahedron == Bloch Octahedron
    octa = space.quotient( 
        (conv(a,b), conv(c,d)),
        (conv(a,b), conv(e,f)),
    )
    assert octa.dim == 3

    # tensor square of this is the separable states
    N = 6
    space = Conv.free(N*N)
    lookup = {(i,j):space[i+N*j] for i in range(N) for j in range(N)}

    rels = []
    for i in range(N):
        a, b, c, d, e, f = [lookup[i,j] for j in range(N)]
        rels.append( (conv(a,b), conv(c,d)) )
        rels.append( (conv(a,b), conv(e,f)) )

        a, b, c, d, e, f = [lookup[j,i] for j in range(N)]
        rels.append( (conv(a,b), conv(c,d)) )
        rels.append( (conv(a,b), conv(e,f)) )

    assert len(rels) == 24
    other = space.quotient(*rels)
    assert other.dim == 15

    print(len(space.gens))
    print(len(set(space.gens)))

    # -------------------------------------------------

    names = [(pm, a, b) for pm in (1,-1) for a in "XYZ" for b in "XYZ"]
    print(names)
    N = len(names)

    space = Conv.free(N)
    lookup = {name:space[i] for (i,name) in enumerate(names)}

    mixed = []
    for a in "XYZ":
      for b in "XYZ":
        l = lookup[(+1, a, b)]
        r = lookup[(-1, a, b)]
        mixed.append( conv(l,r) )
    l = mixed[0]
    rels = [(l,r) for r in mixed[1:]]

    other = space.quotient(*rels)
    print(other.dim)



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
        fn = eval(name)
        fn()

    else:
        test()


    t = time() - start_time
    print("OK! finished in %.3f seconds\n"%t)




