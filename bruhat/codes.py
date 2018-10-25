#!/usr/bin/env python3

"""
Build Reed-Muller codes.
"""


import numpy

from solve import array2, zeros2, dot2, shortstr, rank, find_kernel, span
from solve import linear_independent, parse, pseudo_inverse, eq2, rand2
from argv import argv
from util import choose, cross


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
        return is_morthogonal(self.G, m)

    def is_triorthogonal(self):
        return is_morthogonal(self.G, 3)

    def weight_enum(self, A, B):
        G = self.G
        the_op = None
        for v in span(G):
            op = A if v[0]==0 else B
            for vi in v[1:]:
                op = op * (A if vi==0 else B)
            the_op = op if the_op is None else the_op + op
        return the_op

    def tensor_enum(self, A, B):
        G = self.G
        the_op = None
        for v in span(G):
            op = A if v[0]==0 else B
            for vi in v[1:]:
                op = op @ (A if vi==0 else B)
            the_op = op if the_op is None else the_op + op
        return the_op



def is_morthogonal(G, m):
    k = len(G)
    assert m>=2
    if m>2 and not is_morthogonal(G, m-1):
        return False
    items = list(range(k))
    for idxs in choose(items, m):
        v = G[idxs[0]]
        for idx in idxs[1:]:
            v = v * G[idx]
        if v.sum()%2 != 0:
            return False
    return True


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

    #print(code)
    #print("d =", code.get_distance())
    #code.dump()

    code = code.puncture(3)

    #print(code)
    #code.dump()
    #print("d =", code.get_distance())


    for m in range(2, 6):
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

        if p.is_triorthogonal():
            G = p.G
            #print(shortstr(G))
            A = list(span(G))
            A = array2(A)
            print(is_morthogonal(A, 3))


def test_triorth():
    code = reed_muller(1, 5)
    code = code.puncture(0)
    code.dump()

    print(code.is_triorthogonal())
    A = array2(list(span(code.G)))
    print(is_morthogonal(A, 2))

    k = len(A)

    for i in range(k):
      for j in range(i+1, k):
        u = A[i]
        v = A[j]
        x = (u*v).sum() % 2
        if x == 0:
            continue
        #print(shortstr(u))
        #print(shortstr(v))
        #print()

    for a in range(k):
      for b in range(a+1, k):
       for c in range(b+1, k):
        u = A[a]
        v = A[b]
        w = A[c]
        x = (u*v*w).sum() % 2
        #if x:
            #print(a, b, c)



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


def test_dual():

    from vec import Space, Hom, Map

    import element

    #ring = element.Z
    ring = element.Q
    one = ring.one

    space = Space(2, ring)
    hom = Hom(space, space)

    I = Map.from_array([[1, 0], [0, 1]], hom)
    X = Map.from_array([[0, 1], [1, 0]], hom)
    Z = Map.from_array([[1, 0], [0, -1]], hom)

    assert X+X == (2*X)
    assert X*Z == -Z*X

    assert (X@X) * (Z@Z) == (Z@Z) * (X@X)
    assert (I@X@X) * (I@Z@Z) == (I@Z@Z) * (I@X@X)

    if argv.code == "repitition":
        G = parse("111")
        H = parse("""
        11.
        .11""")
    elif argv.code == "steane":
        G = parse("""
        1111...
        .1.11.1
        ..11.11
        """)
    else:
        return

    code = Code(G)
    if argv.dual:
        code = code.get_dual()

    dual = code.get_dual()

    code.dump()

    #W = lambda A, B : (A@A@A + B@B@B)
    #WD = lambda A, B : (A@A@A + B@B@A + B@A@B + A@B@B)
    W = code.tensor_enum
    WD = dual.tensor_enum

    a = one/(2**len(code.G))
    b = one/(2**len(dual.G))
    A = WD(I, X)
    B = W(I, Z)
    assert A*B == B*A
    AA = W(I+X, I-X)
    BB = WD(I+Z, I-Z)
    assert a*AA == A
    assert b*BB == B
    assert a*b*AA*BB == A*B
    #print(W(I+X, I-X))
    #print(WD(I, X))
    #print(W(I, Z))

    print("OK")



def search():
    # Bravyi, Haah, 1209.2426v1 sec IX.
    # https://arxiv.org/pdf/1209.2426.pdf

    m = argv.get("m", 3)
    n = argv.get("n", 6)
    k = argv.get("k", None) # number of odd-weight rows

    # these are the variables N_x
    xs = list(cross([(0, 1)]*m))
    N = len(xs)

    lhs = []
    rhs = []

    # bi-orthogonality
    for a in range(m):
      for b in range(a+1, m):
        v = zeros2(N)
        for i, x in enumerate(xs):
            if x[a] == x[b] == 1:
                v[i] = 1
        if v.sum():
            lhs.append(v)
            rhs.append(0)

    # tri-orthogonality
    for a in range(m):
      for b in range(a+1, m):
       for c in range(b+1, m):
        v = zeros2(N)
        for i, x in enumerate(xs):
            if x[a] == x[b] == x[c] == 1:
                v[i] = 1
        if v.sum():
            lhs.append(v)
            rhs.append(0)

    # dissallow columns with weight <= 1
    for i, x in enumerate(xs):
        if sum(x)<=1:
            v = zeros2(N)
            v[i] = 1
            lhs.append(v)
            rhs.append(0)

    if k is not None:
      # constrain to k number of odd-weight rows
      assert 0<=k<m
      for a in range(m):
        v = zeros2(N)
        for i, x in enumerate(xs):
          if x[a] == 1:
            v[i] = 1
        lhs.append(v)
        if a<k:
            rhs.append(1)
        else:
            rhs.append(0)

    A = array2(lhs)
    rhs = array2(rhs)
    #print(shortstr(A))

    K = array2(list(find_kernel(A)))
    #print(K)
    #print( dot2(A, K.transpose()))
    #sols = []
    #for v in span(K):
    best = None
    density = 1.0
    trials = argv.get("trials", 1024)
    count = 0
    for trial in range(trials):
        u = rand2(len(K), 1)
        v = dot2(K.transpose(), u)
        assert dot2(A, v).sum()==0
        #if v.sum() != n:
        #    continue
        assert v[0]==0
        Gt = []
        for i, x in enumerate(xs):
            if v[i]:
                Gt.append(x)
        if not Gt:
            continue
        Gt = array2(Gt)
        G = Gt.transpose()
        assert is_morthogonal(G, 3)
        if G.shape[1]<m:
            continue

        #print(shortstr(G))
#        for g in G:
#            print(shortstr(g), g.sum())
#        print()

        _density = float(G.sum()) / (G.shape[0]*G.shape[1])
        if best is None or _density < density:
            best = G
            density = _density

        if 0:
            #sols.append(G)
            Gx = even_rows(G)
            assert is_morthogonal(Gx, 3)
            if len(Gx)==0:
                continue
            GGx = array2(list(span(Gx)))
            assert is_morthogonal(GGx, 3)

        count += 1

    print("found %d solutions" % count)

    G = best
    #print(shortstr(G))
    for g in G:
        print(shortstr(g), g.sum())
    print()
    print("density:", density)
    

    if 0:
        B = pseudo_inverse(A)
        v = dot2(B, rhs)
        print("B:")
        print(shortstr(B))
        print("v:")
        print(shortstr(v))
        assert eq2(dot2(B, v), rhs) 


def even_rows(G):
    Gx = []
    for u in G:
        #print(shortstr(u), u.sum()%2)
        parity = u.sum()%2
        if parity==0:
            Gx.append(u)
    Gx = array2(Gx)
    return Gx


def get_code():

    name = argv.get("code")
    code = None
    if name == "toric":
        G = parse("""
        1.1.....
        .1...1..
        11.11...
        .111..1.
        1...11.1
        """) # l=2 toric code X logops + X stabs 

    elif name == "toric3":
        G = parse("""
        .....1.....1.....1
        ............1.1.1.
        .111..........1...
        ...111..........1.
        1.....11...1......
        ..1....111........
        ....1....111......
        ......1.....11...1
        ........1....111..
        ..........1....111
        """) # l=3 toric code X logops + X stabs 

    elif name == "steane":
        G = parse("""
        1111111
        1111...
        .1.11.1
        ..11.11
        """)

    elif name == "haah":
        # triorthogonal matrix
        G = parse("""
        1111111.......
        .......1111111
        1.1.1.11.1.1.1
        .11..11.11..11
        ...1111...1111
        """)
        G = parse("""
        1.1.1.11.1.1.1
        .11..11.11..11
        ...1111...1111
        """)

        G = parse("""
        ..............1111111111111111
        ......11111111........11111111
        ..1111....1111....1111....1111
        .1..11..11..11..11..11..11..11
        11.1.1.1.1.1.1.1.1.1.1.1.1.1.1
        """)

    elif name == "rm":
        r = argv.get("r", 1)
        m = argv.get("m", 5)
        code = reed_muller(r, m)
        if argv.puncture:
            code = code.puncture(0)

    else:
        return None

    if code is None:
        code = Code(G)

    if argv.dual:
        code = code.get_dual()

    return code


def triortho():
    code = get_code()

    code.dump()
    print(code)

    Gx = []
    for u in code.G:
        print(shortstr(u), u.sum()%2)
        parity = u.sum()%2
        if parity==0:
            Gx.append(u)
    Gx = array2(Gx)

    print("is_triorthogonal:", code.is_triorthogonal())

    A = array2(list(span(Gx)))
    print("span(Gx) is_morthogonal(2):", is_morthogonal(A, 2))
    print("span(Gx) is_morthogonal(3):", is_morthogonal(A, 3))

    return

    G = code.G

#    A = array2(list(span(G)))
#    poly = {}
#    for v in A:
#        w = v.sum()
#        poly[w] = poly.get(w, 0) + 1
#    print(poly)

    k, n = G.shape

    if 0:
        from comm import Poly
        a = Poly({(1,0):1})
        b = Poly({(0,1):1})
        poly = Poly.zero(2)
        for v in span(G):
            w = v.sum()
            term = Poly({(n-w,0) : 1}) * Poly({(0,w) : 1})
            poly = poly + term
        print(poly)

    # print higher genus weight enumerator
    genus = argv.get("genus", 1)
    assert 1<=genus<=4
    N = 2**genus
    idxs = list(cross([(0,1)]*genus))

    cs = {} # _coefficients : map exponent to coeff
    for vs in cross([list(span(G)) for _ in range(genus)]):
        key = [0]*N
        for i in range(n):
            ii = tuple(v[i] for v in vs)
            idx = idxs.index(ii)
            key[idx] += 1
        key = tuple(key)
        cs[key] = cs.get(key, 0) + 1
    #print(cs)
    keys = list(cs.keys())
    keys.sort()
    print(idxs)
    for key in keys:
        print(key, cs[key])
            


if __name__ == "__main__":

    name = argv.next()
    if name is None:
        test()
    else:
        fn = eval(name)
        fn()



