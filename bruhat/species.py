#!/usr/bin/env python3

"""
Combinatorial species.

See chapter 4 here:
    COMBINATORIAL SPECIES AND LABELLED STRUCTURES 
    Brent Abraham Yorgey 
    https://www.cis.upenn.edu/~sweirich/papers/yorgey-thesis.pdf
"""

from types import GeneratorType
import string
letters = list(string.ascii_lowercase)

from bruhat.util import factorial, all_perms, all_ders, cross, choose
from bruhat.theta import divisors
from bruhat import series
from bruhat.series import Series, ring
from bruhat.argv import argv


# Do we bless the empty set with a partition,
# even though it doesn't deserve to have one?
EMPTY_SET_HAS_A_PARTITION = True
# Yes: these partitions are really surjections (up to isomorphism).

if argv.no_bless:
    print("EMPTY_SET_HAS_A_PARTITION = False")
    EMPTY_SET_HAS_A_PARTITION = False


def iterlen(I):
    return len(tuple(I))


def all_subsets(els):
    els = tuple(els)
    n = len(els)
    if n==0:
        yield ()
        return
    if n==1:
        yield ()
        yield els
        return
    bits = [(0,1)]*n
    for select in cross(bits):
        items = tuple(els[i] for i in range(n) if select[i])
        yield items


def all_cycles(els):
    els = tuple(els)
    n = len(els)
    if n<3:
        yield els
        return
    if n==3:
        yield els
        yield (els[0], els[2], els[1])
        return
    for perm in all_perms(els[1:]):
        yield (els[0],) + perm



# See: https://oeis.org/A121860
#    For example, inequivalent representatives of the a(4) = 8 matrices are:
#      [1 2 3 4]
#    .
#      [1 2] [1 2] [1 3] [1 3] [1 4] [1 4]
#      [3 4] [4 3] [2 4] [4 2] [2 3] [3 2]
#    .
#      [1]
#      [2]
#      [3]
#      [4]

def all_rects_shape(els, m, n):
    # A rect of shape (m,n) 
    # is a tuple of length m, of length n tuple's.
    assert m*n == len(els), (m, n, len(els))
    els = tuple(els)
    if m==1:
        row = els
        yield (row,)
        return
    if n==1:
        col = tuple((e,) for e in els)
        yield col
        return

    found = set()
    top = els[0]
    rest = els[1:]
    for row in choose(rest, n-1):
        rest1 = [e for e in rest if e not in row]
        for col in choose(rest1, m-1):
            rest2 = [e for e in rest1 if e not in col]
            for sub in all_perms(rest2):
                rect = [(top,)+row]
                k = 0
                for i in col:
                    row1 = (i,) + sub[k:k+n-1]
                    k += n-1
                    rect.append(row1)
                rect = tuple(rect)
                assert rect not in found
                found.add(rect)
                yield rect


def all_rects(els):
    n = len(els)
    if not n:
        return
    for d in divisors(n):
        #print(n, d, len(els))
        nrows, ncols = d, n//d
        for rect in all_rects_shape(els, nrows, ncols):
            yield rect

def transpose(rect):
    if not rect:
        return rect
    m = len(rect) # rows
    n = len(rect[0]) # cols
    t = tuple(tuple(rect[i][j] for i in range(m)) for j in range(n))
    return t

assert transpose(((0,1),(2,3))) == ((0,2),(1,3))



def all_binary_trees(els):
    if not els:
        return
    els = list(els)
    if len(els) == 1:
        yield els[0]
    #els = [set(e) for e in els]
    found = set()
    n = len(els)
    for i in range(n):
      for j in range(i+1, n):
        a, b = els[i], els[j]
        if str(b)<str(a):
            a, b = b, a
        pair = (a, b)
        rest = list(els)
        rest.pop(j)
        rest.pop(i)
        rest.insert(0, pair)
        for tree in all_binary_trees(rest):
            if tree not in found:
                yield tree
                found.add(tree)
        
    
#for tree in all_binary_trees([1,2,3,4]):
#    print(tree)

assert list(all_binary_trees([1,2,3])) == [((1, 2), 3), ((1, 3), 2), ((2, 3), 1)]


# -------------------------------------------------------


def all_obinary_trees(els):
    if not els:
        return
    els = list(els)
    if len(els) == 1:
        yield els[0]
    #els = [set(e) for e in els]
    found = set()
    n = len(els)
    for perm in all_perms(els):
        for i in range(n):
          for j in range(i+1, n):
            a, b = perm[i], perm[j]
            pair = (a, b)
            rest = list(perm)
            rest.pop(j)
            rest.pop(i)
            rest.insert(0, pair)
            for tree in all_obinary_trees(rest):
                if tree not in found:
                    yield tree
                    found.add(tree)
        
    
assert iterlen(all_obinary_trees([1,2,3])) == 12
assert iterlen(all_obinary_trees([1,2,3,4])) == 120


# -------------------------------------------------------


class Set(object): # copied from bruhat.rel
    def __init__(self, items=[]):
        if type(items) is int:
            assert items <= len(letters)
            items = letters[:items]
        #if items:
        #    tp = type(items[0])
        #    for item in items:
        #        assert type(item) is tp

        items = list(items)
        try:
            items.sort() # eeeeck: watch this !
        except TypeError:
            items.sort(key = str) # suck on this python3 
        except AssertionError:
            items.sort(key = str) # wup

        self.items = tuple(items)
        self.set_items = set(items)
        assert len(self.set_items) == len(self.items), self.items

    def __str__(self):
        return "{%s}" % (', '.join(str(x) for x in self))

    def __repr__(self):
        return "Set(%s)"%(str(list(self.items)))

    @classmethod
    def promote(cls, items):
        if isinstance(items, Set):
            return items
        #if type(items) is GeneratorType:
        return Set(items)

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

    def __contains__(self, item):
        return item in self.set_items

    def __eq__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items == b.items

    def __ne__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items != b.items

    def __lt__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items < b.items

    def __le__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items < b.items

    def __hash__(self):
        return hash(self.items)

    def __add__(a, b):
        # categorical coproduct
        items = [(0, ai) for ai in a] + [(1, bi) for bi in b]
        return Set(items)

    def union(a, b):
        items = a.set_items.union(b.set_items)
        return Set(items)

    def __mul__(a, b):
        # categorical product
        items = [(ai, bi) for ai in a for bi in b]
        return Set(items)
        
    def all_partitions(self):
        # A partition is a Set of non-empty subsets of self
        # such that: disjoint, non-empty (?), covering.
        items = self.items
        if not items:
            if EMPTY_SET_HAS_A_PARTITION:
                yield empty
            return
        #assert items
        if len(items) == 1:
            yield Set((self,))
            return
        if len(items) == 2:
            yield Set((self,))
            a, b = self
            q = Set([Set([a]), Set([b])])
            yield q
            return
        #head = Set((Set(items[:1]),))
        head = Set(items[:1])
        tail = Set(items[1:])
        #head = items[0]
        #tail = items[1:]
        #print("head:", head)
        #print("tail:", tail)
        for partition in tail.all_partitions():
            #print("partition:", partition)
            partition = partition.items # a tuple of Set's of Set's
            #print("%s + %s" % ((head,), partition))
            yield Set((head,) + partition)
            for i in range(len(partition)):
                left = partition[:i]
                right = partition[i+1:]
                middle = partition[i].union(head) # a Set
                yield Set(left + (middle,) + right)

    def all_parts2(self):
        "all the ways to break self into two subsets"
        items = self.items
        n = len(items)
    #    if n==0:
    #        yield items, items
    #        return
    #
    #    if n==1:
    #        yield [], items
    #        yield items, []
    #        return
    #
    #    if n==2:
    #        yield [], items
    #        yield [items[0]], [items[1]]
    #        yield [items[1]], [items[0]]
    #        yield items, []
    #        return
    
        bits = [(0, 1)]*n
        for idxs in cross(bits):
            left = Set([items[i] for i in range(n) if idxs[i]==0])
            right = Set([items[i] for i in range(n) if idxs[i]==1])
            yield left, right

empty = Set()
star = Set(["*"])

assert len(list(Set(0).all_parts2())) == 1
assert len(list(Set(1).all_parts2())) == 2
assert len(list(Set(3).all_parts2())) == 8
assert len(list(Set(3).all_partitions())) == 5


# -------------------------------------------------------

class Species(object):
    def __init__(self, f, name=None):
        self.construct = f
        self.name = name or self.__class__.__name__

    def __repr__(self):
        return self.name

    def __str__(self):
        return "%s%s"%(self.__class__.__name__, self.sequence())

    def __getitem__(self, U):
        return self.construct(U)

    def diff(F, degree=star):
        return DiffSpecies(F, degree)

    def point(F):
        return PointSpecies(F)

    def __add__(F, G):
        return AddSpecies(F, G)

    def mul(F, G):
        "Cartesian/Hadamard product"
        return MulSpecies(F, G)

    def __mul__(F, G):
        "Partitional/Cauchy product"
        return DotSpecies(F, G)
    dot = __mul__

    def __call__(F, G):
        return ComposeSpecies(F, G)

    def dirichlet(F, G):
        return DirichletSpecies(F, G)
    __matmul__ = dirichlet # why not?

    def sequence(F, n=6):
        items = []
        for i in range(n):
            U = Set(i)
            items.append(iterlen(F[U]))
        return items


class Number(Species):
    def __init__(self, n):
        self.n = n

    def __getitem__(self, U):
        if len(U):
            return empty
        return range(self.n)


class _UnaryOpSpecies(Species):
    def __init__(self, F):
        self.F = F


class DiffSpecies(_UnaryOpSpecies):
    def __init__(self, F, degree=star):
        self.degree = star # could be any set here
        self.F = F

    def __getitem__(self, U):
        U = Set.promote(U)
        U = U + self.degree # disjoint union
        return self.F[U]


class PointSpecies(_UnaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        for s in F[U]:
          for u in U: # the point
            yield (s, u)


class _BinaryOpSpecies(Species):
    def __init__(self, F, G):
        self.F = F
        self.G = G


class AddSpecies(_BinaryOpSpecies):
    def __getitem__(self, U):
        #return list(self.F[U]) + list(self.G[U])
        # Disjoint union
        for left in self.F[U]:
            yield (0, left)
        for right in self.G[U]:
            yield (1, right)
    

class MulSpecies(_BinaryOpSpecies):
    def __getitem__(self, U):
        return Set.promote(self.F[U]) * Set.promote(self.G[U])


class Data(object):

    pass

#    def __eq__(self, other):
#        return self.data == other.data


class Pair(Data):
    def __init__(self, a, b):
        self.data = (a, b)



class DotSpecies(_BinaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        G = self.G
        items = []
        for idx, (left, right) in enumerate(U.all_parts2()):
            left = Set.promote(F[left])
            right = Set.promote(G[right])
            for item in left*right:
                items.append((idx, item))
        return Set(items)


class ComposeSpecies(_BinaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        G = self.G
#        print("ComposeSpecies.__getitem__", F, G, U)
        items = []
        for p in U.all_partitions():
            # the partition p is a Set of Set's (subsets of U)
#            print("partition:", p)
            factors = []
            for V in p:
#                print("V = %s, len(V)=%s, G[V] = %s" % (V, len(V), G[V]))
                factors.append(tuple(G[V]))
            for struct in F[p]:
#                print("F[p]:", struct)
                #item = Set.mul_many([struct] + factors)
                for item in cross(factors):
                    yield (struct, item)
                #items.append(item)

        #return Set(items)

    
class DirichletSpecies(_BinaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        G = self.G

        for rect in all_rects(U):
            rows = Set(rect)
            cols = Set(transpose(rect))
            ts = list(G[cols])
            for s in F[rows]:
                for t in ts:
                    yield s, t


# ----------------------------------------------------

Zero = Species(lambda items : Set([]), "Zero")
One = Species(lambda items : (Set([items]) if len(items)==0 else empty), "One")
X = Species(lambda items : (Set([items]) if len(items)==1 else empty), "X")
E = Species(lambda items : Set([items]), "E")
E_plus = Species(lambda items : (Set([items]) if len(items)>0 else empty), "E_plus")
E_2 = Species(lambda items : (Set([items]) if len(items)==2 else empty), "E_2")
Point = Species(lambda items : Set.promote(items), "Point")
List = Species(lambda items : Set(all_perms(items)), "List")
Perm = List
Cycle = Species(lambda items : Set(all_cycles(items)))
Der = Species(lambda items : Set(all_ders(items)), "Der")
Par = Species(lambda items : Set(items.all_partitions()), "Par")
BinaryTree = Species(all_binary_trees, "BinaryTree")
OrderedBinaryTree = Species(all_obinary_trees, "OrderedBinaryTree")
Pow = Species(all_subsets, "Pow")


# ----------------------------------------------------
    
class GeneratingFunction(Series):
    def __init__(self, species):
        Series.__init__(self)
        self.species = species

class EGF(GeneratingFunction):
    def __getitem__(self, idx):
        U = Set(idx)
        FU = self.species[U]
        FU = list(FU)
        return ring.one*len(FU) / factorial(idx)

class OGF(GeneratingFunction):
    def __getitem__(self, idx):
        U = Set(idx)
        FU = self.species[U]
        FU = list(FU)
        return ring.one*len(FU)


# ----------------------------------------------------

def test_eq(F, G, n=6, skip_empty=False):
    start = 1 if skip_empty else 0
    for i in range(start, n):
        U = Set(i)
        if iterlen(F[U]) != iterlen(G[U]):
            return False
    return True



def test():

    for n in range(0):
        U = Set(n)
        print("U =", U)
        print("\tZero:", Zero[U])
        print("\tX:", X[U])
        print("\tE:", E[U])
        print("\tPoint:", Point[U])

    assert EGF(Zero).eq(series.zero)
    assert EGF(One).eq(series.one)
    assert EGF(X).eq(series.x)
    assert EGF(E).eq(series.exp)
    #print(EGF(E))

    f = EGF(List) # = 1/(1-x)
    for i in range(5):
        assert f[i] == 1

    assert EGF(X+X).eq(series.x + series.x)

    XX = X.dot(X)
    assert EGF(XX).eq(series.x * series.x)

    assert test_eq(X*E, Point)
    assert test_eq(List.diff(), List*List)
    
    assert test_eq(E.diff(), E)
    if EMPTY_SET_HAS_A_PARTITION:
        assert test_eq(Par.diff(), E*Par)
    else:
        assert test_eq(Par.diff(), E*Par + E)

    
    if EMPTY_SET_HAS_A_PARTITION:
        assert Par.sequence(6) == [1, 1, 2, 5, 15, 52]
        assert E(E_plus).sequence(6) == [1, 1, 2, 5, 15, 52] # whoops
        assert E(X).sequence(4) == [1, 1, 1, 1]
        assert E(E).sequence(4) == [1, 1, 2, 5]
        assert Perm.sequence(5) == (E.dot(Der)).sequence(5)

    else:
        # misses the first term... boo
        assert Par.sequence(6) == [0, 1, 2, 5, 15, 52]
        assert E(E_plus).sequence(6) == [0, 1, 2, 5, 15, 52]
        assert E(X).sequence(4) == [0, 1, 1, 1]
        assert E(E).sequence(4) == [0, 1, 2, 5]

    def r(n):
        result = 0
        for d in divisors(n):
            result += factorial(n) // (factorial(d) * factorial(n//d))
        return result

    for i in range(1, 9):
        assert iterlen(all_rects(list(range(i)))) == r(i)

    # https://oeis.org/A323295
    # ways to fill a matrix with the first n positive integers.
    LL = List.dirichlet(List)
    assert LL.sequence(6) == [0, 1, 4, 12, 72, 240]

    EE = E.dirichlet(E)
    assert EE.sequence(9) == [0, 1, 2, 2, 8, 2, 122, 2, 1682] # same as r(i) above

    # List == One + X.dot(List)
    R = One + X*List
    assert List.sequence(6) == R.sequence(6)

    # https://oeis.org/A001147
    # See also, Stanley, Vol 2, example 5.2.6
    assert BinaryTree.sequence(7) == [0, 1, 1, 3, 15, 105, 945]
    # From Stanley, Vol 2, page 178
    # EGF(x) == 1 - sqrt(1-2*x)

    assert test_eq(Cycle.point(), Perm, skip_empty=True)

    assert test_eq(Pow, E*E)
    assert Pow.sequence(6) == [1, 2, 4, 8, 16, 32]
    
    items = [X, E, List, Der, Cycle, BinaryTree]
    for F in items:
        assert test_eq(F.point(), Point.mul(F))
        assert test_eq(F.point(), X*F.diff())

        assert test_eq(F*One, One*F)
        assert test_eq(F*One, F)

        assert test_eq(F*Zero, Zero*F)
        assert test_eq(F*Zero, Zero)

        assert test_eq(F+F, Number(2)*F)
        assert test_eq(F+F+F, Number(3)*F)

        assert test_eq(F(X), X(F), skip_empty=True)
        assert test_eq(F(X), F, skip_empty=True)  # identity

        assert test_eq(F@X, F, skip_empty=True)
        assert test_eq(F@X, F(X), skip_empty=True)

    for F in items:
      for G in items:

        # From "Introduction to the Theory of Species of Structures" 
        # Francois Bergeron, Gilbert Labelle, and Pierre Leroux 
        # http://bergeron.math.uqam.ca/wp-content/uploads/2013/11/book.pdf

        # Exercise 2.38
        assert test_eq( (F+G).point(), F.point()+G.point() )
        assert test_eq( (F*G).point(), F.point()*G+F*G.point() )

        # Exercise 2.36
        if EMPTY_SET_HAS_A_PARTITION:
            assert test_eq( (F(G)).point(), F.diff()(G)*G.point() )
            assert test_eq( 
                E(G) * (F.diff() + G.diff()*F),
                (E(G) * F).diff() )

        # Exercise 2.44
        assert test_eq(F.mul(G), G.mul(F)) # comm
        for H in items:
            assert test_eq(F.mul(G+H), F.mul(G)+F.mul(H)) # comm
        # etc. etc.

    for F in items:
        assert test_eq(F.mul(E), F) # identity
        assert test_eq(E.mul(F), F) # identity


    for F in items:
      for G in items:

        assert test_eq( (F+G).diff()  , F.diff() + G.diff() )
        assert test_eq( (F*G).diff()  , F.diff()*G + F*G.diff() )

        if EMPTY_SET_HAS_A_PARTITION:
            assert test_eq( (F(G)).diff() , F.diff()(G) * G.diff())
        #print((F(G)).diff().sequence(), 
        #    (F.diff()(G) * G.diff()).sequence())


        assert test_eq(F*G, G*F) # comm
        assert test_eq(F@G, G@F) # comm
        assert test_eq((F@G).point(), F.point()@G.point())

        for H in items:

            assert test_eq( F(G(H)), (F(G))(H) ) # assoc
            assert test_eq( (F+G)(H),  F(H) + G(H) ) # dist
            if EMPTY_SET_HAS_A_PARTITION:
                assert test_eq( (F*G)(H),  F(H) * G(H) ) # dist

            assert test_eq(F*(G*H), (F*G)*H) # assoc
            assert test_eq(F*(G+H), F*G + F*H) # dist

            assert test_eq(F@(G@H), (F@G)@H) # assoc
            assert test_eq(F@(G+H), F@G + F@H) # dist

            #assert test_eq(F@(G*H), (F@G) * (F@H)) # nope!
            #assert test_eq((F*G)@H, (F@H) * (G@H)) # nope!

    print("OK")


    
def main():

    print(BinaryTree.sequence(7)) # 0, 1, 1, 3, 15, 105, 945

    F = X + E_2(BinaryTree)
    print(F.sequence(7)) # 0, 1, 1, 3, 15, 105, 945

    # https://oeis.org/A001813
    # Quadruple factorial numbers: a(n) = (2n)!/n!. 
    #print(OrderedBinaryTree.sequence(6)) # 0, 1, 2, 12, 120, 1680

    P = OrderedBinaryTree
    #print(P.sequence(6))
    #Q = X + X*Cycle(P)
    Q = X*List(P)
    print(Q.sequence(6))

    #F = OrderedBinaryTree * OrderedBinaryTree
    #print(F.sequence(6)) # 0, 0, 2, 12, 120, 1680

    return

    # Conjecture: these are "wavefronts" on a BinaryTree.
    F = BinaryTree(BinaryTree)
    #print(F.sequence(7)) # == [0, 1, 2, 9, 63, 600, 7245]

    F = BinaryTree(E)

    # Number of planted binary phylogenetic trees with n labels. 
    # http://oeis.org/A006677
    assert F.sequence(6) == [0, 1, 2, 7, 41, 346]

    FF = F.dirichlet(F)

    print(F.sequence(7))
    print(FF.sequence(7))

if __name__ == "__main__":

    if argv.test:
        test()
    else:
        main()












