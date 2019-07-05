#!/usr/bin/env python3

"""
Combinatorial species.

See chapter 4 here:
    COMBINATORIAL SPECIES AND LABELLED STRUCTURES 
    Brent Abraham Yorgey 
    https://www.cis.upenn.edu/~sweirich/papers/yorgey-thesis.pdf
"""

import string
letters = list(string.ascii_lowercase)

from bruhat.util import factorial, allperms, allders, cross, choose
from bruhat.theta import divisors
from bruhat import series
from bruhat.series import Series, ring



def iterlen(I):
    return len(tuple(I))


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
            for sub in allperms(rest2):
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
        items.sort() # eeeeck: watch this !
        self.items = tuple(items)
        self.set_items = set(items)
        assert len(self.set_items) == len(self.items)

    def __str__(self):
        return "{%s}" % (', '.join(str(x) for x in self))

    def __repr__(self):
        return "Set(%s)"%(str(list(self.items)))

    @classmethod
    def promote(cls, items):
        if isinstance(items, Set):
            return items
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

assert len(list(Set(3).all_parts2())) == 8
assert len(list(Set(3).all_partitions())) == 5


# -------------------------------------------------------

class Species(object):
    def __init__(self, f, name=None):
        self.construct = f
        self.name = name

    def __str__(self):
        return self.name or self.__class__.__name__+"(?)"

    def __getitem__(self, U):
        return self.construct(U)

    def __add__(F, G):
        return AddSpecies(F, G)

    def __mul__(F, G):
        "Cartesian/Hadamard product"
        return MulSpecies(F, G)

    def dot(F, G):
        "Partitional/Cauchy product"
        return DotSpecies(F, G)

    def __call__(F, G):
        return ComposeSpecies(F, G)

    def dirichlet(F, G):
        return DirichletSpecies(F, G)

    def sequence(F, n):
        items = []
        for i in range(n):
            U = Set(i)
            items.append(iterlen(F[U]))
        return items


class BinaryOpSpecies(Species):
    def __init__(self, F, G):
        self.F = F
        self.G = G


class AddSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        return self.F[U] + self.G[U]
    

class MulSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        return self.F[U] * self.G[U]


class DotSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        G = self.G
        items = []
        for idx, (left, right) in enumerate(U.all_parts2()):
            for item in F[left] * G[right]:
                items.append((idx, item))
        return Set(items)


class ComposeSpecies(BinaryOpSpecies):
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

    
class DirichletSpecies(BinaryOpSpecies):
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
Point = Species(lambda items : Set.promote(items), "Point")
List = Species(lambda items : Set(allperms(items)), "List")
Perm = List
Der = Species(lambda items : Set(allders(items)), "Der")
Par = Species(lambda items : Set(items.all_partitions()), "Par")
BinaryTree = Species(all_binary_trees, "BinaryTree")


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

    assert X.dot(E).sequence(6) == Point.sequence(6)

    assert Par.sequence(6) == [0, 1, 2, 5, 15, 52]
    assert E(E_plus).sequence(6) == [0, 1, 2, 5, 15, 52]

    assert E(X).sequence(4) == [0, 1, 1, 1] # misses the first term... boo

    assert Perm.sequence(5) == (E.dot(Der)).sequence(5)

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
    R = One + X.dot(List)
    assert List.sequence(6) == R.sequence(6)

    # https://oeis.org/A001147
    # See also, Stanley, Vol 2, example 5.2.6
    assert BinaryTree.sequence(7) == [0, 1, 1, 3, 15, 105, 945]
    # From Stanley, Vol 2, page 178
    # EGF(x) == 1 - sqrt(1-2*x)


    F = BinaryTree(BinaryTree)
    #print(F.sequence(7)) # == [0, 1, 2, 9, 63, 600, 7245]
    # Conjecture: these are "wavefronts" on a BinaryTree.

    


if __name__ == "__main__":

    test()










