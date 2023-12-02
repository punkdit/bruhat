#!/usr/bin/env python3

"""
Here we make some Categories where the _objects are
natural numbers: 0, 1, 2, ...

"""

from bruhat.expr import Solver, Variable, Expr, Term
from bruhat.argv import argv


class Cat(Solver):
    def variable(self, name, tgt, src):
        return Gen(self, name, tgt, src)

class Morphism(object): # Mixin
    def __init__(self, tgt, src):
        self.tgt = int(tgt)
        self.src = int(src)

    def __mul__(self, other):
        cat = self.solver
        assert cat is other.solver
        return Mul(cat, self, other)


class Gen(Morphism, Variable):
    def __init__(self, cat, name, tgt, src):
        Variable.__init__(self, cat, name)
        Morphism.__init__(self, tgt, src)

    def __str__(self):
        return "%s(%d <-- %d)"%(self.name, self.tgt, self.src)


class Mul(Morphism, Term):
    def __init__(self, cat, left, right):
        op = cat.get_op("*", 2)
        Term.__init__(self, cat, op, (left, right), inline=True)
        assert isinstance(left, Morphism)
        assert isinstance(right, Morphism)
        assert left.src == right.tgt
        Morphism.__init__(self, left.tgt, right.src)
        if isinstance(left, Mul):
            a, b = left.args
            c = right
            #print("Mul", self)
            self.equate(a*(b*c))
        #elif isinstance(right, Mul):
        #    a = left
        #    b, c = right.args
        #    self.equate((a*b)*c)




def test_assoc():
    terms = "((a*b)*c)*d (a*(b*c))*d a*((b*c)*d) a*(b*(c*d)) (a*b)*(c*d)".split()
    for lhs in terms:
      for rhs in terms:
        if lhs == rhs:
            continue
        cat = Cat()
        variable = cat.variable
        a, b, c, d = [variable(s, 1, 1) for s in 'abcd']
        l = eval(lhs, locals())
        r = eval(rhs, locals())
        assert l == r

    cat = Cat()
    variable = cat.variable
    a, b, c, d, e = [variable(s, 1, 1) for s in 'abcde']

    lhs = (((a*b)*c)*d)*e
    rhs = a*(b*(c*(d*e)))
    assert lhs == rhs


def test_inj():
    cat = Cat()

    N = 5
    rows = []
    lookup = {}
    for i in range(N):
      row = []
      for j in range(i):
        #name = "inj(%d,%d)"%(i, j)
        #name = "inj_%d_%d"%(i, j)
        name = "i%d%d"%(i, j)
        op = cat.variable(name, i, i-1)
        #exec("%s = op"%(name), globals(), locals())
        lookup[i,j] = op
        row.append(op)
      rows.append(row)

    for i in range(N):
      for j in range(i):
        op = rows[i][j]
        print(op, end=" ")
      print()

    pop = lambda i,j: lookup[i,j]

    for tgt in range(1, N):
      for i in range(tgt):
       for j in range(tgt-1):
        e = lookup[tgt, i] * lookup[tgt-1, j]

    for src in range(2, N):
     for i in range(src):
      for j in range(src-1):
        lhs = pop(src, i) * pop(src-1, j)
        if i<j:
            rhs = pop(src, j-1) * pop(src-1, i)
        elif i>j:
            rhs = pop(src, j) * pop(src-1, i-1)
        else:
            assert i==j
            rhs = pop(src, j+1) * pop(src-1, i)
        #print(lhs, "==", rhs)
        lhs.equate(rhs)

    assert pop(3,0)*pop(2,0) == pop(3,1)*pop(2,0)
    assert pop(4,3)*pop(3,0) == pop(4,0)*pop(3,2)
    assert pop(4,1)*pop(3,0)*pop(2,0) == pop(4,1)*pop(3,1)*pop(2,0)

    assert pop(4,1)*pop(3,0) != pop(4,3)*pop(3,1)
    #assert pop(4,1)*pop(3,0)*pop(2,0) == pop(4,3)*pop(3,1)*pop(2,0)
    assert pop(4,1)*pop(3,0) != pop(4,3)*pop(3,1) # FAIL
    assert pop(4,1)*pop(3,0)*pop(2,0) != pop(4,3)*pop(3,1)*pop(2,0) # FAIL

        

def test():
    #test_assoc()
    test_inj()

    
    


if __name__ == "__main__":

    from time import sleep, time
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

    else:
        fn = eval(name)
        fn()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))


