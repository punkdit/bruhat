#!/usr/bin/env python

"""
previous version: higher.py
see also: _theory_z3.py

"""

from time import time
start_time = time()


class Sort(object):
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return self.name
    def __repr__(self):
        return "Sort(%r)"%(self.name,)


# Variable's are universally quantified,
# Const's are existentially quantified.

class Expr(object):
    DEBUG = False
    depth = 0
    def info(self, *msg):
        if not self.DEBUG:
            return
        indent = " "*(Expr.depth*2 - 1)
        if indent:
            print(indent, *msg)
        else:
            print(*msg)
    def push(self):
        if self.DEBUG:
            Expr.depth += 1
    def pop(self):
        if self.DEBUG:
            Expr.depth -= 1

    def __init__(self, theory, key, sort, is_const):
        self.theory = theory
        self.key = key # hashable
        self.sort = sort
        self.is_const = is_const
        self.parent = None

    def find_root(self):
        expr = self
        count = 0
        while expr.parent is not None:
            assert expr.parent is not expr
            expr = expr.parent
            count += 1
            assert count < 100
        return expr

    def find_root_rebuild(self):
        expr = self.find_root()
        expr = expr.rebuild()
        return expr

    def rewrite(lhs, rhs):
        assert isinstance(rhs, Expr), rhs
        lhs.info("Expr.rewrite", lhs, "-->", rhs)
        lhs.push()
        #lhs = lhs.find_root()
        #rhs = rhs.find_root_rebuild()
        lhs = lhs.rebuild()
        rhs = rhs.rebuild()
        lhs = lhs.find_root()
        rhs = rhs.find_root()
        if lhs is not rhs:
            lhs.info("Expr.rewrite %s --> %s" % (lhs, rhs))
            #if str(rhs) == "((h*g)*Y.identity)":
            #    assert 0
            lhs.parent = rhs
        lhs.pop()

    def __hash__(self):
        assert 0, "use .key"

    def __eq__(self, other):
        assert isinstance(other, Expr), other
        lhs = self.find_root_rebuild()
        rhs = other.find_root_rebuild()
        #lhs = self.find_root()
        #rhs = other.find_root()
        return lhs is rhs

    def __mul__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("*", self.sort, other.sort)
        if op is None:
            raise AttributeError("operator %s*%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __lshift__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("<<", self.sort, other.sort)
        if op is None:
            raise AttributeError("operator %s<<%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __matmul__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("@", self.sort, other.sort)
        if op is None:
            raise AttributeError("operator %s@%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __getattr__(self, attr):
        theory = self.theory
        op = theory.get_operator(attr, self.sort)
        if op is None:
            raise AttributeError("attr %r not found"%(attr,))
        return op(self)

    @staticmethod
    def match(src, tgt, send=None, depth=0):
        if send is None:
            send = {}
        indent = "  "*depth
        #print(indent+"match:", src, "->", tgt, send)
        if src.sort != tgt.sort:
            return None
        if send:
            src = src.substitute(send) # ?
        if isinstance(src, Variable):
            if src.key in send and send[src.key] is not tgt:
                return None # <------------ return
            send[src.key] = tgt
            #print(indent+"match: Const")
            return send # <-------------- return
        if isinstance(tgt, Variable):
            return None # <-------------- return
        if src is tgt:
            return send
        if isinstance(src, Const) or isinstance(tgt, Const):
            return None
        assert isinstance(src, Term), (src, tgt)
        assert isinstance(tgt, Term), (src, tgt)
        if src.op != tgt.op:
            return None # fail
        for left, right in zip(src.args, tgt.args):
            send = Expr.match(left, right, send, depth=depth+1) # recurse
            if send is None:
                return None # fail
        return send


class Variable(Expr):
    def __init__(self, theory, name, sort):
        self.name = name
        self.sort = sort
        key = id(self)
        Expr.__init__(self, theory, key, sort, False)

    def __str__(self):
        return self.name
        #return "%s:%s"%(self.name, self.sort)
    __repr__ = __str__

    def rebuild(self):
        assert 0, "i'm a Variable!"

    def find(self, v):
        assert isinstance(v, Const)
        return self is v

    def substitute(self, send):
        return send.get(self.key, self)

    def all_vars(self):
        return {self.key : self}


class Const(Expr):
    def __init__(self, theory, name, sort):
        self.name = name
        self.sort = sort
        key = id(self)
        Expr.__init__(self, theory, key, sort, True)

    def __str__(self):
        return self.name
        #return "%s:%s"%(self.name, self.sort)
    __repr__ = __str__

    def rebuild(self):
        #assert self.parent is None
        return self

    def find(self, v):
        assert isinstance(v, Const)
        return self is v

    def substitute(self, send):
        return send.get(self.key, self)

    def all_vars(self):
        return {self.key : self}


class Operator(object):
    def __init__(self, theory, name, rsort, argsorts, inline=False, postfix=False):
        self.theory = theory
        self.name = name
        self.rsort = rsort
        self.argsorts = argsorts
        self.inline = inline
        self.postfix = postfix

    def __str__(self):
        return self.name

    def __repr__(self):
        return "%s(%r)"%(self.__class__.__name__, self.name)

    def __call__(self, *args):
        theory = self.theory
        argsorts = self.argsorts
        assert len(args) == len(argsorts)
        for (arg, sort) in zip(args, argsorts):
            assert arg.sort == sort
        term = theory.Term(self, *args, 
            sort=self.rsort, inline=self.inline, postfix=self.postfix)
        return term


class Term(Expr):
    def __init__(self, theory, op, *args, sort=None, inline=False, postfix=False):
        assert isinstance(op, Operator)
        self.op = op
        self.args = args
        self.inline = inline
        self.postfix = postfix
        key = (op,) + tuple(arg.key for arg in args)
        is_const = True
        for arg in args:
            if not arg.is_const:
                is_const = False
        Expr.__init__(self, theory, key, sort, is_const)

    def __str__(self):
        args = self.args
        op = self.op
        if self.inline:
            assert len(args) == 2
            l, r = args
            s = "(%s%s%s)" % (l, op, r)
        elif self.postfix:
            assert len(args) == 1
            l, = args
            s = "%s.%s"%(l, op)
        elif len(args) == 0:
            s = str(op) # Nullary
        else:
            s = "%s%s"%(op, args)
        return s
    __repr__ = __str__

    count = 0
    def rebuild(self):
        #assert self.parent is None
        self.info("Expr.rebuild", self, id(self))
        self.push()
        if self.DEBUG:
            assert self.depth < 10
            #assert self.count < 10
            #self.count += 1
        expr = self
        op = self.op
        for arg in self.args:
            assert arg is not self
            self.info("Expr.rebuild\t", arg, "-->", arg.parent, id(arg.parent))
        #args = [arg.find_root_rebuild() for arg in self.args] # recurse
        args = [arg.rebuild() for arg in self.args] # recurse
        term = self.theory.Term(op, *args, 
            sort=self.sort, inline=self.inline, postfix=self.postfix)
        if term is not expr:
            #term = term.find_root_rebuild() # recurse !
            term = term.rebuild() # recurse !
            assert term is not expr
            assert term.parent is None
            self.info("Expr.rebuild =", term)
            assert expr.parent is None
            expr.parent = term # equate these
        self.pop()
        #assert term.parent is None
        return term

    def __getitem__(self, idx):
        return self.args[idx]

    def find(self, v):
        assert isinstance(v, Const)
        for arg in self.args:
            if arg.find(v):
                return True
        return False

    def substitute(self, send):
        op = self.op
        args = [arg.substitute(send) for arg in self.args]
        theory = self.theory
        term = theory.Term(op, *args, 
            sort=self.sort, inline=self.inline, postfix=self.postfix)
        return term

    def all_vars(self):
        vs = {}
        for arg in self.args:
            for v in arg.all_vars().values():
                vs[v.key] = v
        return vs


class Rewrite(object):
    def __init__(self, theory, src, tgt):
        assert isinstance(src, Expr)
        assert isinstance(tgt, Expr)
        assert theory is src.theory
        assert theory is tgt.theory
        assert src.sort == tgt.sort
        assert not src.is_const, 'makes no sense..'
        self.src = src
        self.tgt = tgt

    def __str__(self):
        return "%s -> %s" % (self.src, self.tgt)

    def match(self, expr, directional=True):
        assert isinstance(expr, Expr)
        src, tgt = self.src, self.tgt
        #print("Rewrite.match", src, "-->", expr)
        send = Expr.match(src, expr)
        #print("Rewrite.match =", send)
        if send is not None:
            tgt = tgt.substitute(send)
            return tgt


class Theory(object):
    DEBUG = False
    def info(self, *msg):
        if self.DEBUG:
            print(*msg)

    def __init__(self):
        self.cache = {} # cache all constant expressions
        self.sigs = {}
        self.eqns = []
        self.lookup = {}

    def dump(self):
        vs = list(self.cache.values())
        print("Theory.dump")
        for v in vs:
            print('\t', v, end=" ")
            v = v.parent
            while v:
                print("-->", v, end=" ")
                v = v.parent
            print()
        print()

    def Const(self, name, sort, expr=None):
        if expr is None:
            expr = Const(self, name, sort)
        #solver, lookup = self.solver, self.lookup
        #v = solver.get_const(name)
        #lookup[expr.key] = v
        self.info("Const:", expr)
        return expr

    def Variable(self, name, sort):
        expr = Variable(self, name, sort)
        self.info("Variable:", expr)
        return expr

    def rewrite(self, lhs, rhs):
        lhs.rewrite(rhs)

    def apply_rewrite(self, eqn, expr):
        while 1:
            other = eqn.match(expr)
            if other is None:
                return
            assert other.key != expr.key, "continue?"
            self.info("Theory.apply_rewrite:", expr, "-->", other)
            self.rewrite(expr, other)
            expr = other

    def Term(self, op, *args, sort=None, inline=False, postfix=False):
        cache = self.cache
        expr = Term(self, op, *args, sort=sort, inline=inline, postfix=postfix)
        if not expr.is_const:
            return expr # <------------- return

        if expr.key in cache:
            return cache[expr.key] # <------------- return

        # we made a new Term
        cache[expr.key] = expr

        rewrite = False
        args = list(args)
        for i in range(len(args)):
            arg = args[i].find_root()
            if arg is not args[i]:
                rewrite = True
                args[i] = arg
        if rewrite:
            term = self.Term(op, *args, sort=sort, inline=inline, postfix=postfix)
            self.rewrite(expr, term)

        self.info("Theory.Term:", expr)
        for eqn in self.eqns:
            self.apply_rewrite(eqn, expr) # may recurse !

        return expr

    def add_operator(self, op):
        sig = (op.name, tuple(op.argsorts))
        assert sig not in self.sigs
        self.sigs[sig] = op

    def Operator(self, name, rsort, argsorts, inline=False, postfix=False):
        op = Operator(self, name, rsort, argsorts, inline=inline, postfix=postfix)
        if argsorts:
            self.add_operator(op)
        return op

    def get_operator(self, name, *argsorts):
        sig = (name, argsorts)
        return self.sigs.get(sig, None)

    def Nullary(self, name, sort):
        op = self.Operator(name, sort, ())
        expr = op()
        return expr

    def Rewrite(self, lhs, rhs):
        if lhs.is_const:
            assert rhs.is_const
            self.rewrite(lhs, rhs) # now these constant's are equal
            return
        eqn = Rewrite(self, lhs, rhs)
        self.info("theory.Rewrite:", eqn)
        self.eqns.append(eqn)

        for expr in list(self.cache.values()):
            assert expr.is_const, expr
            self.apply_rewrite(eqn, expr)
        return eqn

    def Equation(self, lhs, rhs):
        self.Rewrite(lhs, rhs)
        self.Rewrite(rhs, lhs)

    #def equate(self, lhs, rhs):
    #    solver, lookup = self.solver, self.lookup
    #    lhs, rhs = lookup[lhs.key], lookup[rhs.key]
    #    solver.equate(lhs, rhs)

    #def is_equal(self, lhs, rhs):
    #    solver, lookup = self.solver, self.lookup
    #    lhs, rhs = lookup[lhs.key], lookup[rhs.key]
    #    return solver.is_equal(lhs, rhs)



def distinct(items):
    n = len(items)
    for i in range(n):
      for j in range(i+1, n):
        if items[i] == items[j]:
            return False
    return True
        

def test_base_theory():

    theory = Theory()
    sort = Sort("sort")

    a, b, c, d, e, f, g = [theory.Const(name, sort) for name in 'abcdefg']
    t, u, v, w, x, y, z = [theory.Variable(name, sort) for name in 'tuvwxyz']
    #one = theory.Nullary("1", sort)
    one = theory.Const("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Const("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    for (src, tgt) in [
        (u, u*v),
        (u*v, v*u),
        (u*v, v*(w*w)),
        ((u*v)*w, (x*y)*(z*t)),
    ]:
        send = Expr.match( src, tgt )
        assert send
        assert src.substitute(send).key == tgt.key

    send = Expr.match( u*v, v )
    assert send is None

    send = Expr.match( (u*v)*w, u*v )
    assert send is None

    eqn = Rewrite(theory, (u*v)*w, u*(v*w))

    lhs = (a*b)*(c*d)
    rhs = eqn.match(lhs)
    assert rhs is a*(b*(c*d)), rhs

    lhs = one*u
    rhs = one*(one*v)
    assert Expr.match(lhs, rhs)

    eqn = Rewrite(theory, one*u, u)
    lhs = one*(one*a)
    rhs = eqn.match(lhs)
    assert rhs is one*a

    _sort = Sort("_sort")
    _a = theory.Const("_a", _sort)
    assert Expr.match(u, _a) is None


def test_theory():

    theory = Theory()
    sort = Sort("monoid")

    a, b = [theory.Const(name, sort) for name in 'ab']
    x, y = [theory.Variable(name, sort) for name in 'xy']
    #one = theory.Nullary("one", sort)
    one = theory.Const("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a != b
    assert a*b != b*a

    theory.Rewrite( one*x, x )
    theory.Rewrite( x*one, x )

    assert a*one == a
    assert one*b == b
    assert distinct([a, b, a*b, b*a])
    #Theory.DEBUG = Expr.DEBUG = True
    #theory.dump()
    #lhs = a*one*one*a
    #theory.dump()

    assert a*one*one*a == a*a
    assert a*b != b*a
    assert one*one == one
    assert ((one*one)*a) == (one*a)

    for u in [one, a, b, a*b, b*a]:
      for v in [one, a, b, a*b, b*a]:
        assert u*one*v == u*v, "%s != %s" % (u*one*v, u*v)


def test_monoid_theory():

    theory = Theory()
    sort = Sort("monoid")

    a, b, c, d = [theory.Const(name, sort) for name in 'abcd']
    u, v, w = [theory.Variable(name, sort) for name in 'uvw']
    one = theory.Const("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Const("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    assert (a*b)*c != a*(b*c)

    #Theory.DEBUG = Expr.DEBUG = True
    theory.Equation( (u*v)*w, u*(v*w) )
    #theory.dump()
    assert (a*b)*c == a*(b*c)

    lhs = ((a*b)*c)*d
    assert lhs == (a*b)*(c*d)
    assert lhs == a*(b*(c*d))
    assert lhs == a*((b*c)*d)
    assert lhs == (a*(b*c))*d

    theory.Rewrite(one*u, u)
    theory.Rewrite(u*one, u)

    assert a*one == a

    assert a*one*one*a == a*a
    assert a*b != b*a

    assert distinct([a, b, c, d])


def test_category_theory():

    theory = Theory()

    cell0 = Sort("cell0") # object's
    cell1 = Sort("cell1") # morphism's

    Operator = theory.Operator
    Const = theory.Const
    Variable = theory.Variable
   
    Operator("identity", cell1, [cell0], postfix=True)
    Operator("*", cell1, [cell1, cell1], inline=True)

    X, Y, Z, W = [Const(name, cell0) for name in "XYZW"]
    f, g, h = [Const(name, cell1) for name in "fgh"]

    vX, vY, vZ, vW = [Variable(name, cell0) for name in "XYZW"]
    vf, vg, vh = [Variable(name, cell1) for name in "fgh"]

    assert X.identity is X.identity
    assert X.identity == X.identity
    assert X.identity != Y.identity


    #      h      g      f    
    #  W <--- Z <--- Y <--- X
    # 

    assert g*f != h*g
    theory.Equation( (vh*vg)*vf, vh*(vg*vf) )
    assert g*f != h*g

    for cell in [X, Y, Z, W]:
        theory.Rewrite( vf*cell.identity, vf )
        theory.Rewrite( cell.identity*vf, vf )

    assert g*f != h*g
    assert g*f != h*g
    assert distinct([X, Y, Z, W])
    assert distinct([X.identity, Y.identity, Z.identity, W.identity])

    assert g*f != h*g

    assert (h*g)*Y.identity == h*g
    assert g*f != h*g

    assert X.identity * X.identity == X.identity

    #theory.dump()

    assert (h*g)*f == h*(g*f)

    #Theory.DEBUG = Expr.DEBUG = True

    # add an isomorphism ... ?
    i = Const("i", cell1)
    j = Const("j", cell1)
    theory.Rewrite(i*j, X.identity)
    theory.Rewrite(j*i, X.identity)

    # now we do the following chain:
    #  (f*i)*j ----> f*(i*j) ----> f*X.identity -------> f

    assert f*i*j == f

    # look at this! all kinds of chains...
    #theory.dump()



if __name__ == "__main__":

    print("\n\n")
    test_base_theory()
    test_theory()
    test_monoid_theory()
    test_category_theory()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







