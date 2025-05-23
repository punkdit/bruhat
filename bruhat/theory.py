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

class Debug(object):
    DEBUG = False
    depth = 0
    def info(self, *msg):
        if not self.DEBUG:
            return
        indent = " "*(self.__class__.depth*2 - 1)
        if indent:
            print(indent, *msg)
        else:
            print(*msg)
    def push(self):
        if self.DEBUG:
            self.__class__.depth += 1
    def pop(self):
        if self.DEBUG:
            self.__class__.depth -= 1


class Expr(Debug):

    def __init__(self, theory, key, sort, is_const, size):
        self.theory = theory
        self.key = key # hashable
        self.sort = sort
        self.is_const = is_const
        self.size = size
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
            if lhs.size < rhs.size: # argh, what am i doing?
                lhs, rhs = rhs, lhs
            assert lhs.size >= rhs.size
            lhs.parent = rhs
        lhs.pop()

    def __hash__(self):
        assert 0, "use .key"

    def __eq__(self, other):
        assert isinstance(other, Expr), other
        #assert self.sort is other.sort # too strict...
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
            raise AttributeError("binary operator %s*%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __lshift__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("<<", self.sort, other.sort)
        if op is None:
            raise AttributeError("binary operator %s<<%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __matmul__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("@", self.sort, other.sort)
        if op is None:
            raise AttributeError("binary operator %s@%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __getattr__(self, attr):
        theory = self.theory
        op = theory.get_operator(attr, self.sort)
        if op is None:
            raise AttributeError("unary operator %s.%s not found"%(self.sort, attr,))
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
        Expr.__init__(self, theory, key, sort, False, 1)

    def __str__(self):
        return self.name
        #return "%s:%s"%(self.name, self.sort)
    __repr__ = __str__

    def rebuild(self):
        assert 0, "%s is a Variable!"%(self,)

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
        Expr.__init__(self, theory, key, sort, True, 1)

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
        size = 0
        for arg in args:
            if not arg.is_const:
                is_const = False
            size += arg.size
        Expr.__init__(self, theory, key, sort, is_const, size)

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
        #if self.DEBUG:
        #    assert self.depth < 10
        #    #assert self.count < 10
        #    #self.count += 1
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
            #expr.theory.info("%s substitute %s"%(tgt, send))
            tgt = tgt.substitute(send)
            return tgt


class Theory(Debug):

    def __init__(self):
        self.cache = {} # cache all constant expressions
        self.sigs = {}
        self.eqns = []
        self.lookup = {}

    def dump(self, expr=None):
        # print the chains of equivelant Expr's
        assert expr is None or isinstance(expr, Expr)
        vs = list(self.cache.values())
        keys = dict((v.key, v) for v in vs)
        for v in vs:
            if v.parent is None:
                continue
            key = v.parent.key
            if key in keys:
                del keys[key]
        chains = []
        for v in keys.values():
            chain = [v]
            v = v.parent
            while v:
                chain.append(v)
                v = v.parent
            chains.append(chain)

        print("Theory.dump:", expr if expr else "")
        for chain in chains:
            if expr is None or expr in chain:
                print("\t", " --> ".join(str(v) for v in chain))
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

        self.push()
        self.info("Theory.Term:", expr)

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

        for eqn in self.eqns:
            self.apply_rewrite(eqn, expr) # may recurse !

        self.pop()

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


def test_group_theory():

    theory = Theory()
    sort = Sort("group")

    a, b, c, d = [theory.Const(name, sort) for name in 'abcd']
    u, v, w = [theory.Variable(name, sort) for name in 'uvw']
    one = theory.Const("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)
    theory.Operator("inv", sort, [sort], postfix=True)

    assert a is not b
    assert a is not theory.Const("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    assert (a*b)*c != a*(b*c)

    theory.Equation( (u*v)*w, u*(v*w) )
    theory.Rewrite( u*u.inv, one )
    theory.Rewrite( u.inv*u, one )
    theory.Rewrite( u.inv.inv, u )
    theory.Rewrite( (u*v).inv, v.inv * u.inv )
    theory.Rewrite( one.inv, one )

    #theory.dump()
    assert (a*b)*c == a*(b*c)

    assert a.inv.inv == a
    assert one.inv == one

    #Theory.DEBUG = Expr.DEBUG = True
    #theory.dump()

    ab = a*b
    assert ab.inv * ab == one
    assert ab.inv == b.inv * a.inv

    # Coxeter group ?
    theory.Rewrite( a*a, one )
    theory.Rewrite( b*b, one )
    #theory.Rewrite( a*b*a*b*a*b, one )
    #theory.Equation( a*b*a, b*a*b )
    theory.Rewrite( a*b*a, b*a*b )

    assert distinct([one, a, b, a*b, b*a, a*b*a])
    #Theory.DEBUG = Expr.DEBUG = True
    #a*b*a*b == b*a
    #theory.dump()
    #assert a*b*a*b == b*a # FAIL
    #assert a*b*a*b*a*b == one # FAIL


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

    # it actually turns out to be this one:
    #  f*X.identity -------> f ---> ((f*i)*j)
    #theory.dump()

    assert f*(i*j) == f


def build_category_theory(cell0=None, cell1=None):

    theory = Theory()

    if cell0 is None:
        cell0 = Sort("cell0") # object's
        cell1 = Sort("cell1") # morphism's
        theory.cell0 = cell0
        theory.cell1 = cell1

    Operator = theory.Operator
    Const = theory.Const
    Variable = theory.Variable
    Rewrite = theory.Rewrite
    Equation = theory.Equation
   
    Operator("identity", cell1, [cell0], postfix=True)
    Operator("*", cell1, [cell1, cell1], inline=True)
    Operator("inv", cell1, [cell1], postfix=True)
    Operator("src", cell0, [cell1], postfix=True)
    Operator("tgt", cell0, [cell1], postfix=True)

    # build theory
    l = Variable("l", cell0)
    m = Variable("m", cell0)
    n = Variable("n", cell0)
    o = Variable("o", cell0)
    f = Variable("f", cell1) # m <-- l
    g = Variable("g", cell1) # n <-- m
    h = Variable("h", cell1) # o <-- n

    Rewrite( l.identity*l.identity, l.identity )
    Rewrite( f*l.identity, f )
    Rewrite( m.identity*f, f )
    Rewrite( (g*f).src, f.src )
    Rewrite( (g*f).tgt, g.tgt )
    Rewrite( l.identity.src, l )
    Rewrite( l.identity.tgt, l )
    Equation( (h*g)*f, h*(g*f) )

    Rewrite( l.identity.inv, l.identity )
    Rewrite( f*f.inv, f.tgt.identity )
    Rewrite( f.inv*f, f.src.identity )
    Rewrite( f.inv.inv, f )
    Rewrite( f.inv.src, f.tgt )
    Rewrite( f.inv.tgt, f.src )
    #Equation( (f*f).src, (f*f).tgt ) # ??

    Equation( (g*f).inv, f.inv * g.inv )

    return theory


def test_rewrite_category_theory():

    theory = build_category_theory()

    Const = theory.Const
    cell0 = theory.cell0
    cell1 = theory.cell1

    X, Y, Z, W = [Const(name, cell0) for name in "XYZW"]
    f, g, h = [Const(name, cell1) for name in "fgh"]
    iso = Const("iso", cell1)

    assert h*(g*f) == (h*g)*f
    assert f*X.identity == f
    assert Y.identity*f == f
    assert f*f.src.identity == f
    assert f.tgt.identity*f == f

    assert (g*f).src == f.src
    assert (g*f).tgt == g.tgt

    assert (h*g*f).src == f.src
    assert (h*g*f).tgt == h.tgt

    assert X.identity.src == X
    assert X.identity.tgt == X

    assert X.identity.inv == X.identity
    assert f.src.identity.inv == f.src.identity

    assert distinct([f, g, h])

    assert iso * iso.inv == iso.tgt.identity
    assert f * iso * iso.inv == f


def build_bicategory_theory():

    cell0 = Sort("cell0") # object's
    cell1 = Sort("cell1") # morphism's
    cell2 = Sort("cell2") # 2-morphism's

    theory = build_category_theory(cell1, cell2)
    theory.cell0 = cell0
    theory.cell1 = cell1
    theory.cell2 = cell2

    Operator = theory.Operator
    Const = theory.Const
    Variable = theory.Variable
    Rewrite = theory.Rewrite
    Equation = theory.Equation
   
    Operator("identity", cell1, [cell0], postfix=True)
    Operator("<<",       cell1, [cell1, cell1], inline=True)
    Operator("<<",       cell2, [cell2, cell2], inline=True)

    Operator("src", cell0, [cell1], postfix=True)
    Operator("tgt", cell0, [cell1], postfix=True)

    Operator("lunitor",  cell2, [cell1],               postfix=True)
    Operator("runitor",  cell2, [cell1],               postfix=True)
    reassoc = Operator("reassoc",  cell2, [cell1, cell1, cell1], )
    theory.reassoc = reassoc


    # -----------------------------------------
    # build bicategory theory

    l = Variable("l", cell0)
    m = Variable("m", cell0)
    n = Variable("n", cell0)
    o = Variable("o", cell0)

    A = A0 = Variable("A0", cell1) # m <--A-- l
    A1     = Variable("A1", cell1) # m <--A-- l
    A2     = Variable("A2", cell1) # m <--A-- l
    B = B0 = Variable("B0", cell1) # n <--B-- m
    B1     = Variable("B1", cell1) # n <--B-- m
    B2     = Variable("B2", cell1) # n <--B-- m
    C = C0 = Variable("C0", cell1) # o <--C-- n
    C1     = Variable("C1", cell1) # o <--C-- n
    D      = Variable("D0", cell1) # p <--D-- o

    f = f0 = Variable("f0", cell2) #  A1 <----- A0
    f1     = Variable("f1", cell2) #  A2 <----- A1
    g = g0 = Variable("g0", cell2) #  B1 <----- B0
    g1     = Variable("g1", cell2) #  B2 <----- B1
    h      = Variable("n1", cell2) #  C1 <----- C0

    Equation( f.src.src, f.tgt.src )
    Equation( f.src.tgt, f.tgt.tgt )

    Rewrite( l.identity.src, l )
    Rewrite( l.identity.tgt, l )

    Rewrite( (B<<A).src, A.src )
    Rewrite( (B<<A).tgt, B.tgt )

    #Equation( (g<<f).src, g.src<<f.src )
    Rewrite( (g<<f).src, g.src<<f.src )

    Equation( (B<<A).identity, B.identity << A.identity )

    Rewrite( (g1*g0) << (f1*f0) , (g1 << f1) * (g0 << f0) )
    Rewrite( (g1 << f1) * (g0 << f0), (g1*g0)<<(f1*f0) )

    #Equation( (g<<f).inv, g.inv << f.inv ) # i don't think it's needed...?

    # unitors
    Rewrite( A.lunitor.tgt, A )
    Rewrite( A.lunitor.src, A.tgt.identity << A )
    Rewrite( A.runitor.tgt, A )
    Rewrite( A.runitor.src, A << A.src.identity )

    # naturality
    Equation( f * A0.lunitor, f.tgt.lunitor * (A0.tgt.identity.identity << f) )
    Equation( f * A0.runitor, f.tgt.runitor * (f << A0.src.identity.identity) )
    Equation( l.identity.lunitor, l.identity.runitor )

    # reassoc
    Rewrite( reassoc(C, B, A).tgt, C<<(B<<A) )
    Rewrite( reassoc(C, B, A).src, (C<<B)<<A )

    # naturality
    #Equation( reassoc(C1, B1, A1) * ( (h<<g)<<f ), ( h << (g<<f) ) * reassoc(C0, B0, A0) ) # fails because of free variables!
    Equation( reassoc(h.tgt, g.tgt, f.tgt) * ( (h<<g)<<f ), ( h << (g<<f) ) * reassoc(h.src, g.src, f.src) )

    # pentagon equation:
    # lhs:                rhs:
    #  ((D<<C)<<B)<<A     ((D<<C)<<B)<<A
    #                      (D<<C)<<(B<<A)
    #  (D<<(C<<B))<<A
    #  D<<((C<<B)<<A)
    #  D<<(C<<(B<<A))       D<<(C<<(B<<A))
    lhs = reassoc(D, C, B) << A.identity
    lhs = reassoc(D, C<<B, A) * lhs
    lhs = D.identity << reassoc(C, B, A) * lhs
    rhs = reassoc(D<<C, B, C)
    rhs = reassoc(D, C, B<<A) * rhs
    Equation(lhs, rhs)

    # triangle equation(s)
    Equation( B.runitor << A.identity, (B.identity << A.lunitor)*reassoc(B, B.src.identity, A) )
    Rewrite( (B<<A).lunitor, B.lunitor << A.identity ) # do we need this?
    Rewrite( (B<<A).runitor, B.identity << A.runitor )

    return theory


def test_rewrite_bicategory_theory():

    theory = build_bicategory_theory()
    cell0, cell1, cell2 = theory.cell0, theory.cell1, theory.cell2
    Const = theory.Const
    reassoc = theory.reassoc

    Equation = theory.Equation
    Rewrite = theory.Rewrite

    # -----------------------------------------
    # test theory

#    l = Const("l", cell0)
#    m = Const("m", cell0)
#    n = Const("n", cell0)
#    o = Const("o", cell0)
#
#    A = A0 = Const("A0", cell1) # m <--A-- l
#    A1     = Const("A1", cell1) # m <--A-- l
#    A2     = Const("A2", cell1) # m <--A-- l
#    B = B0 = Const("B0", cell1) # m <--B-- l
#    B1     = Const("B1", cell1) # m <--B-- l
#    B2     = Const("B2", cell1) # m <--B-- l

    f = f0 = Const("f0", cell2) #  A1 <----- A0
    f1     = Const("f1", cell2) #  A2 <----- A1
    f2     = Const("f2", cell2) #  A3 <----- A2
    f22    = Const("f22", cell2)#  A3 <----- A2
    g = g0 = Const("g0", cell2) #  B1 <----- B0
    g1     = Const("g1", cell2) #  B2 <----- B1
    g2     = Const("g2", cell2) #  B2 <----- B1
    h = h0 = Const("h0", cell2) #  C1 <----- C0
    k = k0 = Const("k0", cell2) #  D1 <----- D0
    E = Const("E", cell1)

    iso = Const("iso", cell2)

    A1, A0, A = f0.tgt, f0.src, f0.src
    A2, _A1   = f1.tgt, f1.src
    A3, _A2   = f2.tgt, f2.src

    B1, B0, B = g0.tgt, g0.src, g0.src
    B2, _B1   = g1.tgt, g1.src
    B3, _B2   = g2.tgt, g2.src

    C1, C0, C = h.tgt, h.src, h.src
    D1, D0, D = k.tgt, k.src, k.src

    assert distinct([f0, f1, f2, f22, g0, g1, g2, h0, k0])
    assert distinct([A0, A1, B0, B1, C0, C1, D0, D1])

    m, l = A0.tgt, A0.src
    assert A1.tgt == m
    assert A1.src == l

    Equation( A0.tgt, B0.src )
    assert A0.tgt == B0.src

    Equation( A1, _A1 )
    Equation( A2, _A2 )
    Equation( B1, _B1 )
    Equation( B2, _B2 )

    assert f2*(f1*f0) == (f2*f1)*f0
    assert f*A0.identity == f
    assert A1.identity*f == f
    assert f*f.src.identity == f
    assert f.tgt.identity*f == f

    assert A.identity.src == A
    assert A.identity.tgt == A

    assert (f1*f0).src == f0.src
    assert (f1*f0).tgt == f1.tgt

    assert A.identity.inv == A.identity
    assert f.src.identity.inv == f.src.identity

    assert iso * iso.inv == iso.tgt.identity
    assert f * iso * iso.inv == f

    A.lunitor * A.lunitor.inv == A.identity
    #theory.dump( A.lunitor * A.lunitor.inv )
    assert A.lunitor * A.lunitor.inv == A.identity
    assert f * f.src.lunitor == f.tgt.lunitor * (f.tgt.tgt.identity.identity << f)

    lhs = A1.lunitor * (A1.tgt.identity.identity << f) * A0.lunitor.inv
    assert lhs == f

    assert B.identity << A.identity == (B<<A).identity
    assert B.lunitor << A.identity != (B<<A).identity
    assert B.lunitor << A.identity == (B<<A).lunitor

    BA = B<<A
    CB = C<<B
    assert BA.identity.src == BA

    # triangle equation
    assert B.runitor << A.identity == (B.identity << A.lunitor) * reassoc(B, B.src.identity, A)
    assert CB.runitor << A.identity == (CB.identity << A.lunitor) * reassoc(CB, B.src.identity, A)

    I = l.identity
    assert I == I.src.identity
    assert I.runitor << I.identity == (I.identity << I.lunitor) * reassoc(I, I, I)
    #assert (B.runitor << A.identity)*reassoc(B, B.src.identity, A).inv == B.identity << A.lunitor # FAIL

    #Theory.DEBUG = True
    #Expr.DEBUG = True

    # naturality
    assert reassoc(C1, B1, A1) * ( (h<<g)<<f ) == ( h << (g<<f) ) * reassoc(C0, B0, A0)
    assert A.identity.tgt == A
    #assert reassoc(C1, B1, A) * ( (h<<g)<<A.identity ) == ( h << (g<<A.identity) ) * reassoc(C0, B0, A0) # FAIL


    # -----------------------

    f21 = f2*f1
    f210 = f21*f0
    f10 = f1*f0
    rhs = f2*f10
    assert f210 == (f2*f1)*f0
    assert f210 == f2*(f1*f0)

    assert f0.sort == cell2

    assert f2 != f22
    assert f2*f1 != f22*f1

    ff = f1*f0
    assert ff == f1*f0

    assert A1.identity*f0 == f0*A0.identity
    assert (f2*f1)*f0 == f2*(f1*f0)

    assert (g0<<f0).src == g0.src << f0.src

    lhs = g0 << f0
    i = g0.src.identity << f0.src.identity
    assert lhs == lhs * i

    i = (g0.src << f0.src).identity
    assert lhs == lhs * i

    gf = g<<f
    lhs = gf * gf.src.lunitor 
    rhs = gf.tgt.lunitor * (gf.src.tgt.identity.identity << gf)
    assert lhs == rhs

    assert A1.lunitor * A1.lunitor.inv * f == f
    assert A1.lunitor is A1.lunitor

    gf = g0<<f0
    assert gf == g0<<f0

    assert gf.src == B0<<A0
    assert gf * BA.identity == gf

    lhs = (g1*g0) << (f1*f0)
    rhs = (g1<<f1) * (g0<<f0)
    assert lhs == rhs

    g210 = g2*g1*g0
    f210 = f2*f1*f0
    lhs = g210 << f210

    rhs = (g2<<f2) * (g1<<f1) * (g0<<f0)
    assert lhs == rhs

    rrhs = (g2<<f2) * ((g1<<f1) * (g0<<f0))
    assert rhs == rrhs

    def test_identity(f):
        assert f * f.src.identity == f
        assert f.tgt.identity * f == f
    test_identity(g<<f)
    test_identity((g*g.src.lunitor)<<f)
    test_identity((g1*g0)<<f)

    def test_assoc(f2, f1, f0):
        assert (f2*f1)*f0 == f2*(f1*f0)
    test_assoc(f2, f1, f0)
    test_assoc(g2<<f2, g1<<f1, g0<<f0)

    def test_iso(f, name):
        #i = cat.Iso(f.src, f.src, name)
        i = Const(name, cell2)
        Equation( i.src, f.src )
        Equation( i.tgt, f.src )
        assert i * i.inv == i.src.identity # FAIL
        #assert i * i.inv == f.src.identity # FAIL
        assert f * i * i.inv == f

    test_iso(f, "iso_f")
    test_iso(g<<f, "iso_gf")
    test_iso((g1*g0)<<f, "iso_ggf")
    test_iso(g0<<(f1*f0), "iso_gff")
    test_iso((g1*g0)<<(f1*f0), "iso_gff")
    test_iso((g*g.src.lunitor)<<f, "iso_glf")

    # naturality of unitors
    def test_unitors(f):
        lhs = f * f.src.lunitor 
        rhs = f.tgt.lunitor * (f.src.tgt.identity.identity << f)
        assert lhs == rhs
        lhs = f * f.src.runitor 
        rhs = f.tgt.runitor * (f << f.src.src.identity.identity)
        assert lhs == rhs

    test_unitors(f)
    test_unitors(g)
    test_unitors(g<<f)

    lhs = g * g.src.lunitor 
    rhs = g.tgt.lunitor * (g.src.tgt.identity.identity << g)
    assert lhs == rhs
    lhs = lhs << f
    rhs = rhs << f
    assert lhs == rhs

    gf = g<<f
    assert gf.src == g.src << f.src
    #theory.dump( gf.src )
    #theory.dump( gf.src.lunitor )
    assert gf.src.lunitor == (g.src << f.src).lunitor
    I = g.src.tgt.identity
    assert gf.src.lunitor.src == I << (g.src << f.src)

    if 0:
        reassoc = ((I<<g.src)<<f.src).reassoc
        #test_assoc((g*g.src.lunitor) << f, reassoc.inv, reassoc)
        lhs = gf * gf.src.lunitor 
        assert lhs == ((g*g.src.lunitor) << f) * reassoc.inv
        assert reassoc.inv * reassoc == reassoc.src.identity
        assert lhs*reassoc == ((g*g.src.lunitor) << f) * reassoc.inv * reassoc
        test_assoc((g*g.src.lunitor) << f, reassoc.inv, reassoc) # FAIL !!!
        assert lhs*reassoc == ((g*g.src.lunitor) << f) * (reassoc.inv * reassoc)
        assert lhs*reassoc == ((g*g.src.lunitor) << f) * reassoc.src.identity
        test_identity( ((g*g.src.lunitor) << f) )
        assert lhs*reassoc == ((g*g.src.lunitor) << f)
        assert lhs == ((g.tgt.lunitor*(I.identity << g)) << f) * reassoc.inv
        rhs = gf.tgt.lunitor * (gf.src.tgt.identity.identity << gf)
        #assert rhs * reassoc 
    
        #assert lhs == rhs # FAIL 

    # _associators
    lhs = (C << B) << A
    rhs = C << (B << A)
    assert reassoc(C, B, A).src == lhs
    assert reassoc(C, B, A).tgt == rhs

    cell = (D<<C) << (B<<A)
    assert reassoc(D, C, B<<A).src == cell

    # naturality of reassoc
    hgf = (h0 << g0) << f0
    assert reassoc(h0.tgt, g0.tgt, f0.tgt) * hgf == (h0<<(g0<<f0)) * reassoc(h0.src, g0.src, f0.src)

    # triangle equation
    def test_triangle(B, A):
        lhs = B.runitor << A.identity
        rhs = B.identity << A.lunitor
        rhs = rhs * reassoc(B, B.src.identity, A)
        #rhs = rhs * lhs.src.reassoc
        assert lhs == rhs
    test_triangle(B, A)
    test_triangle(C<<B, A)
    test_triangle(C, B<<A)

    # fooling around with identity's
    l = Const("l", cell0)
    I = l.identity
    i = I.identity
    assert I.lunitor == I.runitor
    mul = I.lunitor
    comul = mul.inv
    assert comul == I.runitor.inv
    assert mul * comul == i
    ii = i<<i
    #assert comul * mul == ii # FAIL
    assert comul * mul == mul.src.identity
    assert mul.src == I<<I
    #assert mul.src.identity == (I<<I).identity # FAIL
    theory.dump(mul.src)
    theory.dump(mul.src.identity)
    #assert mul.src.identity == (I.identity<<I.identity) # FAIL
    lhs = mul*(mul << i)
    rhs = mul*(i<<mul) 
    a = reassoc(I, I, I)
    rhs = rhs * a
    #assert lhs == rhs # FAIL
    #assert (mul << i) * a.inv * (i << comul) == ii # FAIL
    #assert (i<<mul) * a * (comul << i) == ii # FAIL

    # pentagon equation
    def test_pentagon(D, C, B, A):
        lhs = reassoc(D, C, B) << A.identity
        lhs = reassoc(D, C<<B, A) * lhs
        lhs = D.identity << reassoc(C, B, A) * lhs
        rhs = reassoc(D<<C, B, C)
        rhs = reassoc(D, C, B<<A) * rhs
        assert lhs == rhs
    test_pentagon(D, C, B, A)
    test_pentagon(E<<D, C, B, A)
    test_pentagon(E, D<<C, B, A)

    assert distinct([f0, f1, f2, f22, g0, g1, g2, h0, k0])
    assert distinct([A0, A1, B0, B1, C0, C1, D0, D1])

    # Eckman - Hilton

    m = Const("m", cell0)
    U = m.identity
    u = U.identity
    f3 = Const("f3", cell2)
    g3 = Const("g3", cell2)
    Equation( f3.src, U )
    Equation( f3.tgt, U )
    Equation( g3.src, U )
    Equation( g3.tgt, U )
    assert f3!=g3

    assert f3*u == f3 == u*f3
    assert g3*u == g3 == u*g3
    assert f3<<g3 == (f3*u) << (u*g3)
    assert         (f3*u) << (u*g3) == (f3<<u)*(u<<g3)
    #assert                           (f3<<u)*(u<<g3) == ???
    #assert f3*g3 == g3*f3 # FAIL
    # not yet...

    

if __name__ == "__main__":

    print("\n\n")
    test_base_theory()
    test_theory()
    test_monoid_theory()
    test_group_theory()
    test_rewrite_category_theory()
    test_rewrite_bicategory_theory()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







