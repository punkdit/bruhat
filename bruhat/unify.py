#!/usr/bin/env python3

"""
Minimal theorem proving using term unification...

Failures:
    tried to use __eq__ for Expr's to denote boolean equality Op( ),
    but == operator precedence is too low, it's impossible to
    get this right "a==b & b==c" it get's turned into a chain of
    equality tests: "a==(b&b)==c" and then "a==b&b and b&b==c" which
    cannot be overloaded because of the "and".
"""




class Expr(object):

    # use .key for syntactic comparison, hash, etc.
    def __init__(self, key):
        self.key = key

    def __and__(self, other):
        assert isinstance(other, Expr)
        return Sequent([self, other])

    def eq(self, other):
        if isinstance(other, Sequent):
            return Sequent([self]) == other # yikes
        assert isinstance(other, Expr)
        return eq(self, other)

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def __mul__(self, other):
        return mul(self, other)

    def then(self, other, comment=""):
        assert isinstance(other, Expr)
        return Sequent([self], other, comment)

    @staticmethod
    def unify(lhs, rhs, depth=0):
        indent = "  "*depth
        #print(indent+"unify:", lhs, "&", rhs)
        if lhs.key==rhs.key:
            #print(indent, "<== {}")
            return {}
        if isinstance(lhs, Variable):
            if lhs.key!=rhs.key and rhs.find(lhs):
                return None # fail
            #print(indent, "<== {%s : %s}"%(lhs, rhs))
            return {lhs.key : rhs}
        elif isinstance(rhs, Variable):
            if lhs.key!=rhs.key and lhs.find(rhs):
                return None # fail
            #print(indent, "<== {%s : %s}"%(rhs, lhs))
            return {rhs.key : lhs}
        assert isinstance(lhs, Apply)
        assert isinstance(rhs, Apply)
        if lhs.op != rhs.op:
            return None # fail
        send = {}
        for left, right in zip(lhs.args, rhs.args):
            _send = Expr.unify(left, right, depth+1)
            if _send is None:
                return None # fail
            for _k,_v in _send.items():
                v = send.get(_k)
                if v is _v:
                    continue # agreement
                assert _v is not None
                if v is not None and v.key == _v.key:
                    continue # agreement
                elif v is not None: # collision
                    return None # fail
                # update send dict
                for (k,v) in list(send.items()):
                    send[k] = v.subs({_k : _v})
                # add new key:value pair
                send[_k] = _v
        return send



class Op(object):
    def __init__(self, name, arity=0, infix=False):
        assert type(name) is str
        self.name = name
        self.arity = arity
        self.infix = infix

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __call__(self, *args):
        assert len(args) == self.arity, "%s%s fail" % (self, args)
        return Apply(self, *args)

eq = Op("==", 2, True)
mul = Op("*", 2, True)


class Apply(Expr):
    def __init__(self, op, *args):
        assert len(args) == op.arity
        self.op = op
        self.args = args
        key = (op.name,)+tuple(arg.key for arg in args)
        Expr.__init__(self, key)

    def __str__(self):
        op = self.op
        if op.arity==0:
            return op.name
        elif op.arity==2 and op.infix:
            s = "(%s%s%s)"%(self.args[0], op, self.args[1])
        elif op.arity==1 and op.infix:
            s = "%s.%s"%(self.args[0], op,)
        else:
            s = "%s(%s)"%(op, ', '.join(str(a) for a in self.args))
        return s
    __repr__ = __str__

    def __eq__(self, other):
        assert 0

    def all_vars(self):
        vs = set()
        for arg in self.args:
            vs.update(arg.all_vars())
        return vs

    def find(self, v):
        assert isinstance(v, Variable)
        for arg in self.args:
            if arg.find(v):
                return True
        return False

    def subs(self, send):
        op = self.op
        args = [arg.subs(send) for arg in self.args]
        return op(*args)


class Variable(Expr):
    def __init__(self, name, tp):
        assert type(name) is str
        self.name = name
        self.tp = tp
        #key = (tp, name) # we don't need to be this strict.. ??
        key = name 
        Expr.__init__(self, key)

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __getattr__(self, name):
        assert name[:2] != "__", name
        op = self.tp.get_op(name, 1, infix=True, partial=True)
        return op(self)

    def all_vars(self):
        return {self.name}

    def find(self, v):
        assert isinstance(v, Variable)
        return self is v

    def subs(self, send):
        return send.get(self.name, self)


class Sequent(object):
    def __init__(self, items=[], rhs=None, comment=""):
        # XXX use another class when rhs is None ???
        assert rhs is None or isinstance(rhs, Expr)
        for expr in items:
            assert isinstance(expr, Expr)
        self.items = list(items) # maintain the order given
        self.rhs = rhs
        # now build a canonical form:
        items = list(set(item.key for item in items)) # uniq
        items.sort(key=str) # ordered
        items = tuple(items)
        rhs = rhs.key if rhs is not None else None
        self.key = (items, rhs)
        self.comment = comment

    def __and__(self, other):
        assert self.rhs is None
        if isinstance(other, Sequent):
            assert other.rhs is None
            return Sequent(self.items + other.items)
        assert isinstance(other, Expr)
        return Sequent(self.items + [other])

    def then(self, other, comment=""):
        assert isinstance(other, Expr)
        assert self.rhs is None, str(self)
        return Sequent(self.items, other, comment)

#    def __eq__(self, other):
#        # Crazy hack to get __eq__ to re-associate !
#        # See: https://docs.python.org/3/reference/expressions.html#operator-precedence
#        if isinstance(other, Expr):
#            other = Sequent([other])
#        assert self.rhs is None
#        assert other.rhs is None
#        lhs = self.items
#        rhs = other.items
#        assert len(lhs)
#        assert len(rhs)
#        return Sequent(lhs[:-1] + [lhs[-1]==rhs[0]] + rhs[1:])

    #def __hash__(self):
    #    return hash(self.key)

    def __str__(self):
        items = ','.join(str(e) for e in self.items)
        s = "%s ⊢ %s"%(items, self.rhs)
        if self.comment:
            s = s + " # " + self.comment
        return s

yes = Sequent()
always = yes.then


class Type(object):
    def __init__(self, name):
        self.name = name
        self.has = Op("%s.has"%name, 1) # unary Op, membership
        self.vars = {} # _namespace for Variable's
        self.ops = {}
        self.axioms = []
        self.add_equality()

    def add(self, seq):
        self.axioms.append(seq)

    def add_equality(tp):
        a, b, c = tp.a, tp.b, tp.c
        tp.add( tp.has(a).then(a.eq(a), "(eq-reflexive)") )
        tp.add( (tp.has(a) & tp.has(b) & a.eq(b)).then(b.eq(a), "(eq-symmetric)") )
        tp.add( (tp.has(a) & tp.has(b) & tp.has(c) & a.eq(b) & b.eq(c)).then(a.eq(c), "(eq-transitive)") )

    def add_extensive(self, op):
        if op.arity==0:
            return # these are already extensive 
        # extensivity of eq
        lhs = yes
        us = [self.get_var("u_%d"%i) for i in range(op.arity)]
        vs = [self.get_var("v_%d"%i) for i in range(op.arity)]
        for u,v in zip(us, vs):
            lhs = lhs & self.has(v) & self.has(u) & eq(u, v)
        self.add( lhs.then( eq(op(*vs), op(*us)), "(%s-eq-extensive)"%op.name ) )

    def add_op(self, op, partial=False):
        assert isinstance(op, Op)

        if not partial:
            vs = [self.get_var("v_%d"%i) for i in range(op.arity)]
            lhs = yes
            for v in vs:
                lhs = lhs & self.has(v)
            self.add( lhs.then( self.has(op(*vs)), "(%s-endo)"%op.name ) )
        self.add_extensive(op)

    def __str__(self):
        return self.name

    def get_var(self, name):
        assert name[:2] != "__", name
        ns = self.vars
        v = ns.get(name)
        if v is None:
            v = Variable(name, self)
            ns[name] = v
        return v
    __getattr__ = get_var

    def get_op(self, name, arity=0, infix=False, partial=False):
        ops = self.ops
        if name in ops:
            op = ops[name]
            assert op.arity == arity
            return op
        op = Op(name, arity, infix)
        ops[name] = op
        self.add_op(op, partial)
        return op

    def const(self, name):
        op = self.get_op(name)
        return op()

    def dump(self):
        print("%s:"%self.name)
        for seq in self.axioms:
            print("\t%s"%seq)


def main():

    # ----------------------------------------

    tp = Type("T")

    w, x, y, z = tp.w, tp.x, tp.y, tp.z

    assert repr(w) == 'w'

    f = Op("f", 3)
    g = Op("g", 1)

    lhs = f(x, y, g(z))
    rhs = f(x, y, w)

    assert str(lhs) == "f(x, y, g(z))"
    assert lhs.all_vars() == set("xyz")

    send = Expr.unify(lhs, rhs)
    assert str(send) == "{'w': g(z)}", str(send)

    send = Expr.unify(x, g(x))
    assert send is None

    lhs = f(x, y, z)
    rhs = f(y, x, w)

    assert str(Expr.unify(lhs, rhs)) == "{'x': x, 'y': x, 'z': w}"

    # ----------------------------------------

    # ----------------------------------------
    # Equality 

    tp = Type("T")

    # ----------------------------------------
    # Monoid 

    tp = Type("M")
    f, g, h = tp.f, tp.g, tp.h
    one = tp.const("1")
    tp.add_op(mul)
    tp.dump()

    if 0:
        print( yes.then( tp.has(one) ) )
        print( (tp.has(f) & tp.has(g)).then( tp.has( f*g ) ) )
        print( (tp.has(f) & tp.has(g) & tp.has(h) ).then( eq( (f*g)*h, f*(g*h) ) ) )
        print( tp.has(f).then( (f*one).eq(f) ) )
        print( tp.has(f).then( (one*f).eq(f) ) )
    
    # ----------------------------------------
    # Category ... 

    tp = Type("C")
    A, B, C = tp.A, tp.B, tp.C
    tp.add( tp.has(A).then( tp.has( A.src ) ) )
    tp.add( (tp.has(A) & tp.has(B) & A.tgt.eq(B.src) ).then( tp.has( B*A ) ) )
    tp.add_extensive(mul) # composition is extensive while also being a partial function...

    tp.dump()

    # etc.

    # TODO: equality is extensive
    # TODO: proof's , etc.
    

    


if __name__ == "__main__":

    main()

    print("OK\n")

    



