#!/usr/bin/env python3

"""
Minimal theorem proving using term unification...
"""

from argv import argv



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
        return Op.eq(self, other)

    def __eq__(self, other):
        assert 0, "fail"

    def __mul__(self, other):
        return Op.mul(self, other)

    def then(self, other):
        assert isinstance(other, Expr)
        return Sequent([self], other)

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

Op.eq = Op("==", 2, True)
Op.mul = Op("*", 2, True)


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
        elif op.infix:
            assert op.arity == 2
            s = "(%s%s%s)"%(self.args[0], op, self.args[1])
        else:
            s = "%s(%s)"%(op, ', '.join(str(a) for a in self.args))
        return s
    __repr__ = __str__

    def __eq__(self, other):
        assert 0

    def getvars(self):
        vs = set()
        for arg in self.args:
            vs.update(arg.getvars())
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


class Type(object):
    def __init__(self, name):
        self.name = name
        self.has = Op("%s.has"%name, 1) # unary Op, membership
        self.ns = {} # _namespace

    def __str__(self):
        return self.name

    def __getattr__(self, name):
        assert name[:2] != "__", name
        ns = self.ns
        v = ns.get(name)
        if v is None:
            v = Variable(name)
            ns[name] = v
        return v

    def const(self, name):
        return Op(name, 0)()


class Variable(Expr):
    unopcache = {} # 
    def __init__(self, name, tp=None):
        assert type(name) is str
        self.name = name
        self.tp = tp
        key = name
        Expr.__init__(self, key)

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __getattr__(self, name):
        assert name[:2] != "__", name
        op = self.unopcache.get(name)
        if op is None:
            op = Op(name, 1)
            self.unopcache[name] = op
        return op(self)

#    def __eq__(self, other):
#        assert isinstance(other, Expr)
#        if self.tp is None:
#            assert 0
#            return NotImplemented
#        return self.tp.eq(other)

    def getvars(self):
        return {self.name}

    def find(self, v):
        assert isinstance(v, Variable)
        return self is v

    def subs(self, send):
        return send.get(self.name, self)


class Sequent(object):
    def __init__(self, items=[], rhs=None):
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

    def __and__(self, other):
        assert self.rhs is None
        if isinstance(other, Sequent):
            assert other.rhs is None
            return Sequent(self.items + other.items)
        assert isinstance(other, Expr)
        return Sequent(self.items + [other])

    def then(self, other):
        assert isinstance(other, Expr)
        assert self.rhs is None, str(self)
        return Sequent(self.items, other)

    def __eq__(self, other):
        # Crazy hack to get __eq__ to re-associate !
        # See: https://docs.python.org/3/reference/expressions.html#operator-precedence
        if isinstance(other, Expr):
            other = Sequent([other])
        assert self.rhs is None
        assert other.rhs is None
        lhs = self.items
        rhs = other.items
        assert len(lhs)
        assert len(rhs)
        return Sequent(lhs[:-1] + [lhs[-1]==rhs[0]] + rhs[1:])

    #def __hash__(self):
    #    return hash(self.key)

    def __str__(self):
        items = ','.join(str(e) for e in self.items)
        s = "%s ⊢ %s"%(items, self.rhs)
        return s


#class Type(object):
#    def __init__(self, seqs=[]):
#        self.seqs = list(seqs)
#        has = Op("has", 1)
#        eq = Op("eq", 2)
#        a, b, c = [Variable(v) for v in 'abc']
#        self.append( Sequent([has(a)], eq(a, a)) )
#        self.append( Sequent([has(a), has(b), eq(a, b)], eq(b, a)) )
#        self.append( Sequent([has(a), has(b), has(c), eq(a, b), eq(b, c)], eq(a, c)) )
#        self.has = has
#        self.eq = eq
#
##    def has(self, var):
##        return self.op(var)
#
#    def append(self, seq):
#        self.seqs.append(seq)
#
#
#class Monoid(Type):
#    def __init__(self):
#        Type.__init__(self)
#        mul = Op("mul", 2)
#        one = Op("one", 0)()
#        has = self.has
#        eq = self.eq
#        f, g, h = [Variable(v) for v in 'fgh']
#        self.append( Sequent([], has(one)) )
#        self.append( Sequent([has(f), has(g)], has(mul(f, g))) )
#        self.append( Sequent([has(f), has(g), has(h)], eq( mul(f, mul(g, h)), mul(mul(f, g), h))) )
#        self.append( Sequent([has(f)], eq(f, mul(1, f))) )
#        self.append( Sequent([has(f)], eq(f, mul(f, 1))) )
#    

        

def main():

    # ----------------------------------------

    tp = Type("T")

    w, x, y, z = tp.w, tp.x, tp.y, tp.z

    f = Op("f", 3)
    g = Op("g", 1)

    lhs = f(x, y, g(z))
    rhs = f(x, y, w)

    assert str(lhs) == "f(x, y, g(z))"
    assert lhs.getvars() == set("xyz")

    send = Expr.unify(lhs, rhs)
    assert str(send) == "{'w': g(z)}"

    send = Expr.unify(x, g(x))
    assert send is None

    lhs = f(x, y, z)
    rhs = f(y, x, w)

    assert str(Expr.unify(lhs, rhs)) == "{'x': x, 'y': x, 'z': w}"

    # ----------------------------------------

    eq = Op.eq
    yes = Sequent()

    # ----------------------------------------
    # Equality 

    tp = Type("T")
    a, b, c = tp.a, tp.b, tp.c
    seq = tp.has(a).then(a.eq(a))
    seq = (tp.has(a) & tp.has(b) & a.eq(b)).then(b.eq(a))
    seq = (tp.has(a) & tp.has(b) & tp.has(c) & a.eq(b) & b.eq(c)).then(a.eq(c))

    # ----------------------------------------
    # Monoid 

    tp = Type("M")
    f, g, h = tp.f, tp.g, tp.h
    one = tp.const("1")
    print( yes.then( tp.has(one) ) )
    print( (tp.has(f) & tp.has(g)).then( tp.has( f*g ) ) )
    print( (tp.has(f) & tp.has(g) & tp.has(h) ).then( eq( (f*g)*h, f*(g*h) ) ) )
    print( tp.has(f).then( (f*one).eq(f) ) )
    print( tp.has(f).then( (one*f).eq(f) ) )
    
    # ----------------------------------------
    # Category ... 

    tp = Type("C")
    A, B, C = tp.A, tp.B, tp.C
    print( tp.has(A).then( tp.has( A.src ) ) )
    print( (tp.has(A) & tp.has(B) & A.tgt.eq(B.src) ).then( tp.has( B*A ) ) )

    # etc.

    # TODO: equality is extensive
    # TODO: proof's , etc.
    

    


if __name__ == "__main__":

    main()

    print("OK\n")

    



