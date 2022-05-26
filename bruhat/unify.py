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
    def __init__(self, tp, key):
        assert isinstance(tp, Theory)
        self.tp = tp
        self.key = key

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def __and__(self, other):
        assert isinstance(other, Expr)
        return Sequent([self, other])

    # --------------------------------------------------------
    # _Operators

    def __getattr__(self, name):
        assert name[:2] != "__", name
        #op = self.tp.get_op(name, 1, infix=True, partial=True)
        op = self.tp.get_op(name)
        return op(self)

    def eq(self, other):
        if isinstance(other, Sequent):
            return Sequent([self]) == other # yikes
        assert isinstance(other, Expr)
        return self.tp.eq(self, other)

    def __mul__(self, other):
        return self.tp.get_op("*")(self, other)

    def __add__(self, other):
        return self.tp.get_op("+")(self, other)

    def __sub__(self, other):
        return self.tp.get_op("-")(self, other)

    # etc. etc.

    def __matmul__(self, other):
        return self.tp.get_op("@")(self, other)

    def __lshift__(self, other):
        return self.tp.get_op("<<")(self, other)

    def __rshift__(self, other):
        return self.tp.get_op(">>")(self, other)

    # --------------------------------------------------------

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
    def __init__(self, tp, name, arity=0, infix=True):
        assert isinstance(tp, Theory)
        assert type(name) is str
        self.tp = tp
        self.name = name
        self.arity = arity
        self.infix = infix

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __call__(self, *args):
        #print("Op.__call__", self, args)
        assert len(args) == self.arity, "%s%s fail" % (self, args)
        return Apply(self, *args)


class Apply(Expr):
    def __init__(self, op, *args):
        assert len(args) == op.arity
        self.op = op
        self.args = args
        key = (op.name,)+tuple(arg.key for arg in args)
        Expr.__init__(self, op.tp, key)

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
    def __init__(self, tp, name):
        assert type(name) is str
        self.name = name
        self.tp = tp
        #key = (tp, name) # we don't need to be this strict.. ??
        key = name 
        Expr.__init__(self, tp, key)

    def __str__(self):
        return self.name
    __repr__ = __str__

    def all_vars(self):
        return {self.name}

    def find(self, v):
        assert isinstance(v, Variable)
        return self is v

    def subs(self, send):
        return send.get(self.name, self)


class Sequent(object):
    def __init__(self, src=[], tgt=None, comment=""):
        # XXX use another class when tgt is None ???
        assert tgt is None or isinstance(tgt, Expr)
        for expr in src:
            assert isinstance(expr, Expr)
        self.src = list(src) # maintain the order given
        self.tgt = tgt
        # now build a canonical form:
        src = list(set(item.key for item in src)) # uniq
        src.sort(key=str) # ordered
        src = tuple(src)
        tgt = tgt.key if tgt is not None else None
        self.key = (src, tgt)
        self.comment = comment

    def __and__(self, other):
        assert self.tgt is None
        if isinstance(other, Sequent):
            assert other.tgt is None
            return Sequent(self.src + other.src)
        assert isinstance(other, Expr)
        return Sequent(self.src + [other])

    def then(self, other, comment=""):
        assert isinstance(other, Expr)
        assert self.tgt is None, str(self)
        return Sequent(self.src, other, comment)

    def all_vars(self):
        vs = set()
        for e in self.src:
            vs.update(e.all_vars())
        if self.tgt is not None:
            vs.update(self.tgt.all_vars())
        return vs

    #def __hash__(self):
    #    return hash(self.key)

    def __str__(self):
        src = ','.join(str(e) for e in self.src)
        s = "%s ⊢ %s"%(src, self.tgt)
        if self.comment:
            s = s + " # " + self.comment
        return s

yes = Sequent()
always = yes.then


class Theory(object):
    def __init__(self, name):
        self.name = name
        self.vars = {} # _namespace for Variable's
        self.ops = {}
        self.seqs = []

#        self.has = Op(self, "%s.has"%name, 1) # unary Op, membership
#        self.eq = Op(self, "==", 2, True)
#        self.mul = Op(self, "*", 2, True)
##        self.add = Op(self, "+", 2, True)
#        self.tensor = Op(self, "@", 2, True)
#        self.lshift = Op(self, "<<", 2, True)
#        self.rshift = Op(self, ">>", 2, True)

        self.add_equality()

    def add(self, seq):
        assert isinstance(seq, Sequent)
        self.seqs.append(seq)

    def __getitem__(self, idx):
        return self.seqs[idx]

    def add_equality(tp):
        a, b, c = tp.a, tp.b, tp.c
        #tp.has = tp.add_op("%s.has"%tp.name, 1)
        tp.has = Op(tp, "%s.has"%tp.name, 1, False)
        tp.eq = Op(tp, "==", 2)
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
            lhs = lhs & self.has(v) & self.has(u) & self.eq(u, v)
        self.add( lhs.then( self.eq(op(*vs), op(*us)), "(%s-eq-extensive)"%op.name ) )

    def add_op(self, name, arity=0, infix=True, partial=False):
        assert name not in self.ops
        op = Op(self, name, arity, infix)
        self.ops[name] = op
        if not partial:
            vs = [self.get_var("v_%d"%i) for i in range(op.arity)]
            lhs = yes
            for v in vs:
                lhs = lhs & self.has(v)
            self.add( lhs.then( self.has(op(*vs)), "(%s-endo)"%op.name ) )
        self.add_extensive(op)
        return op

    def const(self, name):
        op = self.add_op(name)
        return op()

    def __str__(self):
        return self.name

    def get_var(self, name):
        assert name[:2] != "__", name
        ns = self.vars
        v = ns.get(name)
        if v is None:
            v = Variable(self, name)
            ns[name] = v
        return v
    __getattr__ = get_var

    def get_op(self, name):
        op = self.ops[name]
        return op

#    def get_op(self, name, arity=0, infix=False, partial=False):
#        ops = self.ops
#        if name in ops:
#            op = ops[name]
#            assert op.arity == arity
#            return op
#        op = Op(self, name, arity, infix)
#        ops[name] = op
#        self.add_op(op, partial)
#        return op

    def deduce(self, idx, jdx):
        left = self[idx]
        right = self[jdx]
        lvars = left.all_vars()
        rvars = right.all_vars()
        assert not (lvars & rvars), "TODO"
        for i, expr in enumerate(right.src):
            send = Expr.unify(left.tgt, expr)
            if send is None:
                continue
            src = [e.subs(send) for e in left.src + right.src[:i] + right.src[i+1:]]
            tgt = right.tgt.subs(send)
            seq = Sequent(src, tgt, "[%s] & [%s]"%(idx, jdx))
            self.add(seq)
            return seq

    def dump(self):
        print("%s:"%self.name)
        for idx, seq in enumerate(self.seqs):
            print("[%d]\t%s"%(idx, seq))


def main():

    # ----------------------------------------

    tp = Theory("T")

    w, x, y, z = tp.w, tp.x, tp.y, tp.z

    assert repr(w) == 'w'

    f = Op(tp, "f", 3, False)
    g = Op(tp, "g", 1, False)

    lhs = f(x, y, g(z))
    rhs = f(x, y, w)

    assert str(lhs) == "f(x, y, g(z))", str(lhs)
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

    tp = Theory("T")

    # ----------------------------------------
    # Monoid 

    tp = Theory("M")
    f, g, h = tp.f, tp.g, tp.h
    one = tp.const("1")
    tp.add_op("*", 2)

    seq = tp.deduce(3, 4)
    seq = tp.deduce(3, 6)

    tp.dump()

    # ----------------------------------------
    # Category, see Freyd-Scedrov 

    tp = Theory("C")
    tp.add_op("*", 2)
    tp.add_op("src", 1)
    tp.add_op("tgt", 1)
    eq = tp.eq

    A, B, C = tp.A, tp.B, tp.C
    tp.add( tp.has(A).then( tp.has( A.src ) ) )
    tp.add( tp.has(A).then( tp.has( A.tgt ) ) )
    tp.add( (tp.has(A) & tp.has(B) & A.src.eq(B.tgt) ).then( tp.has( A*B ) ) )
    tp.add( (tp.has(A) & tp.has(B) & tp.has(A*B) ).then( A.src.eq(B.tgt) ) )
    tp.add( tp.has(A).then( eq( A.src.tgt, A.src ) ) )
    tp.add( tp.has(A).then( eq( A.tgt.src, A.tgt ) ) )
    tp.add( tp.has(A).then( eq( A.tgt*A, A ) ) )
    tp.add( tp.has(A).then( eq( A*A.src, A ) ) )
    tp.add( (tp.has(A) & tp.has(B) & tp.has(A*B) ).then(
        eq( (A*B).tgt, (A*B.tgt).tgt )) )
    tp.add( (tp.has(A) & tp.has(B) & tp.has(A*B) ).then(
        eq( (A*B).src, (A.src*B).src )) )


    tp.dump()


    


if __name__ == "__main__":

    main()

    print("OK\n")

    



