#!/usr/bin/env python3

"""
Term Unification...
"""

from argv import argv


class Expr(object):

    @staticmethod
    def unify(lhs, rhs, depth=0):
        indent = "  "*depth
        #print(indent+"unify:", lhs, "&", rhs)
        if lhs==rhs:
            #print(indent, "<== {}")
            return {}
        if isinstance(lhs, Variable):
            if lhs!=rhs and rhs.find(lhs):
                return None # fail
            #print(indent, "<== {%s : %s}"%(lhs, rhs))
            return {lhs : rhs}
        elif isinstance(rhs, Variable):
            if lhs!=rhs and lhs.find(rhs):
                return None # fail
            #print(indent, "<== {%s : %s}"%(rhs, lhs))
            return {rhs : lhs}
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
                if v == _v:
                    continue # agreement
                if v is not None: # collision
                    return None # fail
                # update send dict
                for (k,v) in list(send.items()):
                    send[k] = v.subs({_k : _v})
                # add new key:value pair
                send[_k] = _v
        return send



class Variable(Expr):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __hash__(self):
        #return hash(self.name)
        return id(self)

    def getvars(self):
        return {self}

    def find(self, v):
        assert isinstance(v, Variable)
        return self==v

    def subs(self, send):
        return send.get(self, self)


class Op(object):
    def __init__(self, name, arity=0):
        self.name = name
        self.arity = arity

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __call__(self, *args):
        assert len(args) == self.arity, "%s%s fail" % (self, args)
        return Apply(self, *args)

    def __hash__(self):
        return id(self)


class Apply(Expr):
    def __init__(self, op, *args):
        self.op = op
        self.args = args

    def __str__(self):
        return "%s(%s)"%(self.op, ', '.join(str(a) for a in self.args))
    __repr__ = __str__

    def __eq__(self, other):
        if other is None:
            return False
        assert isinstance(other, Expr)
        return str(self) == str(other) # lazy me

    def __hash__(self):
        return hash(str(self)) # um... 

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


class Sequent(object):
    def __init__(self, lhs, rhs):
        assert isinstance(rhs, Expr)
        for expr in lhs:
            assert isinstance(expr, Expr)
        self.lhs = set(lhs)
        self.rhs = rhs
        lhs = list(lhs)
        lhs.sort(key=str)
        lhs = tuple(lhs)
        self.key = (lhs, rhs)

    def __eq__(self, other):
        return self.lhs==other.lhs and self.rhs==other.rhs

    def __hash__(self):
        return hash(self.key)

    def __str__(self):
        lhs = ','.join(str(e) for e in self.key[0])
        s = "%s ⊢ %s"%(lhs, self.rhs)
        return s


class Type(object):
    def __init__(self, seqs=[]):
        self.seqs = list(seqs)
        has = Op("has", 1)
        eq = Op("eq", 2)
        a, b, c = [Variable(v) for v in 'abc']
        self.append( Sequent([has(a)], eq(a, a)) )
        self.append( Sequent([has(a), has(b), eq(a, b)], eq(b, a)) )
        self.append( Sequent([has(a), has(b), has(c), eq(a, b), eq(b, c)], eq(a, c)) )
        self.has = has
        self.eq = eq

#    def has(self, var):
#        return self.op(var)

    def append(self, seq):
        self.seqs.append(seq)


class Monoid(Type):
    def __init__(self):
        Type.__init__(self)
        mul = Op("mul", 2)
        one = Op("one", 0)()
        has = self.has
        eq = self.eq
        f, g, h = [Variable(v) for v in 'fgh']
        self.append( Sequent([], has(one)) )
        self.append( Sequent([has(f), has(g)], has(mul(f, g))) )
        self.append( Sequent([has(f), has(g), has(h)], eq( mul(f, mul(g, h)), mul(mul(f, g), h))) )
        self.append( Sequent([has(f)], eq(f, mul(1, f))) )
        self.append( Sequent([has(f)], eq(f, mul(f, 1))) )
    
        

def main():

    w, x, y, z = [Variable(c) for c in 'wxyz']

    f = Op("f", 3)
    g = Op("g", 1)

    lhs = f(x, y, g(z))
    rhs = f(x, y, w)

    assert lhs == f(x, y, g(z))
    assert lhs != f(x, y, z)
    assert str(lhs) == "f(x, y, g(z))"
    assert lhs.getvars() == {x, y, z}

    send = Expr.unify(lhs, rhs)
    assert send == {w : g(z)}

    send = Expr.unify(x, g(x))
    assert send is None

    lhs = f(x, y, z)
    rhs = f(y, x, w)

    #print(Expr.unify(lhs, rhs))

    T = Type()

    M = Monoid()
    for seq in M.seqs:
        print(seq)


if __name__ == "__main__":

    main()

    print("OK\n")

    



