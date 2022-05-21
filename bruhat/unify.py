#!/usr/bin/env python3

"""
Term Unification...
"""

from argv import argv


class Expr(object):
    pass


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
        assert len(args) == self.arity
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
        _send = unify(left, right, depth+1)
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

    send = unify(lhs, rhs)
    assert send == {w : g(z)}

    send = unify(x, g(x))
    assert send is None

    lhs = f(x, y, z)
    rhs = f(y, x, w)

    print(unify(lhs, rhs))




if __name__ == "__main__":

    main()

    print("OK\n")

    



