#!/usr/bin/env python3

"""
Minimal theorem _proving using term unification...

Failures:
    tried to use __eq__ for Expr's to denote boolean equality Op( ),
    but == operator precedence is too low, it's impossible to
    get this right "a==b & b==c" it get's turned into a chain of
    equality tests: "a==(b&b)==c" and then "a==b&b and b&b==c" which
    cannot be overloaded because of the "and".
"""

import heapq
from random import randint
from string import ascii_letters

from bruhat.util import all_perms


class Expr(object):

    # use .key for syntactic comparison, hash, etc.
    def __init__(self, tp, key):
        assert isinstance(tp, Theory)
        self.tp = tp
        self.key = key

    def __eq__(self, other):
        if other is None:
            return False
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
    def unify(lhs, rhs, send={}, depth=0):
        indent = "  "*depth
#        print(indent+"unify:", lhs, "&", rhs)
        if send:
            lhs = lhs.subs(send)
            rhs = rhs.subs(send)
        if isinstance(rhs, Variable):
            lhs, rhs = rhs, lhs
        if isinstance(lhs, Variable):
            if lhs!=rhs and rhs.find(lhs):
                return None # fail (recursive)
#            print(indent, "<== {%s : %s}"%(lhs, rhs))
            assert lhs not in send
            if lhs != rhs:
                send = dict(send)
                for k,v in list(send.items()):
                    send[k] = v.subs({lhs:rhs})
                send[lhs] = rhs
            return send
        assert isinstance(lhs, Apply)
        assert isinstance(rhs, Apply)
        if lhs.op != rhs.op:
            return None # fail
        for left, right in zip(lhs.args, rhs.args):
            send = Expr.unify(left, right, send, depth=depth+1)
            if send is None:
                return None # fail
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
        return Apply(self, args)


class Apply(Expr):
    def __init__(self, op, args):
        assert len(args) == op.arity
        self.op = op
        self.args = args
        key = (op.name,)+args
        Expr.__init__(self, op.tp, key)
        self.size = sum([arg.size for arg in args], 1)

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

#    def specialize(self, other, iso=False):
#        # can we substitute into self to get other?
#        if isinstance(other, Variable):
#            return None # fail
#        assert isinstance(other, Apply)
#        if other.op != self.op:
#            return None # fail
#        send = {}
#        for l,r in zip(self.args, other.args):
#            _send = l.specialize(r, iso)
#            if _send is None:
#                return None # fail
#            for (_k,_v) in _send.items():
#                if _k in send and send[_k] != _v:
#                    return None # collision; fail
#                send[_k] = _v
#        return send

    def specialize(self, other, send={}, iso=False):
        # can we substitute into self to get other?
        if isinstance(other, Variable):
            return None # fail
        assert isinstance(other, Apply), type(other)
        if other.op != self.op:
            return None # fail
        for l,r in zip(self.args, other.args):
            send = l.specialize(r, send, iso)
            if send is None:
                return None # fail
        return send

    def subs(self, send): # hotspot 
        op = self.op
        args = [arg.subs(send) for arg in self.args]
        return op(*args)


class Variable(Expr):
    size = 1
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
        return {self}

    def find(self, v):
        assert isinstance(v, Variable)
        return self is v

    def subs(self, send):
        return send.get(self, self)

#    def specialize(self, other, iso=False):
#        # can we substitute into self to get other? 
#        if not iso or isinstance(other, Variable):
#            return {self:other}

    def specialize(self, other, send={}, iso=False):
        # can we substitute into self to get other? 
        if not iso or isinstance(other, Variable):
            value = send.get(self)
            if value == other:
                return send
            if value is None:
                send = dict(send)
                send[self] = other
                return send

#    def equiv(self, other):
#        # reversible specialize
#        if isinstance(other, Variable):
#            return {self : other}


class Sequent(object):
    def __init__(self, src=[], tgt=None, comment=None, parent=None):
        # XXX use another class when tgt is None ???
        assert tgt is None or isinstance(tgt, Expr)
        uniqsrc = set()
        self.src = [] # maintain the given order
        size = 0
        for expr in src:
            assert isinstance(expr, Expr)
            if expr not in uniqsrc:
                self.src.append(expr)
                uniqsrc.add(expr)
                size += expr.size
        if tgt is not None:
            size += tgt.size
        self.tgt = tgt
        # now build a canonical form:
        uniqsrc = list(uniqsrc)
        uniqsrc.sort(key=str) # ordered
        uniqsrc = tuple(uniqsrc)
        self.uniqsrc = uniqsrc
        self.key = (uniqsrc, tgt)
        self.comment = comment
        self.parent = parent
        self.size = size

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def weak_eq(self, other):
        #print("weak_eq:")
        #print("\t%s"%(self,))
        #print("\t%s"%(other,))
        if len(self.uniqsrc) != len(other.uniqsrc):
            return 
        assert other.tgt is not None
        subs = self.tgt.specialize(other.tgt, iso=True)
        if subs is None:
            return 
        n = len(self.uniqsrc)
        items = list(range(n))
        #print("weak_eq", subs)
        for perm in all_perms(items): # XXX slow & stupid, fix me
            #print("\t", perm)
            _subs = dict(subs)
            for i,j in enumerate(perm):
                lhs, rhs = self.uniqsrc[i], other.uniqsrc[j]
                _subs = lhs.specialize(rhs, _subs, iso=True)
                #print("weak_eq: specialize", lhs, rhs, _subs)
                if _subs is None:
                    break
            else:
                return _subs

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

    def src_vars(self):
        vs = set()
        for e in self.src:
            vs.update(e.all_vars())
        return vs

    def subs(self, send, comment=None):
        if comment is None:
            comment = self.comment
        return Sequent([e.subs(send) for e in self.src], self.tgt.subs(send), comment)

    #def __hash__(self):
    #    return hash(self.key)

    def __str__(self):
        src = ','.join(str(e) for e in self.src)
        s = "%s ⊢ %s"%(src, self.tgt)
        if self.comment:
            s = s + " # " + self.comment
        return s

    def rename(self):
        vs = self.all_vars()
        if not vs:
            return self
        tp = self.tgt.tp
        vs = list(vs)
        vs.sort(key = str)
        send = {}
        for i,v in enumerate(vs):
            send[v] = tp.get_var(ascii_letters[i])
        seq = self.subs(send)
        return seq

    def can_deduce(left, right):
        left = left.tgt
        lvars = left.all_vars()
        rvars = right.all_vars()
        ivars = lvars & rvars
        if ivars:
            tp = left.tp
            send = {v:tp.get_var(v.name+"_l") for v in ivars}
            left = left.subs(send)
        for i, expr in enumerate(right.src):
            send = Expr.unify(left, expr)
            if send is not None:
                return True

    def deduce(left, right, comment=None, verbose=False): # hotspot
        tp = left.tgt.tp
        lvars = left.all_vars()
        rvars = right.all_vars()
        if (lvars & rvars):
            ivars = lvars & rvars
            lsend = {}
            rsend = {}
            for v in ivars:
                lsend[v] = tp.get_var(v.name + "_l")
                rsend[v] = tp.get_var(v.name + "_r")
            left = left.subs(lsend)
            right = right.subs(rsend)
            lvars = left.all_vars()
            rvars = right.all_vars()
#        print("deduce")
#        print("\t%s"%(left,))
#        print("\t%s"%(right,))
        for i, expr in enumerate(right.src):
            send = Expr.unify(left.tgt, expr)
            if send is None:
                continue
#            print("\tunify %s & %s ==> %s" % (left.tgt, expr, send))
            src = [e.subs(send) for e in left.src + right.src[:i] + right.src[i+1:]]
            tgt = right.tgt.subs(send)
            seq = Sequent(src, tgt, comment, parent=(left, right))
            seq = seq.rename()
            yield seq


yes = Sequent()
always = yes.then


class Theory(object):
    def __init__(self, name):
        self.name = name
        self.vars = {} # _namespace for Variable's
        self.ops = {}
        self.seqs = []
        self.found = set()
        self.add_equality()

    def add(self, seq):
        assert isinstance(seq, Sequent)
        if seq in self.found:
            return
        self.seqs.append(seq)
        self.found.add(seq)

    def __getitem__(self, idx):
        return self.seqs[idx]

    def __len__(self):
        return len(self.seqs)

    def add_equality(tp):
        a, b, c = tp.a, tp.b, tp.c
        #tp.has = tp.add_op("%s.has"%tp.name, 1)
        tp.has = Op(tp, "%s.has"%tp.name, 1, False)
        tp.eq = Op(tp, "==", 2)
        #tp.add( tp.has(a).then(a.eq(a), "(eq-reflexive)") )
        #tp.add( (tp.has(a) & tp.has(b) & a.eq(b)).then(b.eq(a), "(eq-symmetric)") )
        #tp.add( (tp.has(a) & tp.has(b) & tp.has(c) & a.eq(b) & b.eq(c)).then(a.eq(c), "(eq-transitive)") )
        tp.add( yes.then(a.eq(a), "(eq-reflexive)") )
        tp.add( a.eq(b).then(b.eq(a), "(eq-symmetric)") )
        tp.add( (a.eq(b) & b.eq(c)).then(a.eq(c), "(eq-transitive)") )

    def add_extensive(self, op):
        if op.arity==0:
            return # these are already extensive 
        # extensivity of eq
        lhs = yes
        us = [self.get_var(ascii_letters[i]+"0") for i in range(op.arity)]
        vs = [self.get_var(ascii_letters[i]+"1") for i in range(op.arity)]
        for u,v in zip(us, vs):
            #lhs = lhs & self.has(v) & self.has(u) & self.eq(u, v)
            lhs = lhs & self.eq(u, v)
        self.add( lhs.then( self.eq(op(*vs), op(*us)), "(%s-eq-extensive)"%op.name ) )

    def add_op(self, name, arity=0, infix=True, partial=False):
        assert name not in self.ops
        op = Op(self, name, arity, infix)
        self.ops[name] = op
        #if not partial:
        #    vs = [self.get_var(ascii_letters[i]) for i in range(op.arity)]
        #    lhs = yes
        #    for v in vs:
        #        lhs = lhs & self.has(v)
        #    self.add( lhs.then( self.has(op(*vs)), "(%s-endo)"%op.name ) )
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

    def deduce(self, idx, jdx, verbose=False):
        left = self[idx]
        right = self[jdx]
        comment = "[%s] & [%s]"%(idx, jdx)
        return left.deduce(right, comment, verbose=verbose)

    def doprove(self, goal):
        print("goal:", goal)
        for i, seq in enumerate(self):
            send = seq.tgt.specialize(goal.tgt)
            if send is None:
                continue
            print("\tsubgoal:", seq)
            print("\t", send)
            print("\tsubgoal:", seq.subs(send, "from [%d]"%i))
            print()

    def dump(self):
        print("%s:"%self.name)
        for idx, seq in enumerate(self.seqs):
            print("[%d]\t%s"%(idx, seq))


def test():

    # ----------------------------------------

    tp = Theory("T")

    w, x, y, z = tp.w, tp.x, tp.y, tp.z

    assert repr(w) == 'w'

    f = Op(tp, "f", 3, False)
    g = Op(tp, "g", 1, False)

    lhs = f(x, y, g(z))
    rhs = f(x, y, w)

    assert str(lhs) == "f(x, y, g(z))", str(lhs)
    assert lhs.all_vars() == {x, y, z}
    assert not lhs.specialize(rhs)
    assert rhs.specialize(lhs) == {x:x, y:y, w:g(z)}
    assert not lhs.specialize(rhs, iso=True)
    assert not rhs.specialize(lhs, iso=True)

    send = Expr.unify(lhs, rhs)
    assert str(send) == "{w: g(z)}", str(send)

    send = Expr.unify(x, g(x))
    assert send is None

    lhs = f(x, y, z)
    rhs = f(y, x, w)

    send = Expr.unify(lhs, rhs)
    assert send == {y: x, w: z}, send

    assert lhs.specialize(rhs, iso=True) == {x:y, y:x, z:w}

    # ----------------------------------------

    # unify (a_l==(a_l*1)) & (1==a_r) ==> {a_l: 1, a_r: (a_l*1)}
    a_l, a_r = tp.a_l, tp.a_r
    eq = tp.eq
    one = tp.const("1")
    tp.add_op("*", 2)

    lhs = eq(a_l, a_l*one)
    rhs = eq(one, a_r)
    send = Expr.unify(lhs, rhs)
    assert send == {a_l:one, a_r:one*one}

    # ----------------------------------------

    # ----------------------------------------
    # Equality 

    tp = Theory("T")
    #tp.dump()

    for seq in tp:
        other = seq.subs( {tp.a : tp.z} )
        send = seq.weak_eq(other)
        assert send and send[tp.a] == tp.z
        for _seq in tp:
            send = seq.weak_eq(_seq)
            assert (_seq==seq) == (send is not None)


class Prover(object):
    def __init__(self, tp):
        self.tp = tp
    
class Pair(object):
    def __init__(self, left, right):
        self.left = left
        self.right = right
        self.n = left.size + right.size

    def __lt__(self, other):
        return self.n < other.n


def search(tp, tgt, maxsize=1000, verbose=False):
    if verbose:
        print("="*79)
        print("search: tgt = %s" % tgt)
    found = set(tp.seqs)
    table = list(tp.seqs)
    n = len(table)
    pairs = [Pair(left, right) for left in table for right in table]
    heapq.heapify(pairs)
    count = 0
    while 1:
        count += 1
        if count % 1000 == 0:
            print("heap=%d, found=%d" % (len(pairs), len(found)))
        #print("search:", len(found))
        pair = heapq.heappop(pairs)
        a, b = pair.left, pair.right
        for c in a.deduce(b):
            if c in found:
                continue
            if c.src and not c.src_vars():
                # don't need these... ?
                continue 
            if verbose:
                print("search:", c)
            if c.weak_eq(tgt):
                return c
            found.add(c)
            n = len(table)
            table.append(c)
            for i in range(n):
                for (lhs, rhs) in [(table[i], c), (c, table[i])]:
                    if lhs.can_deduce(rhs):
                        heapq.heappush(pairs, Pair(lhs, rhs))
                #heapq.heappush(pairs, Pair(table[i], c))
                #heapq.heappush(pairs, Pair(c, table[i]))
            if len(found) > maxsize:
                #assert 0, "maxsize reached"
                return 


def show_tree(seq, indent=0):
    print("  "*indent + str(seq))
    if seq is None:
        return
    parent = seq.parent
    if parent:
        show_tree(parent[0], indent+1)
        show_tree(parent[1], indent+1)


def test_monoid_0():
    # ----------------------------------------
    # Monoid 

    tp = Theory("M")
    one = tp.const("1")
    tp.add_op("*", 2)
    eq = tp.eq
    has = tp.has

    a, b, c = tp.a, tp.b, tp.c
    tp.add( (has(a) & has(b) & has(c)).then( eq(a*(b*c), (a*b)*c ), "_assoc" ) )
    tp.add( (has(a)).then( eq(a*one, a) ) )
    tp.add( (has(a)).then( eq(one*a, a) ) )

#    _one = tp.const("_1")
#    tp.add( (has(a)).then( eq(a*_one, a) ) )
#    tp.add( (has(a)).then( eq(_one*a, a) ) )

    #seq = tp.deduce(3, 4)
    #seq = tp.deduce(3, 6)

    tp.dump()

#    lhs = has(a).then(eq(a, a*one))
#    rhs = (has(a) & eq(one, a)).then(eq(one, a))
#    for seq in lhs.deduce(rhs):
#        show_tree(seq)
#        print("fail")
#    return

    #tgt = always(eq(one, _one))
    tgt = has(a).then( eq(one*a, a*one) )

    seq = search(tp, has(a).then( has(a*one) ) )
    tp.add(seq)
    seq = search(tp, has(a).then( has(one*a) ) )
    tp.add(seq)

    tgt = has(a).then( eq(one*a, a*one) )
    seq = search(tp, tgt, 100)
    show_tree(seq)
    
#    for i in range(5):
#        idx = randint(0, len(tp)-1)
#        jdx = randint(0, len(tp)-1)
#        print("deduce", idx, jdx)
#        for seq in tp.deduce(idx, jdx):
#            print("==>", seq)
#            print("==>", seq.size)


def test_monoid():
    # ----------------------------------------
    # Monoid 

    tp = Theory("M")
    one = tp.const("1")
    tp.add_op("*", 2)
    eq = tp.eq
    has = tp.has

    a, b, c = tp.a, tp.b, tp.c
    tp.add( always( eq(a*(b*c), (a*b)*c ), "_assoc" ) )
    tp.add( always( eq(a*one, a) ) )
    tp.add( always( eq(one*a, a) ) )

    tp.dump()

    #tgt = always(eq(one, _one))
    #tgt = always( eq(one*a, a*one) )

    tgt = always( eq(one*a, a*one) )
    seq = search(tp, tgt, 100)
    assert seq is not None
    #show_tree(seq)
    tp.add(seq)

    # uniq identity
    _one = tp.const("_1")
    tp.add( always( eq(a*_one, a) ) )
    tp.add( always( eq(_one*a, a) ) )

    tgt = always( eq(_one*a, a*_one) )
    seq = search(tp, tgt, 1000)
    assert seq is not None
    #show_tree(seq)
    tp.add(seq)

    #return

    tgt = always( eq(one, _one) )
    tgt = always( eq(one*_one, _one) )

    seq = search(tp, tgt, 4000, verbose=True)
    show_tree(seq)

    

def test_freyd_scedrov():

    # ----------------------------------------
    # Category, see Freyd-Scedrov 

    tp = Theory("C")
    tp.add_op("*", 2, partial=True)
    tp.add_op("src", 1)
    tp.add_op("tgt", 1)
    eq = tp.eq

    f, g, h = tp.f, tp.g, tp.h
    for seq in [
        tp.has(f).then( tp.has( f.src ) ),
        tp.has(f).then( tp.has( f.tgt ) ),
        (tp.has(f) & tp.has(g) & f.src.eq(g.tgt) ).then( tp.has( f*g ) ),
        (tp.has(f) & tp.has(g) & tp.has(f*g) ).then( f.src.eq(g.tgt) ),
        tp.has(f).then( eq( f.src.tgt, f.src ) ),
        tp.has(f).then( eq( f.tgt.src, f.tgt ) ),
        tp.has(f).then( eq( f.tgt*f, f ) ),
        tp.has(f).then( eq( f*f.src, f ) ),
        (tp.has(f) & tp.has(g) & tp.has(f*g) ).then( eq( (f*g).tgt, (f*g.tgt).tgt )),
        (tp.has(f) & tp.has(g) & tp.has(f*g) ).then( eq( (f*g).src, (f.src*g).src )),
        (tp.has(f) & tp.has(g) & tp.has(h) & tp.has(f*g) & tp.has(g*h) ).then( eq( f*(g*h), (f*g)*h ), "_assoc" ),
    ]:
        tp.add(seq)
    
    tp.dump()

    # yikes...
    goal = tp.has(f).then(eq(f.tgt.tgt, f.tgt))
    #goal = (tp.has(f) & tp.has(g) & tp.has(h) & tp.has(f*g) & tp.has(g*h)).then(tp.has(f*(g*h)))
    tp.doprove(goal)


    


if __name__ == "__main__":

    test()

    from argv import argv
    name = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%name)
    else:
        fn = eval(name)
        fn()

    print("OK\n")

    



