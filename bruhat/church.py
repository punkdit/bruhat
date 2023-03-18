#!/usr/bin/env python

" lambda calculus "

import string, os
from time import sleep, time
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)


from bruhat.argv import argv


class Term(object):
    def __mul__(self, other):
        return Apply(self, other)
    def __call__(self, *args):
        if args:
            head, tail = args[0], args[1:]
            term = Apply(self, head)
            return term(*tail) # recurse
        term = self.beta()
        return term


class Var(Term):
    counter = 0
    def __init__(self, name=None):
        if name is None:
            name = "v%s"%Var.counter
            Var.counter += 1
        self.name = name
    def __str__(self):
        return self.name
    def __eq__(self, other):
        return isinstance(other, Var) and self.name==other.name
    def free(self):
        yield self
    def subs(self, u, t): # replace u with t
        assert isinstance(u, Var)
        return t if self==u else self
    def beta(self): # beta reduction
        return self


class Lambda(Term):
    def __init__(self, var, term):
        assert isinstance(var, Var)
        assert isinstance(term, Term)
        self.var = var
        self.term = term
    def __str__(self):
        return "(^%s.%s)"%(self.var, self.term)
    def __eq__(self, other):
        if not isinstance(other, Lambda):
            return False
        if self.var == other.var:
            return self.term == other.term
        #return (self.var, self.term) == (other.var, other.term)
        return self.term == other.term.subs(other.var, self.var)
    def free(self):
        var, term = self.var, self.term
        for v in term.free():
            if v!=var:
                yield var
    def subs(self, u, t): # replace u with t
        assert isinstance(u, Var)
        var, term = self.var, self.term
        if u==var:
            return self
        if var in list(t.free()):
            # alpha convert
            var1 = Var()
            self = Lambda(var1, term.subs(var, var1))
            return self.subs(u, t) # recurse
        return Lambda(var, term.subs(u, t))
    def beta(self): # beta reduction
        var, term = self.var, self.term
        term = term.beta()
        return Lambda(var, term)


class Apply(Term):
    def __init__(self, lhs, rhs):
        assert isinstance(lhs, Term)
        assert isinstance(rhs, Term)
        self.lhs = lhs
        self.rhs = rhs
    def __str__(self):
        return "(%s %s)"%(self.lhs, self.rhs)
    def __eq__(self, other):
        if not isinstance(other, Apply):
            return False
        return (self.lhs, self.rhs) == (other.lhs, other.rhs)
    def free(self):
        for v in self.lhs.free():
            yield v
        for v in self.rhs.free():
            yield v
    def subs(self, u, t): # replace u with t
        assert isinstance(u, Var)
        lhs, rhs = self.lhs, self.rhs
        lhs = lhs.subs(u, t)
        rhs = rhs.subs(u, t)
        return Apply(lhs, rhs)
    def beta(self): # beta reduction
        lhs, rhs = self.lhs, self.rhs
        lhs = lhs.beta()
        rhs = rhs.beta()
        if isinstance(lhs, Lambda):
            term = lhs.term.subs(lhs.var, rhs)
            term = term.beta() # << recurse
        else:
            term = Apply(lhs, rhs)
        return term

def eq(lhs, rhs):
    lhs = lhs.beta()
    rhs = rhs.beta()
    return lhs == rhs


def test():

    a, b, c, f, m, n, x, y, z = [Var(c) for c in 'abcfmnxyz']

    fxx = Lambda(x, x)
    fyy = Lambda(y, y)
    fxy = Lambda(x, y)
    fyx = Lambda(y, x)
    assert str(fxx) == '(^x.x)'
    assert fxx == fyy
    assert fxy != fxx
    assert fxy != fyx

    xy = Apply(x, y)

    assert list(fxx.free()) == []
    assert list(xy.free()) == [x, y]

    assert Apply(fxx, x).beta() == x
    assert Apply(fyy, x).beta() == x
    assert Lambda(x, y) == Lambda(z, y)
    assert fxy.subs(y, x) == Lambda(z, x)

    nums = [
        Lambda(f, fxx),
        Lambda(f, Lambda(x, f*x)),
        Lambda(f, Lambda(x, f*(f*x))),
        Lambda(f, Lambda(x, f*(f*(f*x)))),
        Lambda(f, Lambda(x, f*(f*(f*(f*x))))),
    ]

    L = Lambda
    succ = L(n, L(f, L(x, (f * (n*f*x)))))
    plus = L(m, L(n, L(f, L(x, (m*f)*(n*f*x)))))

    assert succ(nums[0]) == nums[1] 
    assert succ(nums[1]) == nums[2] 
    assert succ(nums[3]) == nums[4] 

    assert plus(nums[1], nums[1]) == nums[2]
    assert plus(nums[2], nums[1]) == nums[3]

    for i in range(10):
        nums.append(succ(nums[-1]))

    #mul = L(m, L(n, m*(plus * n) *nums[0])) # fail
    mul = L(m, L(n, L(f, m*(n*f))))
    assert mul(nums[2], nums[3]) == nums[6]
    assert mul(nums[3], nums[3]) == nums[9]

    # fail
    power = L(x, L(y, y*x))
    #print(power(nums[2], nums[2]) == nums[4])
    #print(power(nums[2], nums[3]))
    #print(nums[8])
    #assert power(nums[2], nums[3]) == nums[8]

    I = L(x, x)
    K = L(x, L(y, x))
    S = L(x, L(y, L(z, ((x*z)*(y*z)))))
    SII = S*I*I

    assert SII(a) == a*a
    #assert SII(SII) == SII*SII # infinite recursion

    assert (S(K(a))(SII))(b) == a(b*b)

    assert (S(K(S*I))*K)(a, b) == b*a

    assert (S*K*I)(K*I*S) == I
    assert (K*S)(I(S*K*S*I)) == S
    assert (S*K*I*K)() == K
    


if __name__ == "__main__":

    start_time = time()


    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")



