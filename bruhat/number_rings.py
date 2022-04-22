#!/usr/bin/env python

"""
some quadratic number rings...

See:
https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/quadratic-integers-and-coxeter-groups/CF262D475903A0104145D1294DA80EF9
Quadratic Integers and Coxeter Groups
1999 Norman W. Johnson and Asia Ivić Weiss

"""


from functools import total_ordering


from bruhat.element import GenericElement, Ring
from bruhat.lin import Space, Lin





class NumberRing(Ring):
    def __init__(self, reduct=(-1, 0)):
        self.zero = Number((0, 0), self)
        self.one = Number((1, 0), self)
        self.i = Number((0, 1), self)
        self.reduct = reduct

    def promote(self, a):
        if isinstance(a, Number):
            assert a.tp == self
            return a
        if isinstance(a, int):
            return Number((a, 0), self)
        try:
            a = int(a)
            return Number((a, 0), self)
        except ValueError:
            pass
        except TypeError:
            pass
        return None

    def __truediv__(self, other):
        assert isinstance(other, Ideal)
        return QuotientRing(self, other)

    def norm(self, a):
        a = a.value
        return a[0]**2 + a[1]**2

    def add(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]+b[0], a[1]+b[1]), self)

    def sub(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]-b[0], a[1]-b[1]), self)

    def mul(self, a, b):
        a, b = a.value, b.value
        reduct = self.reduct
        c = a[1]*b[1]
        re = a[0]*b[0] + reduct[0]*c
        im = a[0]*b[1] + a[1]*b[0] + reduct[1]*c
        return Number((re, im), self)

    def neg(self, a):
        a = a.value
        return Number((-a[0], -a[1]), self)

    def mod(self, a, b):
        a, b = a.value, b.value
        return Number((a[0]%b[0], a[1]%b[1]), self)

    def i_find(self):
        yield self.zero
        one, i = self.one, self.i
        n = 1
        while 1:
            x = n*one
            while x.value[0] > 0:
                yield x
                x += i-1
            while x.value[1] > 0:
                yield x
                x += -i-1
            while x.value[0] < 0:
                yield x
                x += -i+1
            while x.value[1] < 0:
                yield x
                x += i+1
            n += 1

    def find(self, count):
        items = self.i_find()
        for i in range(count):
            yield items.__next__()


@total_ordering
class Number(GenericElement):

    def __str__(self):
        a, b = self.value
        if b==0:
            s = str(a)
        elif a==0:
            s = "%d*i"%b
        else:
            s = "(%d+%d*i)"%(a, b)
            s = s.replace("+-", "-")
        return s

    # choose some total_ordering ...
    def __lt__(a, b):
        assert isinstance(b, Number)
        assert a.tp == b.tp
        tp = a.tp
        ra = tp.norm(a)
        rb = tp.norm(b)
        if ra < rb:
            return True
        if ra > rb:
            return False
        # now ra==rb
        a, b = a.value, b.value
        if a[0] > b[0]:
            return True
        if a[0] < b[0]:
            return False
        # now a[0] == b[0]
        if a[1] > b[1]:
            return True
        return False


class Ideal(object):
    "_principle ideal"
    def __init__(self, ring, a):
        assert isinstance(ring, NumberRing)
        assert a.tp == ring
        self.ring = ring
        #assert a.value[0]>0 and a.value[1]>0
        self.a = a
        self.cache = {}
        i = self.i_find()
        items = []
        while len(items) < 10:
            x = i.__next__()
            if x != 0:
                items.append(x)
        self.items = items

    def __contains__(self, b):
        pass

    def i_find(self):
        ring = self.ring
        a = self.a
        for x in ring.i_find():
            yield x*a

    def reduce(self, src):
        ring = self.ring
        cache = self.cache
        if src in cache:
            return cache[src]
        items = self.items
        done = False
        tgt = src
        while not done:
            done = True
            for x in items:
                if tgt + x < tgt:
                    tgt = tgt+x
                    done = False
        cache[src] = tgt
        return tgt


class QuotientRing(Ring):
    def __init__(self, ring, ideal):
        self.ring = ring
        self.ideal = ideal
        self.one = GenericElement(ring.one, self) # wrapped
        self.zero = GenericElement(ring.zero, self) # wrapped
        self.i = GenericElement(ring.i, self) # wrapped

        # let's hope it's finite and not too big.
        self.items = self._get_items()

    def _get_items(self, N=5):
        one, i = self.one, self.i
        items = set()
        for i in range(-N, N):
          for j in range(-N, N):
            items.add(i*one + j*i)
        items = list(items)
        items.sort(key=str)
        return items

    def __getitem__(self, i):
        return self.items[i]

    def promote(self, a):
        ring = self.ring
        if isinstance(a, GenericElement):
            assert a.tp == self
            return a
        if isinstance(a, Number):
            assert a.tp == ring
            return GenericElement(a, self)
        if isinstance(a, int):
            a = ring.promote(a)
            return GenericElement(a, self)
        try:
            a = int(a)
            a = ring.promote(a)
            return GenericElement(a, self)
        except ValueError:
            pass
        except TypeError:
            pass
        return None

    def add(self, a, b):
        ring, ideal = self.ring, self.ideal
        value = ring.add(a.value, b.value)
        value = ideal.reduce(value)
        a = GenericElement(value, self)
        return a

    def sub(self, a, b):
        ring, ideal = self.ring, self.ideal
        value = ring.sub(a.value, b.value)
        value = ideal.reduce(value)
        a = GenericElement(value, self)
        return a

    def mul(self, a, b):
        ring, ideal = self.ring, self.ideal
        value = ring.mul(a.value, b.value)
        value = ideal.reduce(value)
        a = GenericElement(value, self)
        return a

    def get_gl(self, V):
        assert V.n == 2, "todo.."
        lins = []
        for a in self:
         for b in self:
          for c in self:
           for d in self:
               if a*d-b*c != 0:
                   f = Lin(V, V, [[a, b], [c, d]])
                   lins.append(f)
        return lins

    def get_sl(self, V):
        assert V.n == 2, "todo.."
        lins = []
        for a in self:
         for b in self:
          for c in self:
           for d in self:
               if a*d-b*c == 1:
                   f = Lin(V, V, [[a, b], [c, d]])
                   lins.append(f)
        return lins


# does not need hashable operators
def mulclose(gen, verbose=False, maxsize=None):
    ops = list(gen)
    bdy = gen
    while bdy:
        _bdy = []
        for g in bdy:
            for h in gen:
                k = g*h
                for j in ops:
                    if j==k:
                        break
                else:
                    ops.append(k)
                    _bdy.append(k)
        bdy = _bdy
        if verbose:
            print("mulclose:", len(ops))
        if maxsize and len(ops) >= maxsize:
            break
    return ops



def main():

    # ----------------------------
    # Gaussian integers
    R = NumberRing((-1, 0))
    one = R.one
    i = R.i

    assert i*i == -one
    assert i**3 == -i

    assert (1+i)**2 == 2*i
    assert (2+i)*(2-i) == 5
    assert (3+2*i)*(3-2*i) == 13

    items = []
    i = R.i_find()
    while len(items) < 20:
        items.append(i.__next__())
    items.sort()
    #print(' '.join(str(x) for x in items))

    # ----------------------------
    # Eisenstein integers
    R = NumberRing((-1, -1))
    one = R.one
    i = R.i

    assert i**3 == 1
    assert 1 + i + i**2 == 0

    # ----------------------------
    # Kleinian integers
    R = NumberRing((-2, -1))
    one = R.one
    i = R.i

    assert i * (-1-i) == 2
    assert i * (1+i) == -2

    # ----------------------------

    R = NumberRing((-1, 0))
    I = Ideal(R, 2-R.i)
    S = R/I

    zero, one, i = S.zero, S.one, S.i
    
    assert 1+i != 0
    assert 2-i == 0
    assert 2+i != 0
    assert 1-i == 2+i

    V = Space(S, 2)
    SL = S.get_sl(V)
    GL = S.get_gl(V)

    def lin(a, b, c, d):
        a = S.promote(a)
        b = S.promote(b)
        c = S.promote(c)
        d = S.promote(d)
        return Lin(V, V, [[a, b], [c, d]])
        
    # generators for SL
    gen = lin(0, 1, -1, 0), lin(1, 0, 1, 1), lin(1, 0, i, 1)
    G = mulclose(gen)
    assert len(G) == 120



if __name__ == "__main__":

    main()










