#!/usr/bin/env python3

"""
Failed attempt at working with group presentations.
See todd_coxeter.py for a working version of this.
"""

from bruhat.action import mulclose_fast as mulclose

print


class Free(object):

    def __init__(self, *items):
        self.items = items

    def __mul__(self, other):
        left, right = list(self.items), list(other.items)
        while left and right:
            l, r = left[-1], right[0]
            if "~"+l==r or l=="~"+r:
                left.pop(-1)
                right.pop(0)
            else:
                break
        items = left + right
        return Free(*items)

    def __pow__(g, e):
        if e < 0:
            return (~g).__pow__(-e) # recurse
        op = Free()
        for i in range(e):
            op = g*op
        return op

    def __invert__(self):
        items = []
        for item in reversed(self.items):
            if item[0] == "~":
                items.append(item[1:])
            else:
                items.append("~"+item)
        return Free(*items)

    def __eq__(self, other):
        return self.items == other.items

    def __lt__(self, other):
        return self.items < other.items

    def __hash__(self):
        return hash(self.items)

    def __str__(self):
        return "*".join(self.items) or "1"
    __repr__ = __str__





def test():

    I = Free()
    a, b, c = Free("a"), Free("b"), Free("c")

    assert I*a == a
    assert a*b != b*a
    assert ~(a*b) == (~b)*(~a)
    abc = a*b*c
    assert abc * ~abc == I
    assert str(abc) == "a*b*c"

    G = mulclose([a, b, ~a, ~b], maxsize=40)
    G.sort(key = lambda g : str(g))

    gen = [a*a, b*b, (a*b)**3]
    gen = gen + [~g for g in gen]
    H = mulclose(gen, maxsize=200)
    H = H + [g * h * (~g) for g in G for h in H]
    H = set(H)
    #H = set(h*g for g in H for h in H)

    print("OK")




if __name__ == "__main__":

    test()


