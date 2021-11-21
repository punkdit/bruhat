#!/usr/bin/env python3

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



class Cayley(object):
    def __init__(self, table):
        self.table = dict(table)

    def __len__(self):
        return len(self.table)

    def dump(self):
        table = self.table
        vals = list(set(table.values()))
        vals.sort()
        print('     '+''.join("%3d"%i for i in vals))
        print('     '+'-'*3*len(vals))
        for i in vals:
          print("%3d |"%i, end="")
          for j in vals:
            k = table.get((i, j))
            if k is None:
                s = "   "
            else:
                s = "%3d"%k
            print(s, end="")
          print()

    def _deduce(self):
        table = self.table
        vals = list(set(table.values()))
        vals.sort()
        for i in vals:
            found = {}
            for j in vals:
                k = table.get((i, j))
                if k is None:
                    continue
                if k not in found:
                    found[k] = j
                    continue
                j1 = found[k]
                # i*j == i*j1 == k => j==j1
                self.rewrite(j1, j)
                return True
        return False

    def deduce(self):
        while self._deduce():
            pass

    def rewrite(self, i, j):
        # send j to i
        table = self.table
        keys = list(table.keys())
        for ii, jj in keys:
            kk = table[ii, jj]
            changed = False
            if ii == j or jj == j:
                del table[ii, jj]
                if ii == j:
                    ii = i
                if jj == j:
                    jj = i
                changed = True
            if kk == j:
                kk = i
                changed = True
            if changed:
                table[ii, jj] = kk



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

    lookup = dict((g, idx) for (idx, g) in enumerate(G))

    table = {}
    for g in G:
      for h in G:
        k = g*h
        if k not in lookup:
            continue
        table[lookup[g], lookup[h]] = lookup[k]

    table = Cayley(table)
    ops = list(G)
    i = 0
    while i < len(ops):
        g = ops[i]
        j = i+1
        while j < len(ops):
            h = ops[j]
            if (~g)*h in H:
                table.rewrite(lookup[g], lookup[h])
                ops.pop(j)
            else:
                #print((~g)*h)
                j += 1
        i += 1
        print(len(table), len(ops))

    for op in ops:
        print(op)

    table.dump()
    table.deduce()
    print()
    table.dump()
    for i in table.table.values():
        print("%d : %s" % (i, G[i]))

    return


    print("OK")




if __name__ == "__main__":

    test()


