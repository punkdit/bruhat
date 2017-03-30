#!/usr/bin/env python
"""
Find idempotent elements of Burnside ring, apart from 0 and 1.
These exist iff the group is non-solvable.
"""

import sys, os


# parse this multiplication table:
data = """\
  | A  B      C        D        E        F        G        H         I     
--+------------------------------------------------------------------------
A | A  B      C        D        E        F        G        H         I     
B | B  B+G    H        G+H      I        F+I      2*G+I    H+2*I     5*I   
C | C  H      C+H      2*H      E+I      3*H      2*I      2*H+2*I   6*I   
D | D  G+H    2*H      D+H+I    2*I      3*H+I    G+3*I    2*H+4*I   10*I  
E | E  I      E+I      2*I      2*E+2*I  3*I      4*I      6*I       12*I  
F | F  F+I    3*H      3*H+I    3*I      3*F+3*I  5*I      3*H+6*I   15*I  
G | G  2*G+I  2*I      G+3*I    4*I      5*I      2*G+6*I  10*I      20*I  
H | H  H+2*I  2*H+2*I  2*H+4*I  6*I      3*H+6*I  10*I     2*H+14*I  30*I  
I | I  5*I    6*I      10*I     12*I     15*I     20*I     30*I      60*I  """


class Element(object):

    keys = None
    table = None
    def __init__(self, cs=None):
        if cs is None:
            cs = [0]*len(self.keys)
        assert len(cs) == len(self.keys)
        self.cs = cs

    def __str__(self):
        cs = self.cs
        items = []
        for key, c in zip(self.keys, self.cs):
            if c==0:
                continue
            if c==1:
                s = key
            else:
                s = "%s*%s"%(c, key)
            items.append(s)
        s = "+".join(items) or "0"
        #return s
        s = s.replace("+-", "-")
        return s
    __repr__ = __str__

    def __hash__(self):
        return hash(self.__str__())

    def norm(self):
        return sum(abs(c) for c in self.cs)

    def __eq__(self, other):
        return self.cs == other.cs

    def __ne__(self, other):
        return self.cs != other.cs

    def __add__(self, other):
        cs = [a+b for (a,b) in zip(self.cs, other.cs)]
        return Element(cs)

    def __sub__(self, other):
        cs = [a-b for (a,b) in zip(self.cs, other.cs)]
        return Element(cs)

    def __mul__(self, other):
        table = self.table
        a = Element() # zero
        for i, c in enumerate(self.cs):
          for j, d in enumerate(other.cs):
            if c*d==0:
                continue
            product = [c*d*val for val in table[i, j]]
            a = a + Element(product)
        return a

    def __rmul__(self, r):
        return Element([r*c for c in self.cs])

    def __neg__(self):
        return Element([-1*c for c in self.cs])

    def __pow__(self, n):
        cs = [0]*len(self.keys)
        cs[0] = 1
        a = Element(cs)
        for i in range(n):
            a = self * a
        return a

    @classmethod
    def parse(cls, data):
    
        rows = data.split("\n")
    
        header = rows.pop(0)
        rows.pop(0)
    
        keys = header[3:].split()
        print keys
        cls.keys = keys
        n = len(keys)
    
        table = {}
        for row in rows:
            left = row[0]
            for j, d in enumerate(row):
                if j<=2:
                    continue
                right = header[j]
                if right == ' ':
                    continue
                value = row[j:].split()[0]
                vals = value.split("+")
                #print vals
                product = [0 for key in keys]
                for val in vals:
                    if '*' not in val:
                        val = '1*'+val
                    c, val = val.split("*")
                    c = int(c)
                    product[keys.index(val)] = c
                #print product
                #print "%s*%s = %s" % (left, right, product)
                table[keys.index(left), keys.index(right)] = product
        cls.table = table

        basis = []
        for i in range(n):
            cs = [0]*n
            cs[i] = 1
            element = Element(cs)
            basis.append(element)
        return basis


items = Element.parse(data)

if 0:
    A, B, C, D, E, F, G, H, I = items
    print items
    #print (B+C+D-G-2*H+I)**2
    print B*B+2*B*(C+D-G-2*H+I)
    print C*C+2*C*(D-G-2*H+I)
    print D*D+2*D*(-G-2*H+I)
    print G*G+2*(-G)*(-2*H+I)
    print 4*H*H+2*(-2*H)*(I)
    print I*I
    sys.exit(0)


zero = Element()
cs = [0]*len(Element.keys)
cs[0] = 1
one = Element(cs)

#items = set(items+[-1*a for a in items])
items = set(items)

while 1:

    new = set()


    for a0 in items:
      for a1 in items:

        for x in [a0+a1, a0-a1]:
    
            if x==zero or x==one:
                continue
    
            if x*x == x:
                print x
                y = one-x
                print y*y == y
                sys.exit(0)
    
            if x.norm() < 5 and x not in items:
                new.add(x)

    if not new:
        break

    items.update(new)
    print len(items)

    
    
