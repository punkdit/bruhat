#!/usr/bin/env python

"""
abelian galois extensions / cyclotomic number fields 
"""

from time import time
start_time = time()
import string

from sage.all_cmdline import *
from sage import all_cmdline 

from bruhat.argv import argv
from bruhat.smap import SMap
from huygens.namespace import *

yellow = RGB(1.0, 0.8, 0.2)
orange = RGB(252/255, 132/255, 3/255) * (1./0.8)


def main():
    n = argv.get("n")
    if n is not None:
        rows = argv.get("rows")
        cols = argv.get("cols")
        build(n, rows, cols)

    else:
        build(7, 1, 6)
        return
        build(17, 1, 16)
        build(20, 2, 4)
        build(44, 2, 10)
        build(63, 6, 6)

class Table(object):
    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self.values = {}
    def __setitem__(self, key, value):
        i, j = key
        self.values[key] = value
    def __getitem__(self, key):
        return self.values[key]
    @property
    def smap(self):
        s = SMap()
        rows, cols = self.rows, self.cols
        for row in range(rows):
          for col in range(cols):
            value = self.values.get((row, col), ".")
            s[row, col] = value
        return s
    def __str__(self):
        s = self.smap
        return str(s)
    def render(self):
        cvs = Canvas()
        dx, dy = 1., 1.
        m, n = (self.rows), (self.cols)
        values = self.values
        for i in range(m):
         for j in range(n):
            k = values.get((i, j))
            cvs.fill(path.rect(j*dx, -i*dy, dx, -dy), [white if k is None else orange])
        st = st_THick
        for i in range(m+1):
            cvs.stroke(path.line(0, -i*dy, n*dx, -i*dy), st)
        for i in range(n+1):
            cvs.stroke(path.line(i*dx, 0, i*dx, -m*dy), st)
        return cvs



def build(nn, nrows, ncols):
    K = CyclotomicField(nn)
    G = K.galois_group()
    Hs = G.subgroups()
    for H in Hs:
        if len(H) != 18:
            continue
        K0 = H.fixed_field()[0]
        #print(K0)
        #print(K0.defining_polynomial())

    items = []
    for H in Hs:
        gens = H.gens()
        table = Table(nrows, ncols)
        for g in H:
          cs = list(g.cycles()) or [(1,)]
          c = cs[0]
          c = eval(str(c)) # wtf
          for ii, i in enumerate(c):
              i -= 1 # one-based index
              col = i%ncols
              row = i//ncols
              if ii in [0,1]:
                table[row, col] = "X"
        print(table)
        K0 = H.fixed_field()[0]
        print(K0.defining_polynomial())
        print(H.fixed_field())
        print()
        item = (H, table, K0.defining_polynomial())
        items.append(item)

    N = 1
    col = 0
    smap = SMap()
    rows = []
    row = []
    for item in items:
        if len(item[0]) == N:
            smap.paste(item[1].smap, 0, col)
        else:
            rows.append(row)
            row = []
            print(smap, nrows*ncols//N, '\n')
            N = len(item[0])
            smap = SMap()
            col = 0
            smap.paste(item[1].smap, 0, col)
        row.append(item)
        col += ncols+1
    rows.append(row)
    print(smap, nrows*ncols//N, '\n')

    fg = Canvas()
    x = 0
    y = 0
    margin = 0.2
    get_width = lambda bb : 1.2*bb.width + margin
    lookup = {} # map subgroup to cvs
    for row in rows:
        for item in row:
            lookup[item[0]] = item[1].render()
    rows = [[lookup[item[0]] for item in row] for row in rows]
    widths = [sum(get_width(cvs.get_bound_box()) for cvs in row) for row in rows]
    width = max(widths)
    for row in rows:
        x = 0.5 * width * (1/len(row))
        for cvs in row:
            bb = cvs.get_bound_box()
            fg.insert(x, y, cvs)
            #x += get_width(bb)
            x += width * (1/len(row))
        y -= 2.2*bb.height + margin

    bg = Canvas()
    links = get_links(items)
    st = st_THICk + [grey]
    margin = 0.2
    for (i, j) in links:
        Hi, Hj = Hs[i], Hs[j]
        x0, y0 = fg.find_bound_box(lookup[Hi]).south
        x1, y1 = fg.find_bound_box(lookup[Hj]).north
        #bg.stroke(path.line(x0, y0, x1, y1), st)
        y2 = 0.5*(y0+y1)
        bg.stroke(path.curve(
            x0, y0-margin, 
            x0, y2, 
            x1, y2,
            x1, y1+margin,
        ), st)
    cvs = Canvas([bg, fg])
    cvs.writePDFfile("galois_%s.pdf"%nn)



def get_links(items):
    n = len(items)
    links = set()
    for i in range(n):
     for j in range(i+1, n):
        Hi = items[i][0]
        Hj = items[j][0]
        if Hi.is_subgroup(Hj):
            links.add((i, j))
    for i in range(n):
     for j in range(i+1, n):
      for k in range(j+1, n):
        if (i,j) in links and (j,k) in links and (i,k) in links:
            links.remove((i,k))
    return links

def make_dot(items):
    # not great....
    n = len(items)
    links = get_links(items)
    letters = string.ascii_letters
    lines = ["digraph {"]
    for (i,j) in links:
        #lines.append("  %s -> %s;" % (letters[i], letters[j]))
        lhs, rhs = str(items[i][1]), str(items[j][1])
        lhs, rhs = repr(lhs), repr(rhs)
        lhs = lhs.replace("'", '"')
        rhs = rhs.replace("'", '"')
        lhs = lhs.replace("X", '1')
        rhs = rhs.replace("X", '1')
        lines.append("  %s -> %s;" % (lhs, rhs))
    
    lines.append("}")
    f = open("galois.dot", "w")
    f.write('\n'.join(lines))
    f.close()





if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


