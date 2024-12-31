#!/usr/bin/env python

import sys
from os import popen
from string import ascii_lowercase, ascii_uppercase
letters = ascii_uppercase+ascii_lowercase

import numpy
numpy.set_printoptions(threshold=sys.maxsize)

from bruhat.argv import argv
from bruhat.smap import SMap


def load_tom(name="M11"):
    cmd = """
    SizeScreen([20000,20000]);
    tom := TableOfMarks("%s");
    Display(tom);
    quit;
    """%name
    
    f = open("/tmp/tmp.gap", "w")
    print(cmd, file=f)
    f.close()
    
    lines = popen("gap --norepl /tmp/tmp.gap").readlines()
    
    tom = []
    start = False
    for line in lines:
        line = line.strip()
        if line.startswith("1:"):
            start = True
        if not start:
            continue
        if not line:
            break
        flds = line.split(":")[1]
        flds = flds.replace(".", "0")
        flds = flds.split()
        flds = [int(i) for i in flds]
        tom.append(flds)
    
    #print([len(item) for item in tom], len(tom))
    assert len(tom[-1]) == len(tom)
    
    N = len(tom)
    for item in tom:
        item += [0]*(N-len(item))
    
    lrev = lambda items:list(reversed(items))
    
    tom = [lrev(item) for item in reversed(tom)]
    tom = numpy.array(tom)
    return Tom(tom)
    

def get_letter(i, post=0):
    if i >= len(letters):
        i -= len(letters)
        post += 1
        return get_letter(i, post)
    else:
        return letters[i]+("_"+str(post) if post else "")


class Tom:
    def __init__(self, rows, names=None):
        self.rows = rows
        N = len(rows)
        if names is None:
            names = [get_letter(i) for i in range(N)]
        self.names = list(names)

    def str(self):
        lines = []
        for item in self.rows:
            line = " ".join((str(i) if i else ".") for i in item)
            lines.append(line)
        return '\n'.join(lines)

    def __getitem__(self, key):
        i = self.names.index(key)
        return self.rows[i]

    def __len__(self):
        return len(self.rows)

    def __str__(self):
        rows = self.rows
        names = self.names
        N = len(rows)
        smap = SMap()
        for j in range(N):
            smap[0,3*j] = "%3s"%names[j]
            smap[j+1,3*N+2] = names[j]
        for i,row in enumerate(rows):
            for j,c in enumerate(row):
                if j<i:
                    continue
                smap[i+1,3*j] = "%3s"%(c if c else ".")
        return str(smap)

    def get_desc(self, items):
        rows = self.rows
        names = self.names
        N = len(rows)
        assert len(items) == N
        items = items.copy()
        desc = []
        for i in range(N):
            assert items.min() >= 0
            if items[i] == 0:
                continue
            #print(items, "i=", i)
            n = items[i]
            row = rows[i]
            assert n % row[i] == 0
            j = n//row[i]
            items -= j*row
            s = "%d*%s"%(j, names[i]) if j>1 else names[i]
            desc.append(s)
        return "+".join(desc)


    def dump_burnside(self):
        N = len(self)
        names = self.names
        for a in names:
          for b in names:
            A = self[a]
            B = self[b]
            AB = A*B
            #print(a,"*", b)
            #print(AB)
            name = self.get_desc(AB)
            print("%6s"%name, end=" ")
          print()


def test():

    tom = load_tom("A5")
    B = tom["B"]
    C = tom["C"]

    BC = B*C
    #assert tom.get_desc(BC) == "H"

    BG = tom["B"]*tom["G"]
    assert tom.get_desc(BG) == "2*G+I"

    tom.dump_burnside()
    #return


def main():

    #test()
    #tom = load_tom("A5")
    #tom.dump_burnside()

    tom = load_tom("M24")
    N = len(tom)

    rows = tom.rows

    print("maximal subgroups / minimal structures:")
    for i,row in enumerate(rows):
        col = rows[:,i]
        found = [j for j in col if j]
        if len(found)==2:
            print(tom.names[i], row[N-1])
    
    octad = tom.rows[4]
    desc = tom.get_desc(octad*octad)
    print("octad*octad =", desc)

    names = desc.split("+")
    div = octad[N-1]
    for name in names:
        row = tom[name]
        count = row[N-1]
        print(name, count//div)
        assert count%div == 0
    print()

    octern = tom["s_8"]

    desc = tom.get_desc(octern*octern)
    print(desc)

    names = desc.split("+")
    div = octern[N-1]
    total = 0
    for name in names:
        mul = 1
        if "*" in name:
            mul, name = name.split("*")
            mul = int(mul)
        row = tom[name]
        count = mul*row[N-1]
        print(name, count//div)
        total += count // div
        assert count%div == 0
    print()
    print(total)







if __name__ == "__main__":
    from time import time
    start_time = time()

    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))

    

    
