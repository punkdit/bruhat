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
    

class Tom:
    def __init__(self, rows):
        self.rows = rows

    def str(self):
        lines = []
        for item in self.rows:
            line = " ".join((str(i) if i else ".") for i in item)
            lines.append(line)
        return '\n'.join(lines)

    def __getitem__(self, key):
        i = letters.index(key)
        return self.rows[i]

    def __len__(self):
        return len(self.rows)

    def __str__(self):
        rows = self.rows
        N = len(rows)
        smap = SMap()
        for j in range(N):
            smap[0,3*j] = "%3s"%letters[j]
            smap[j+1,3*N+2] = letters[j]
        for i,row in enumerate(rows):
            for j,c in enumerate(row):
                if j<i:
                    continue
                smap[i+1,3*j] = "%3s"%(c if c else ".")
        return str(smap)

    def get_name(self, items):
        rows = self.rows
        N = len(rows)
        assert len(items) == N
        items = items.copy()
        name = []
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
            s = "%d*%s"%(j, letters[i]) if j>1 else letters[i]
            name.append(s)
        return "+".join(name)


def test():

    tom = load_tom("A5")
    B = tom["B"]
    C = tom["C"]

    BC = B*C
    #assert tom.get_name(BC) == "H"

    BG = tom["B"]*tom["G"]
    assert tom.get_name(BG) == "2*G+I"

    #return

def main():

    tom = load_tom("M11")

    N = len(tom)
    names = letters[:N]
    for a in names:
      for b in names:
        A = tom[a]
        B = tom[b]
        AB = A*B
        #print(a,"*", b)
        #print(AB)
        name = tom.get_name(AB)
        print("%6s"%name, end=" ")
      print()






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

    

    
