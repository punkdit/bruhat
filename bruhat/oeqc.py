#!/usr/bin/env python

from time import time
start_time = time()

from bruhat.argv import argv

from bruhat.solve import parse, span, shortstr, array2
from bruhat.isomorph import Tanner, search
from bruhat.equ import Equ, quotient


def preproc():
    data = open(argv.next()).read()
    
    data = data.replace(".", "")
    data = data.replace("3", "")
    data = data.replace(" ", "")
    data = data.replace("\n\n", "\n")
    data = data.replace("[[", "[")
    data = data.replace("]]", "]\n")
    data = data.replace("[", "")
    data = data.replace("]", "")
    
    #data = data.replace("\n\n", ".")
    #data = data.replace("\n", "")
    #data = data.replace(".", "\n")
    
    f = open('oeqc.txt', 'w')
    f.write(data)
    f.close()


def find_isos(H1, H2):
    src = Tanner.build2(H1, H1)
    tgt = Tanner.build2(H2, H2)
    G = [f for f in search(src, tgt)]
    return G

def is_iso(H1, H2):
    src = Tanner.build2(H1, H1)
    tgt = Tanner.build2(H2, H2)
    for f in search(src, tgt):
        return True
    return False


def read_codes():
    print("read_codes")
    data = open("oeqc.txt").read()
    items = data.split("\n\n")
    #print(len(items))
    
    found = set()
    for item in items:
        item = item.strip()
        if not item or item[0] not in "01":
            #print("skipping:", repr(item))
            continue
        try:
            A = parse(item)
        except ValueError:
            print("parse failed:", item)
            raise
        yield A

def get_dist():
    for A in read_codes():
        n = A.shape[1]
        dist = [0 for i in range(n)]
        for v in span(A):
            dist[v.sum()] += 1
        dist = tuple(dist)
        if dist not in found:
            print(dist)
            found.add(dist)

def get_autos():
    H = parse("""
    101001010011
    011010010011
    000111010011
    000000110101
    000000001111
    """)

    rows = [v for v in span(H) if v.sum()]
    V = array2(rows)
    print(shortstr(V))
    count = 0
    for f in find_isos(V, V):
        print(f)
        count += 1
    print(count)


        
def find_iso():
    codes = list(read_codes())
    print(len(codes))

    #codes = codes[:100]
    #codes = [H.tobytes() for H in codes]
    equs = quotient(codes, is_iso, verbose=True)
    found = set()
    for equ in equs:
        found.add(equ.top)
    print()
    print("distinct codes:", len(found))
    for equ in found:
        print(shortstr(equ.items[0]))
        print()
    return


    n = len(codes)
    H1 = codes[0]
    for H2 in codes[1:]:
        if not is_iso(H1, H2):
            print("distinct")
            break
        else:
            print(".", flush=True, end="")
    print()


if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "test_all"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()


