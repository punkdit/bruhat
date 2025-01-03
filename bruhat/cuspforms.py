#!/usr/bin/env python

from os import popen

#from sage.all_cmdline import GL, SL, matrix, Matrix, GF
from sage import all_cmdline as sage

from bruhat.smap import SMap
from bruhat.argv import argv
from bruhat.equ import Equ
from bruhat.algebraic import Algebraic, mulclose, Matrix
from bruhat import algebraic

parse = lambda s : Matrix(algebraic.parse(s.replace(" ", "\n")))


def get_inv(G):
    inv = {}
    I = G.I
    print("inv:", end="", flush=True)
    remain = set(G)
    itemss = []
    while remain:
        #print("|G| =", len(G))
        g = iter(remain).__next__()
        for h in remain:
            if g*h == I:
                inv[g] = h
                inv[h] = g
                break
        else:
            assert 0
        remain.remove(g)
        if g != h:
            remain.remove(h)
        if len(remain)%100 == 0:
            print("[%d]"%len(remain), end="", flush=True)
        
    print()
    assert len(inv) == len(G)
    return inv

def slow_get_conjugacy_classes(G):
    #inv = {g:g.inverse() for g in G}
    inv = get_inv(G)
    remain = set(G)
    itemss = []
    while remain:
        g = iter(remain).__next__()
        remain.remove(g)
        items = [g]
        for h in G:
            k = inv[h]*g*h
            if k in remain:
                remain.remove(k)
                items.append(k)
        itemss.append(items)
    return itemss

def get_cls(G, g):
    gen = G.gen
    found = {g}
    bdy = list(found)
    inv = {h:h.inverse() for h in G.gen}
    while bdy:
        _bdy = []
        for a in gen:
            for g in bdy:
                h = inv[a]*g*a
                if h not in found:
                    found.add(h)
                    _bdy.append(h)
        bdy = _bdy
    return found


def get_conjugacy_classes(G):
    remain = set(G)
    itemss = []
    while remain:
        g = iter(remain).__next__()
        items = get_cls(G, g)
        #for h in items:
        #    remain.remove(h)
        remain.difference_update(items)
        items = list(items)
        itemss.append(items)
    return itemss



def display(items):
    if not items:
        return
    n = len(items[0])
    smap = SMap()
    col = 0
    for A in items:
        smap[0,col] = A.shortstr()
        col += n+1
        if col > 100:
            print(smap, "\n")
            smap = SMap()
            col = 0
    print(smap, "\n")


#%See also:
#%``Classification of the Subgroups of the Two-Qubit Clifford Group''
#%    Eric Kubischta, Ian Teixeira
#% https://arxiv.org/abs/2409.14624

# https://www.maths.usyd.edu.au/u/don/code/Magma/SpConjugacy.pdf
# Conjugacy Classes in Finite Symplectic Groups, 2021
# D. E. Taylor

def main_sp():
    
    n = argv.get("n", 6)
    p = argv.get("p", 2)

    G = Algebraic.Sp(n, p)
    print("G = Sp(%d,%d)"%(n,p))
    #print("form:")
    #print(G.invariant_form)

    assert len(mulclose(G.gen)) == len(G)
    print("|G| =", len(G))

    if p==2:
        cgys = get_conjugacy_classes(G)
    else:
        cgys = slow_get_conjugacy_classes(G)

    print("char classes:", len(cgys))

    cgys.sort(key = len)
    cgys = [[G.to_ziporder(g) for g in cgy] for cgy in cgys] # <----- to_ziporder <---

    I = G.I
    F = sage.GF(p)
    chars = []
    for cgy in cgys:
        g = cgy[0]
        A = sage.matrix(F, n, n, g.A)
        char = A.characteristic_polynomial()
        #print(len(cgy), "-->", char, "\t", char.factor())
        chars.append(char)
        w = g.sum()
        best = g
        for g in cgy:
            if g.sum() < w:
                best = g
                w = g.sum()

        s = str(char.factor())
        s = s.replace(" ", "").replace("*", "")
        print(len(cgy), " & %s"%g.latex(), " & $%s$"%s, " & ? & ? ", r"\\")
        found = {char}
        for h in cgy[1:]:
            B = sage.matrix(F, n, n, h.A)
            dhar = B.characteristic_polynomial()
            assert char == dhar

    return

    for cgy in cgys:
        print(len(cgy))
        cgy.sort(key = lambda A:A.sum())
        #display([g for g in cgy if g.is_upper_triangular()])
        display(cgy)


def main_gap():

    n = argv.get("n", 10)
    p = argv.get("p", 2)
    
    cmd = r"""
    cgys := ConjugacyClasses(Sp(%s,%s));
    for cgy in cgys do m:=Representative(cgy); p:=CharacteristicPolynomial(m); Print(p,"\n"); od;
    quit;
    """ % (n, p)

    print("running gap...")
    print(cmd)
    f = open("/tmp/tmp.gap", "w")
    print(cmd, file=f)
    f.close()

    #lines = popen("gap --norepl /tmp/tmp.gap").readlines()
    data = popen("gap --norepl /tmp/tmp.gap").read()
    lines = data.split("\n")
    lines = [l for l in lines if l.startswith("x_1")]
    lines = [l.replace("x_1", "x").replace("Z(2)^0", "1") for l in lines]
    print("cgys:", len(lines))
    #print(" ".join(lines))

    polys = "x^3+x^2+x+1 x^4+1 x^5+x^4+x+1 x^6+x^4+x^2+1 x^7+x^6+x^5+x^4+x^3+x^2+x+1".split()
    polys += ["x^10+x^8+x^2+1", "x^12+x^8+x^4+1", "x^14+x^12+x^10+x^8+x^6+x^4+x^2+1"]
    for p in polys:
        print(lines.count(p), p)
    #return
    #print(lines.count("x^10+x^8+x^2+1"))
    #print(lines.count("x^12+x^8+x^4+1"))
    #print(lines.count("x^14+x^12+x^10+x^8+x^6+x^4+x^2+1"))





def main():
    n = argv.get("n", 4)
    p = argv.get("p", 2)
    GL = argv.get("GL", True)
    Sp = argv.get("Sp", False)

    if Sp:
        G = Algebraic.Sp(n, p)
        print("G = Sp(%d,%d)"%(n,p))
    elif GL:
        G = Algebraic.GL(n, p)
        print("G = GL(%d,%d)"%(n,p))
    else:
        G = Algebraic.SL(n, p)
        print("G = SL(%d,%d)"%(n,p))

    assert len(mulclose(G.gen)) == len(G)
    print("|G| =", len(G))

    if p==2:
        itemss = get_conjugacy_classes(G)
    else:
        itemss = slow_get_conjugacy_classes(G)

    print("char classes:", len(itemss))

    itemss.sort(key = len)
    print([len(items) for items in itemss], "==", len(G))
    assert sum([len(items) for items in itemss]) == len(G)

    A = Matrix([[1]])
    I2 = (parse("1.\n.1"))
    B = (parse("11\n1."))
    C = (parse("11.\n1.1\n.1."))
    D = (parse("11.\n111\n.1."))
    #u = parse("11.. .1.. ..11 ..1.") # 3360
    AC = A.direct_sum(C)
    AD = A.direct_sum(D)

    I = G.I
    F = sage.GF(p)
    chars = []
    for items in itemss:
        g = items[0]
        A = sage.matrix(F, n, n, g.A)
        char = A.characteristic_polynomial()
        #print(len(items), "-->", char, "\t", char.factor())
        chars.append(char)
        s = str(char.factor())
        s = s.replace(" ", "").replace("*", "")
        print(len(items), " & $%s$"%I.latex(), " & $%s$"%s, " & ? & ? ", r"\\")
        found = {char}
        for h in items[1:]:
            B = sage.matrix(F, n, n, h.A)
            dhar = B.characteristic_polynomial()
            assert char == dhar
            #if dhar not in found:
            #    print("\t", dhar, "\t", dhar.factor())
            #    found.add(dhar)

    return

    for i,items in enumerate(itemss):
        print(len(items), chars[i].factor())
        if len(items)!=1344:
            continue
        for u in items:
            #print("*" if u.sum()<=5 else ".", end="")
            #if u.sum() <= 5:
            if u[0,2]==u[0,3]==u[1,3]==u[2,0]==u[3,0]==u[3,1]==0:
                print(u.shortstr())
                print("$"+u.latex()+"$")
                print()


    
    return

    for items in itemss:
        smap = SMap()
        count = 0
        print(len(items))
        for i, A in enumerate(items):
            u = max( A[k,:k].max() for k in range(1,n) )
            if u:
                continue
            if p==2:
                smap[0,(n+1)*count] = A.shortstr()
            else:
                smap[0,(3*n+1)*count] = str(A) #A.shortstr()
            count += 1
        if count:
            print(smap)
        else:
            print("-"*n)
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



