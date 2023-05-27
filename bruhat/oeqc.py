#!/usr/bin/env python

from time import time
start_time = time()

import numpy

import bruhat.solve

from bruhat.solve import parse, span, shortstr, array2, solve, intersect, row_reduce, zeros2
from bruhat.solve import dot2, identity2, linear_independent
from bruhat.isomorph import Tanner, search
from bruhat.equ import Equ, quotient
from bruhat.action import Perm, Group, mulclose
from bruhat.util import factorize
from bruhat.qcode import QCode
from bruhat.argv import argv

#import bruhat.solve
#bruhat.solve.int_scalar = numpy.int8 # tiny !!
from qupy.ldpc import solve
solve.int_scalar = bruhat.solve.int_scalar
from qupy.ldpc.css import CSSCode


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


def get_autos(H):
    rows = [v for v in span(H) if v.sum()]
    V = array2(rows)
    print(shortstr(V))
    count = 0
    for f in find_isos(V, V):
        #print(f)
        print(".", end="", flush=True)
        count += 1
    print(count)


            
def gap_fmt(perm):
    cs = perm.cycles()
    ss = []
    for c in cs:
        if len(c)==1:
            continue
        c = [i+1 for i in c]
        ss.append(tuple(c))
    return ''.join(str(s) for s in ss)

def gap_code(perms):
    s = [gap_fmt(perm) for perm in perms]
    s = "Group(%s);"%(', '.join(s))
    return s

def get_autos_nauty(H):
    #print("get_autos_nauty", H.shape)
    rows = [v for v in span(H) if v.sum()]
    V = array2(rows)
    m, n = V.shape
    from pynauty import Graph, autgrp
    g = Graph(n+m) # bits + checks
    for bit in range(n):
        checks = [n+check for check in range(m) if V[check, bit]]
        g.connect_vertex(bit, checks)
    for check in range(m):
        bits = [bit for bit in range(n) if V[check, bit]]
        g.connect_vertex(n+check, bits)
    g.set_vertex_coloring([set(range(n)), set(range(n, m+n))])
    aut = autgrp(g)
    N = int(aut[1])
    #print("N =", N)

    gen = aut[0]
    items = list(range(n))
    perms = []
    for perm in gen:
        perm = perm[:n]
        #print(perm)
        perm = Perm(perm, items)
        perms.append(perm)
    #print(gap_code(perms))
    return N, perms

    #G = Group.generate(perms)
    G = mulclose(perms, verbose=True)
    N = len(G) # 322560, (C2 x C2 x C2 x C2) : A8
    print("autos:", N, factorize(N))

def css_autos_nauty(Hx, Hz):
    #print("get_autos_nauty", H.shape)
    Vx = array2([v for v in span(Hx) if v.sum()])
    Vz = array2([v for v in span(Hz) if v.sum()])
    mx, n = Vx.shape
    mz, n = Vz.shape
    from pynauty import Graph, autgrp
    g = Graph(n+mx+mz) # bits + checks
    for bit in range(n):
        checks = [n+check for check in range(mx) if Vx[check, bit]]
        checks += [n+mx+check for check in range(mz) if Vz[check, bit]]
        g.connect_vertex(bit, checks)
    for check in range(mx):
        bits = [bit for bit in range(n) if Vx[check, bit]]
        g.connect_vertex(n+check, bits)
    for check in range(mz):
        bits = [bit for bit in range(n) if Vz[check, bit]]
        g.connect_vertex(n+mx+check, bits)
    g.set_vertex_coloring([set(range(n)), set(range(n, n+mx)), set(range(n+mx,n+mx+mz))])
    print("autgrp...")
    aut = autgrp(g)
    N = int(aut[1])
    #print("N =", N)

    gen = aut[0]
    items = list(range(n))
    perms = []
    for perm in gen:
        perm = perm[:n]
        #print(perm)
        perm = Perm(perm, items)
        perms.append(perm)
    #print(gap_code(perms))
    return N, perms


def get_wenum(H0):
    n = H0.shape[1]
    words = {i:[] for i in range(n+1)}
    dist = [0 for i in range(n+1)]
    for v in span(H0):
        d = v.sum()
        dist[d] += 1
        words[d].append(v)
    dist = tuple(dist)
    return dist

def get_words(H, w):
    words = []
    for v in span(H):
        d = v.sum()
        if d==w:
            words.append(tuple(v))
    return words


def main_m24():
    # find some double cosets of 759 

    # Golay code:
    H0 = parse("""
    1...........11...111.1.1
    .1...........11...111.11
    ..1.........1111.11.1...
    ...1.........1111.11.1..
    ....1.........1111.11.1.
    .....1......11.11..11..1
    ......1......11.11..11.1
    .......1......11.11..111
    ........1...11.111...11.
    .........1..1.1.1..1.111
    ..........1.1..1..11111.
    ...........11...111.1.11
    """)

    # RM [[16,6,4]]
    H1 = parse("""
    11111111........
    ....11111111....
    ........11111111
    11..11..11..11..
    .11..11..11..11.
    """)

    m, n = H0.shape
    words = get_words(H0, 16)
    N = len(words)
    #print(words[0])
    codes = []
    for word in words:
        idxs = [idx for idx in range(n) if word[idx]==1]
        jdxs = [idx for idx in range(n) if word[idx]==0]
        Hx = zeros2(len(H1), n)
        Hx[:, idxs] = H1
        Hz = zeros2(len(H1) + len(jdxs), n)
        Hz[:len(H1), idxs] = H1
        for j, jdx in enumerate(jdxs):
            Hz[len(H1)+j, jdx] = 1
        #print("Hx =")
        #print(shortstr(Hx))
        #print("Hz =")
        #print(shortstr(Hz))
        code = CSSCode(Hx=Hx, Hz=Hz)
        #print(code)
        codes.append(code)
        #break

    counts = {i:0 for i in range(6)}
    for i in range(N):
     for j in range(i+1, N):
        H = intersect(codes[i].Lx, codes[j].Lx)
        counts[len(H)] += 1
        if len(H)>3:
            print(len(H), end=" ", flush=True)
     #print()
    print(counts)

    return

    from bruhat.gset import Perm

    _, perms = get_autos_nauty(H0)
    action = []
    for g in perms:
        print(g)
        idxs = []
        for w in words:
            w = tuple(w[g[i]] for i in range(n))
            idxs.append(words.index(w))
        #print(idxs)
        assert len(set(idxs)) == len(words)
        h = Perm(idxs)
        action.append(h)

    def get_orbit(i, j):
        bdy = [(i, j)]
        orbit = set(bdy)
        while bdy:
            _bdy = []
            for g in action:
                for (i0,j0) in bdy:
                    ij = g[i0], g[j0]
                    if ij not in orbit:
                        orbit.add(ij)
                        _bdy.append(ij)
            bdy = set(_bdy)
        return orbit

    remain = set((i,j) for i in range(N) for j in range(N))
    orbits = []
    while remain:
        i, j = iter(remain).__next__()
        orbit = get_orbit(i, j)
        assert len(orbit)%N == 0
        print("orbit:", len(orbit)//N)
        remain.difference_update(orbit)
        orbits.append(orbit)



def mulclose2(gen, verbose=False, maxsize=None):
    #els = {g.tobytes():g for g in gen}
    els = {shortstr(g):g for g in gen}
    bdy = dict(els)
    changed = True
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = {}
        for A in gen:
            for B in bdy:
                C = dot2(A, bdy[B])
                #key = C.tobytes()
                key = shortstr(C)
                if key not in els:
                    els[key] = C
                    _bdy[key] = C
                    if maxsize and len(els)>=maxsize:
                        if verbose:
                            print()
                        return els
        bdy = _bdy
    if verbose:
        print()
    return list(els.values())


def gap_matrix(A):
    m, n = A.shape
    matrix = []
    for row in A:
        line = []
        for x in row:
            line.append('Z(2)' if x else '0*Z(2)')
        line = "[%s]"%(','.join(line))
        matrix.append(line)
    matrix = "[%s]"%(','.join(matrix))
    return str(matrix)


def is_triorthogonal(H):
    # check triorthogonal
    m, n = H.shape
    for i in range(m):
     for j in range(i+1, m):
        if (H[i]*H[j]).sum() % 2:
            return False
        for k in range(j+1, m):
          if (H[i]*H[j]*H[k]).sum() % 2:
            return False
    return True


#def find_perm_logops(code, perms):
    



def main_triorthogonal():
    H = parse("""
    11111111........
    ....11111111....
    ........11111111
    11..11..11..11..
    .11..11..11..11.
    """)
    assert is_triorthogonal(H)

    nn = 2*n
    _, perms = get_autos_nauty(H)
    G = []
    for g in perms:
        A = zeros2((n, n))
        for i in range(n):
            j = g[i]
            A[i, j] = 1
        #print(shortstr(A))
        #print()
        G.append(A)

    if 0:
        K = mulclose2(G, verbose=True, maxsize=1000000)
        assert len(K) == 322560 # yes
    
    if 0:
        code = QCode.build_css(H, H)
        L0 = code.get_logops()
        print(code.get_params(10))
    
        A = (dot2(L0, L0.transpose()))
        print(shortstr(A))

    code = CSSCode(Hx=H, Hz=H)
    k = code.k
    print(code)

    gen = []

    # Hadamard
    Lx, Lz = code.Lz, code.Lx
    gx = dot2(Lx, code.Lz.transpose())
    gz = dot2(Lz, code.Lx.transpose())
    g = zeros2((2*k, 2*k))
    g[:k, k:] = gx
    g[k:, :k] = gz
    gen.append(g)

    # S-gate
    g = identity2(2*k)
    g[:k, k:] = gx
    gen.append(g)

    for A in G:
        Lx = numpy.dot(code.Lx, A)
        Lz = numpy.dot(code.Lz, A)
        gx = dot2(Lx, code.Lz.transpose())
        gz = dot2(Lz, code.Lx.transpose())
        g = zeros2((2*k, 2*k))
        g[:k, :k] = gx
        g[k:, k:] = gz
        gen.append(g)

    vs = []
    for i, g in enumerate(gen):
        print("m%d := %s;;"%(i, gap_matrix(g)))
        vs.append("m%d"%i)
    print("G := Group(%s);" % (','.join(vs)))

    N = 322560

    if 0:
        K = mulclose2(gen, verbose=True, maxsize=N+1)
        assert len(K) <= N # fail
        assert N%len(K) == 0
        print(len(K), N//len(K))

    #K = mulclose2(gen, verbose=True)

    
def main_rm():
    H = parse("""
    11111111........
    ....11111111....
    ........11111111
    11..11..11..11..
    .11..11..11..11.
    """)
    print(H.shape)
    
    _, perms = get_autos_nauty(H)

    #G = mulclose(perms, verbose=True)
    #N = len(G) # 322560, (C2 x C2 x C2 x C2) : A8
    #print("autos:", N, factorize(N))

    v = '1111............'
    n = len(v)
    m = len(perms[0].items)
    items = list(range(m-n, m))
    print(items)
    perms = [g.restrict(items) for g in perms]
    items = list(range(n))
    assert 0, "BROKEN below: reverse m and n"
    #perms = [Perm(perm.perm, items) for perm in perms]
    #perms = [Perm(dict((k-m+n,v-m+n) for (k,v) in g.perm.items()), items) for g in perms]
    perms = [dict((k-m+n,v-m+n) for (k,v) in g.perm.items()) for g in perms]
    perms = [[g[i] for i in items] for g in perms]
    print(perms)
    orbit = {v}
    bdy = {v}
    while bdy:
        _bdy = set()
        for g in perms:
            for v in bdy:
                v1 = ''.join(v[i] for i in g)
                if v1 not in orbit:
                    _bdy.add(v1)
                    orbit.add(v1)
        bdy = _bdy
    print(len(orbit))
    orbit = list(orbit)
    orbit.sort()
    orbit = "\n".join(orbit)
    J0 = parse(orbit)
    J0 = row_reduce(J0)
    print(shortstr(J0), J0.shape)
    JX = zeros2(*J0.shape, 2)
    JX[:,:,0] = J0
    JZ = zeros2(*J0.shape, 2)
    JZ[:,:,1] = J0
    J = numpy.concatenate((JX, JZ))
    print(J.shape)

    code = QCode.build_gauge(J)
    print(code)
    L = code.get_logops()
    print(L.shape)



def main_golay():
    """
    012345678901234567890123456789
    1...........11...111.1.1
    .1...........11...111.11
    ..1.........1111.11.1...
    ...1.........1111.11.1..
    ....1.........1111.11.1.
    .....1......11.11..11..1
    ......1......11.11..11.1
    .......1......11.11..111
    ........1...11.111...11.
    .........1..1.1.1..1.111
    ..........1.1..1..11111.
    ...........11...111.1.11
    """

    # Golay code
    H0 = parse("""
    1...........11...111.1.1
    .1...........11...111.11
    ..1.........1111.11.1...
    ...1.........1111.11.1..
    ....1.........1111.11.1.
    .....1......11.11..11..1
    ......1......11.11..11.1
    .......1......11.11..111
    ........1...11.111...11.
    .........1..1.1.1..1.111
    ..........1.1..1..11111.
    ...........11...111.1.11
    """)

    if argv.puncture:
        for i in range(1, 5):
            H = H0[i:, i:]
            print(shortstr(H), H.shape)
            code = QCode.build_css(H, H)
            print(code.get_params(max_mk=24))
        return

    n = H0.shape[1]
    words = {i:[] for i in range(n+1)}
    dist = [0 for i in range(n+1)]
    for v in span(H0):
        d = v.sum()
        dist[d] += 1
        words[d].append(v)
    dist = tuple(dist)
    print(dist)

    codes = []
    for word in words[16]:
        rows = []
        for u in words[8]:
            if numpy.alltrue(u*word == u):
                rows.append(u)
        H = array2(rows)
        H = row_reduce(H)
        code = QCode.build_css(H, H)
        codes.append(code)
    print(len(codes))

    counts = {i:0 for i in range(49)}
    L = codes[0].get_logops()
    print(shortstr(L), L.shape)
    left = codes[0].get_all_logops()
    print(left.shape)
    for code in codes:
        #w = intersect(codes[0].flatH, code.flatH)
        right = code.get_all_logops()
        w = intersect(left, right)
        #w1 = intersect(w, codes[0].flatH
        w1 = intersect(codes[0].get_logops(), code.get_logops())
        counts[len(w1)] += 1
        print(len(w1), end=' ', flush=True)
    print()
    print([(v,k) for (k,v) in counts.items() if v])

    return

    Hs = []
    for i in range(12):
        rows = [j for j in range(12) if H0[j, i+12]==0]
        H1 = H0[rows]
        Hs.append(H1)
        #print(shortstr(H1), H1.shape)
        #print()

    for Ha in Hs:
      for Hb in Hs:
        Hab = intersect(Ha, Hb)
        print(Hab.shape[0], end=' ')
      print()

    #get_autos_nauty(H)

def parse_decl(decl):
    decl = decl.split()
    xop = '\n'.join(line.replace('X','1') for line in decl if 'X' in line)
    zop = '\n'.join(line.replace('Z','1') for line in decl if 'Z' in line)
    #print(xop)
    Hx = bruhat.solve.parse(xop)
    Hz = bruhat.solve.parse(zop)
    #Hx = linear_independent(Hx)
    #Hz = linear_independent(Hz)
    #print(Hx)
    return CSSCode(Hx=Hx, Hz=Hz)


def build_code():
    code = parse_decl("""
    .X.XX.X...X.X.X..X...
    X....X.XX..X.X.X..X..
    X..X.XX...XX..X...X..
    .......XXX...X.XX..X.
    X.XXX.X.XX.X.X...X.XX
    .Z.ZZ.Z...Z.Z.Z..Z...
    Z....Z.ZZ..Z.Z.Z..Z..
    Z..Z.ZZ...ZZ..Z...Z..
    .ZZ.Z.......Z....Z..Z
    Z.ZZZ.Z.ZZ.Z.Z...Z.ZZ
    ......Z.....Z.Z..Z...
    Z......ZZ.........Z..
    .......ZZ.......Z..Z.
    ..Z.........Z....Z...
    .Z.....Z..Z.....Z.Z..
    ..Z.................Z
    ..Z......Z.........ZZ
    Z.......Z..Z.Z.......
    .....Z....Z...Z...Z..
    .Z..........Z........
    """) # [[21,1,2]] # Aut=16384

    code = parse_decl("""
    .X.XX.X...X.X.X..X...
    X....X.XX..X.X.X..X..
    X..X.XX...XX..X...X..
    .......XXX...X.XX..X.
    X.XXX.X.XX.X.X...X.XX
    .Z.ZZ.Z...Z.Z.Z..Z...
    Z....Z.ZZ..Z.Z.Z..Z..
    Z..Z.ZZ...ZZ..Z...Z..
    .ZZ.Z.......Z....Z..Z
    Z.ZZZ.Z.ZZ.Z.Z...Z.ZZ
    ......Z.....Z.Z..Z...
    Z......ZZ.........Z..
    .......ZZ.......Z..Z.
    ..Z.........Z....Z...
    .Z.....Z..Z.....Z.Z..
    ..Z......Z.........ZZ
    Z.......Z..Z.Z.......
    .....Z....Z...Z...Z..
    """) # 

    print(code)
    print(code.distance())
    print("is_triorthogonal", is_triorthogonal(code.Hx))

    dump_transverse(code.Hx, code.Lx)

    N, perms = get_autos_nauty(code.Hx)
    print(N, len(perms))

    G = mulclose(perms, verbose=True); print()
    gen = []
    #for g in perms:
    for g in G:
        idxs = [g[i] for i in range(code.n)]
        Hz = code.Hz[:, idxs]
        H = intersect(Hz, code.Hz)
        if len(H) == len(Hz):
            print('.', end='', flush=True)
            gen.append(g)
        #print(len(code.Hz), len(H))
    print(len(gen))

    #G = mulclose(gen, verbose=True)
    #print(len(G))

    #N, perms = get_autos_nauty(code.Hz)
    #print(N, len(perms))

    #N, perms = css_autos_nauty(code.Hx, code.Hz)
    #print(N, len(perms))


def dump_transverse(Hx, Lx):
    from qupy import CSSLO
    Eq, SX,LX,SZ,LZ = CSSLO.CSSCode(Hx, Lx)
    t = 3
    N = 1<<t
    zList,qList, V, K_M = CSSLO.comm_method(Eq, SX, LX, SZ, t, compact=True, debug=False)
    for z,q in zip(zList,qList):
        print("#", CSSLO.CP2Str(2*q,V,N),"=>",CSSLO.z2Str(z,N))
    print()

def build_trinal():

    jdx = argv.get("idx", 0)
    from bruhat.small_triorthogonal import codes
    for idx, code in enumerate(codes):
        Hx, Hz, Lx, Lz, comment = code
        code = CSSCode(Hx=Hx, Hz=Hz, Lx=Lx, Lz=Lz)
        print(idx, code)
        print(code.distance())
        print(code.longstr())
        print("is_triorthogonal:", is_triorthogonal(code.Hx))
        if Lx is None:
            Lx = code.Lx
        dump_transverse(Hx, Lx)
        if idx:
            break


    return

    if 0:
        N, perms = css_autos_nauty(code.Hx, code.Hz)
        print(gap_code(perms))
        print("|G| =", N)
        for perm in perms:
            print(perm)
        #G = mulclose(perms, verbose=True)
        #print(len(G))

    print(code.distance())
    print(code.longstr())

    for v in span(Lz):
        if v.sum() == 0:
            continue
        for w in span(Hz):
            vw = (v+w)%2
            if vw.sum() == 1:
                print(v)
                print(w)
                print(vw)
                assert 0

        
#def find_iso():
#    # broken: must take span of H
#    codes = list(read_codes())
#    print(len(codes))
#
#    #codes = codes[:100]
#    #codes = [H.tobytes() for H in codes]
#    equs = quotient(codes, is_iso, verbose=True)
#    found = set()
#    for equ in equs:
#        found.add(equ.top)
#    print()
#    print("distinct codes:", len(found))
#    for equ in found:
#        print(shortstr(equ.items[0]))
#        print()



def make_colour():
    from qcode import Geometry, get_adj

    key = argv.get("key", (3, 8))
    idx = argv.get("idx", 10)

    gens = (0,), (1,), (2,)
    f, e, v = gens

    name = argv.name

    for idx in [idx]:
        print("idx =", idx)
        geometry = Geometry(key, idx, False)
    
        hgens = []
        graph = geometry.build_graph(hgens=hgens)
    
        graph = graph.compress()
        print(len(graph))
    
        faces = graph.components([(1,), (2,)])
        print("faces:", [len(c) for c in faces], len(faces))
    
        verts = graph.components([(0,), (1,)])
        print("verts:", [len(c) for c in verts], len(verts))

        H = get_adj(faces, verts)
        H = linear_independent(H)
        print(H.shape)
        assert dot2(H, H.transpose()).sum() == 0

        if name:
            f = open(name, 'w')
            print(shortstr(H), file=f)
            f.close()


def make_ramified():
    from qcode import Geometry, get_adj

    key = (3, 8)
    idx = argv.get("idx", 10)

    gens = (0,), (1,), (2,)
    f, e, v = gens

    for idx in [idx]:
        print("idx =", idx)
        geometry = Geometry(key, idx, False)
    
        hgens = [(1,2)*4]
        #hgens = []
        graph = geometry.build_graph(hgens=hgens)
    
        graph = graph.compress()
        print(len(graph))
    
        faces = graph.components([(1,), (2,)])
        print("faces:", [len(c) for c in faces], len(faces))
    
        verts = graph.components([(0,), (1,)])
        print("verts:", [len(c) for c in verts], len(verts))

        colours = graph.components([e, v, f+e+v+e+f])
        print("colours:", [len(c) for c in colours], len(colours))

        J = get_adj(faces, colours)
        print(shortstr(J))
    
        H = get_adj(faces, verts)
        HHt = dot2(H, H.transpose())
        if HHt.sum():
            print("Not self-dual")
            #print(shortstr(HHt))
            #code = CSSCode(Gx=H, Gz=H) # fail..
            #print(code)
            continue
        #print(shortstr(H))
    
        code = CSSCode(Hx=H, Hz=H)
        print(code)
    
        from bruhat.hecke import get_distance
        print(get_distance(H, code.Lx))
        #print(code.distance())

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


