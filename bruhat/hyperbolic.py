#!/usr/bin/env python3

import sys
from random import randint, seed, choice
from time import sleep
from functools import reduce
from operator import mul, add

import numpy
from numpy import concatenate

#from bruhat.gset import Perm, Group, Coset, mulclose # FAIL
from bruhat.action import Perm, Group, Coset, mulclose, close_hom
from bruhat.util import cross
from bruhat.smap import SMap
from bruhat.argv import argv

from bruhat.todd_coxeter import Schreier



    
def make_random_modular():
    ngens = 4
    a, ai, b, bi = range(ngens)
    rels = [ (ai, a), (bi, b), (a,)*2, (a,b,)*3, ] # modular group 
    rels += [ (b,)*7, ] # 2,3,7 triangle group

    for _ in range(10000):
        graph = Schreier(ngens, rels)
        hgens = []
        for i in range(2):
            gen = tuple(randint(0, ngens-1) for k in range(20))
            hgens.append(gen)
        if graph.build(hgens, maxsize=2000):
            if len(graph) < 3:
                continue
            #print("hgens:", hgens)
            gens = graph.get_gens()
            G = mulclose(gens, maxsize=10000)
            if len(G) < 10000:
                print(len(graph), end=" ", flush=True)
                print(len(G))


def make_random_55():

    #seed(4)

    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    i_rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
    ]
    a1 = (b,a)
    b1 = (a,c)
    rel1 = (3*a1+b1)*3 # Bring's curve

    for _ in range(10000):
        rels = list(i_rels)
        for i in range(1):
            #gen = tuple([a, b, c][randint(0, 2)] for k in range(30))
            gen = ()
            for k in range(50):
                gen += choice([(a,b), (a,c), (b,c)])
            #gen = rel1
            rels.append(gen)
        graph = Schreier(ngens, rels)
        if graph.build(maxsize=1040000):
            print(".", flush=True, end="")
            n = len(graph)
            if n <= 1320:
                continue
            #print("rels:", rels)
            #graph.show()
            print(len(graph.neighbors))
            try:
                gens = graph.get_gens()
            except AssertionError:
                print("XXX")
                continue
            G = mulclose(gens, maxsize=10000)
            if len(G) < 10000:
                print(gen)
                print(len(graph), end=" ", flush=True)
                print(len(G))
            print()


def make_random_homogeneous_refl():

    seed(1)

    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    i_rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
    ]

    while 1:
        rels = []
        for i in range(4):
            gen = ()
            for k in range(20):
                gen += choice([(a,b), (a,c), (b,c)])
            rels.append(gen)

        # G/H = 720/24 = 30, N(H)=24
        _rels = [(0, 2, 2, 4, 0, 4, 2, 4, 2, 4, 0, 4, 0, 4, 2,
            4, 0, 2, 0, 4, 2, 4, 0, 4, 2, 4, 0, 2, 0, 2, 0, 2, 2, 4, 2, 4, 2, 4, 2, 4),
            (0, 2, 2, 4, 0, 4, 2, 4, 0, 2, 0, 4, 0, 2, 0, 4, 0, 4,
            0, 2, 0, 2, 0, 2, 0, 4, 0, 4, 2, 4, 0, 4, 0, 2, 0, 2, 2, 4, 0, 2),
            (0, 2, 2, 4, 2, 4, 2, 4, 2, 4, 0, 2, 0, 4, 0, 4, 0, 4,
            0, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2, 2, 4, 0, 4, 0, 2, 0, 4, 0, 2),
            (0, 2, 0, 4, 0, 4, 2, 4, 0, 4, 0, 2, 0, 2, 2, 4, 0, 4,
            0, 4, 0, 2, 0, 4, 0, 2, 0, 2, 2, 4, 0, 2, 2, 4, 2, 4, 0, 4, 0, 2)]

        # G/H = 9720/324 = 30, N(H)=1944
        _rels = [(0, 4, 0, 2, 0, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2,
        4, 0, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 2, 0, 2, 0, 2, 2,
        4, 0, 2, 0, 4, 0, 4), (0, 4, 2, 4, 0, 2, 0, 2, 0, 2,
        2, 4, 0, 2, 0, 2, 2, 4, 2, 4, 0, 4, 2, 4, 0, 4, 2, 4,
        0, 4, 0, 2, 2, 4, 0, 2, 2, 4, 0, 2), (0, 4, 0, 4, 0,
        4, 0, 2, 0, 4, 0, 4, 2, 4, 0, 4, 0, 2, 0, 2, 0, 4, 0,
        2, 0, 2, 0, 2, 0, 4, 0, 4, 0, 4, 2, 4, 0, 2, 2, 4), (2,
        4, 2, 4, 0, 4, 2, 4, 0, 2, 0, 4, 0, 4, 2, 4, 2, 4, 2,
        4, 2, 4, 2, 4, 0, 4, 0, 4, 2, 4, 2, 4, 2, 4, 0, 4, 0, 2, 2, 4)]


        graph = Schreier(ngens, i_rels)
        if not graph.build(rels, maxsize=10400):
            print(choice("/\\"), end="", flush=True)
            continue

        n = len(graph)
        if n <= 24:
            continue

        print("\n[%d]"%n) #, flush=True, end="")
        print(len(graph.neighbors)) # how big did the graph get? compare with maxsize above

        try:
            gens = graph.get_gens()
        except AssertionError:
            print("** FAIL **")
            continue

        N = 10000
        perms = mulclose(gens, maxsize=N)
        if len(perms) >= N:
            print("|G| = ?")
            continue

        print(rels)
        items = list(range(n))
        G = Group(perms, items)
        G.gens = gens
        print("|G| =", len(G))
        print()

        break

    rels = [reduce(mul, [gens[i] for i in rel]) for rel in rels]
    H = Coset(mulclose(rels), items)
    print(len(H))
    print(len(G) / len(H))

    N = []
    for g in G:
        lhs = g*H
        rhs = H*g
        if lhs == rhs:
            N.append(g)
    print(len(N))


def conj_rels(_rels):
    rels = []
    for rel in _rels:
        for i in range(len(rel)):
            r = rel[i:] + rel[:i] # cyclic permutation
            rels.append(r)
    return rels

def reduce_word(_rels, word):
    rels = []
    for rel in _rels:
        for i in range(len(rel)):
            r = rel[i:] + rel[:i] # cyclic permutation
            rels.append(r)
    done = False
    while not done:
        done = True
        n = len(word)
        for i in range(n):
            for rel in rels:
                m = len(rel)
                if word[i:i+m] == rel:
                    word = word[:i] + word[i+m:]
                    done = False
                    break
            else:
                continue
            break
    return word


def make_random_homogeneous():
    # use rotation group

    #seed(1)

    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    i_rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*5, (b,)*2, (c,)*4, (a, b, c),
    ]

    while 1:
        rels = []
        for i in range(argv.get("nwords", 2)):
            gen = tuple(choice([a, b, c]) for k in range(argv.get("wordlen", 100)))
            rels.append(gen)
            gen = (ci,) + gen + (c,)
            rels.append(gen)
            gen = (ci,) + gen + (c,)
            rels.append(gen)

        rels = [reduce_word(i_rels, rel) for rel in rels]

        graph = Schreier(ngens, i_rels)
        if not graph.build(rels, maxsize=10400):
            print(choice("/\\"), end="", flush=True)
            continue

        n = len(graph)
        if n <= argv.get("minsize", 12):
            print('.', end="", flush=True)
            continue

        print("[%d]"%n, flush=True, end="")

        # how big did the graph get? compare with maxsize above
        print("(%d) "%len(graph.neighbors), flush=True, end="")

        try:
            gens = graph.get_gens()
        except AssertionError:
            print("** FAIL **")
            continue

        items = list(range(n))
        print()
        break

        N = argv.get("N", 10000)
        perms = mulclose(gens, maxsize=N)
        if len(perms) >= N:
            #print("|G| = ?")
            continue

        print()
        print(rels)
        G = Group(perms, items)
        G.gens = gens
        print("|G| =", len(G))


        rels = [reduce(mul, [gens[i] for i in rel]) for rel in rels]
        H = Coset(mulclose(rels), items)
        print("|H| =", len(H))
        assert len(G)%len(H) == 0
        print("[H:G] = ", len(G) // len(H))
    
        N = []
        for g in G:
            lhs = g*H
            rhs = H*g
            if lhs == rhs:
                N.append(g)
        assert len(N)%len(H) == 0
        print("[H:N(H)] =", len(N)//len(H))

    print(rels)

    a, ai, b, bi, c, ci = gens
    assert (a**5).is_identity()
    assert (b**2).is_identity()
    assert (c**4).is_identity()
    assert (a*b*c).is_identity()

    print(a.fixed())
    print(b.fixed())
    print(c.fixed())

    print(a.orbits())
    print(b.orbits())
    print(c.orbits())
    #for face in items:
        

def make_homogeneous():
    # use rotation group

    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    i_rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*5, (b,)*2, (c,)*4, (a, b, c),
    ]
    _rels = [(2, 4, 0, 4, 0, 0, 4, 0, 4, 4, 4, 2, 2, 2, 2,
        2, 0, 4, 2, 0, 4, 4, 2, 4, 0, 4, 0, 4, 4, 0, 2, 0, 4,
        2, 4, 0, 2, 0, 4, 4), (0, 0, 0, 2, 2, 4, 2, 0, 4, 2,
        0, 4, 0, 4, 4, 0, 0, 0, 4, 4, 2, 4, 4, 4, 4, 4, 0, 2,
        2, 2, 0, 4, 4, 0, 4, 0, 2, 4, 0, 4), (2, 0, 4, 4, 0,
        2, 4, 4, 2, 0, 0, 2, 0, 2, 2, 2, 2, 2, 0, 4, 2, 2, 0,
        0, 2, 4, 0, 4, 4, 2, 2, 4, 2, 4, 0, 0, 2, 0, 0, 0)]

    # no fixed points, F,E,V = 12,30,16
    rels = [(2, 0, 4, 4, 0, 0, 4, 4, 0, 4, 0, 0, 0, 0, 4,
        2, 0, 4, 2, 0, 0, 4, 4, 0, 4, 2, 0, 4, 4, 2, 0, 2, 0,
        4, 4, 4), (5, 2, 0, 4, 4, 0, 0, 4, 4, 0, 4, 0, 0, 0,
        0, 4, 2, 0, 4, 2, 0, 0, 4, 4, 0, 4, 2, 0, 4, 4, 2, 0,
        2, 0), (5, 5, 2, 0, 4, 4, 0, 0, 4, 4, 0, 4, 0, 0, 0,
        0, 4, 2, 0, 4, 2, 0, 0, 4, 4, 0, 4, 2, 0, 4, 4, 2, 0,
        2, 0, 4), (0, 0, 4, 0, 0, 4, 0, 4, 4, 4, 0, 0, 4, 0,
        4, 2, 0, 0, 0, 0, 4, 0, 4, 0, 4, 4, 0, 4, 4, 2, 4, 4,
        4, 2, 0, 4, 2, 4, 2, 4, 4, 0, 4, 4, 4), (5, 0, 0, 4,
        0, 0, 4, 0, 4, 4, 4, 0, 0, 4, 0, 4, 2, 0, 0, 0, 0, 4,
        0, 4, 0, 4, 4, 0, 4, 4, 2, 4, 4, 4, 2, 0, 4, 2, 4, 2,
        4, 4, 0), (5, 5, 0, 0, 4, 0, 0, 4, 0, 4, 4, 4, 0, 0,
        4, 0, 4, 2, 0, 0, 0, 0, 4, 0, 4, 0, 4, 4, 0, 4, 4, 2,
        4, 4, 4, 2, 0, 4, 2, 4, 2, 4, 4, 0, 4)]


    rels = [reduce_word(i_rels, word) for word in rels]
    print(rels)

    graph = Schreier(ngens, i_rels)
    if not graph.build(rels, maxsize=10400):
        assert 0

    words = graph.get_words()
    for w in words:
        assert reduce_word(conj_rels(i_rels + rels), w) == w
    #    print(w)
    print("words:", len(words))

    gens = graph.get_gens()

    N = argv.get("N", 10000)
    perms = mulclose(gens, maxsize=N)
    if len(perms) >= N:
        print("|G| = ?")
    else:
        print("|G| = ", len(perms))

    a, ai, b, bi, c, ci = gens
    assert (a**5).is_identity()
    assert (b**2).is_identity()
    assert (c**4).is_identity()
    assert (a*b*c).is_identity()

    print(a.fixed())
    print(b.fixed())
    print(c.fixed())

    print("faces:", len(a.orbits()))
    print("edges:", len(b.orbits()))
    print("vertices:", len(c.orbits()))

    from bruhat.disc import render_group
    rels = conj_rels(rels)
    render_group(5, 2, 4, words, rels)

    return
        
#    edges = set()
#    for tile in items:
#        for g in [a, ai, c, ci]:
#            edge = (tile, g(tile))
#            edges.append(edge)

    def intersect(lhs, rhs):
        return bool(set(lhs).intersection(rhs))
    def eq(lhs, rhs):
        for i in lhs:
            if i not in rhs:
                return False
        for i in rhs:
            if i not in lhs:
                return False
        assert len(lhs)==len(rhs)
        return True
    def search(tiles, tiless):
        for rhs in tiless:
            if eq(tiles, rhs):
                return True

    faces = [tuple(face) for face in a.orbits()]
    print("faces:", faces)

    assert intersect(faces[0], faces[0])
    assert not intersect(faces[0], faces[1])

    flookup = {}
    for face in faces:
        for i in face:
            flookup[i] = face
    for face in faces:
        print(faces.index(face), ":", end=" ")
        for i in face:
            print(faces.index(flookup[c(i)]), end=" ")
        print()
            

def make_euclidean():

    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*4, (b,)*2, (c,)*4, 
        (a, c, b), # Note: this relator gives better (more local) "vertices" than (a, b, c) !
    ]

    L = argv.get("L", 2)
    gamma = [(ci, a)*L,]
    H = [(b,)]
    graph = Schreier(ngens, rels + gamma)
    result = graph.build(H, maxsize=10000)
    assert result
    print("n =", len(graph))

#    labels = graph.labels
#    neighbors = graph.neighbors
#    for i, nbd in enumerate(graph.neighbors):
#        if labels[i] == i:
#            print(i, nbd)

    G = graph.get_group()
    print("|G| =", len(G))
    a, ai, b, bi, c, ci = G.gens

    assert (a**4).is_identity()
    assert (b**2).is_identity()
    assert (c**4).is_identity()
    assert (a*c*b).is_identity()

    print("vertices", a.orbits())
    print("edges   ", b.orbits())
    print("faces   ", c.orbits())
    print("fixed vertices", a.fixed())
    print("fixed edges   ", b.fixed())
    print("fixed faces   ", c.fixed())

    print("faces:", len(a.orbits()))
    print("edges:", len(b.orbits()))
    print("vertices:", len(c.orbits()))

    faces = [tuple(o) for o in c.orbits()]

    flookup = {}
    for face in faces:
        for i in face:
            flookup[i] = face
    for face in faces:
        print(faces.index(face), ":", end=" ")
        for i in face:
            print(faces.index(flookup[b(i)]), end=" ")
        print()
            
    I = G.identity
    def get_face(row, col):
        g = reduce(mul, [b*c*c]*row + [c*b*c]*col or [I])
        return g

    print()
    dj = 4
    di = 3
    star = 0
    smap = SMap()
    for row in range(L):
      for col in range(L):
        g = get_face(row, col)
        i, j = di*row, 2*dj*col
        smap[i, j] = str(g(star))
        smap[i, j+3] = str((c*g)(star))
        smap[i+1, j+3] = str((c*c*g)(star))
        smap[i+1, j] = str((c*c*c*g)(star))
    print(smap)


def make_hyperbolic_525():

    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*5, (b,)*2, (c,)*5, 
        (a, c, b), # Note: this relator gives better (more local) "vertices" than (a, b, c) !
        #(a, b, c),
    ]

    rel = (ci, a)*3 # Bring's curve
    gamma = [rel]
    H = [(b,)]
    graph = Schreier(ngens, rels + gamma)
    result = graph.build(H, maxsize=10000)
    assert result
    print("n =", len(graph))

#    labels = graph.labels
#    neighbors = graph.neighbors
#    for i, nbd in enumerate(graph.neighbors):
#        if labels[i] == i:
#            print(i, nbd)

    G = graph.get_group()
    print("|G| =", len(G))
    a, ai, b, bi, c, ci = G.gens

    assert (a**5).is_identity()
    assert (b**2).is_identity()
    assert (c**5).is_identity()
    assert (a*c*b).is_identity()

    print("vertices", a.orbits())
    print("edges   ", b.orbits())
    print("faces   ", c.orbits())
    print("fixed vertices", a.fixed())
    print("fixed edges   ", b.fixed())
    print("fixed faces   ", c.fixed())

    print("faces:", len(a.orbits()))
    print("edges:", len(b.orbits()))
    print("vertices:", len(c.orbits()))

    faces = [tuple(o) for o in c.orbits()]

    flookup = {}
    for face in faces:
        for i in face:
            flookup[i] = face
    for face in faces:
        print(faces.index(face), ":", end=" ")
        for i in face:
            print(faces.index(flookup[b(i)]), end=" ")
        print()


def make_random_524():
    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*5, (b,)*2, (c,)*4, 
        (a, c, b), # Note: this relator gives better (more local) "vertices" than (a, b, c) !
        #(a, b, c),
    ]

    while 1:
        gamma = [tuple(choice([a, b, c]) for i in range(50)) for j in range(3)]
        gamma = [reduce_word(rels, rel) for rel in gamma]
        H = []
    
        graph = Schreier(ngens, rels + gamma)
        result = graph.build(H, maxsize=10000)
        if not result:
            print("/", end="", flush=True)
            continue
        n = len(graph)
        if n >= 120:
            print("[%d]"%n, end="")

        if n > 320:
            break

        print(".", end="", flush=True)
    print()
    print(rel)
    print("n =", n)


def make_524():
    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*5, (b,)*2, (c,)*4, 
        (a, c, b), # Note: this relator gives better (more local) "vertices" than (a, b, c) !
        #(a, b, c),
    ]

    rel = (ai,b,a,b)*3 # n=120, Bring's curve
    #rel = (2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0) # n = 160
    #rel = (2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0) # 320

    gamma = [rel]
    if argv.wrap:
        H = [(b,)]
        H = [(a,a,a,a,b,a)]
        H = [(c,c)]
    else:
        H = []

    graph = Schreier(ngens, rels + gamma)
    result = graph.build(H, maxsize=10000)
    assert result
    print("n =", len(graph))

    return

    G = graph.get_group()
    print("|G| =", len(G))
    a, ai, b, bi, c, ci = G.gens

    assert (a**5).is_identity()
    assert (b**2).is_identity()
    assert (c**4).is_identity()
    assert (a*c*b).is_identity()

    vertices = c.orbits()
    edges = b.orbits()
    faces = a.orbits()

    print("vertices:")
    for idx, item in enumerate(vertices):
        print("\t", idx, item)
    print("edges:   ")
    for idx, item in enumerate(edges):
        print("\t", idx, item)
    print("faces:   ")
    for idx, item in enumerate(faces):
        print("\t", idx, item)

    print("faces:", len(a.orbits()))
    print("edges:", len(b.orbits()))
    print("vertices:", len(c.orbits()))

    print("fixed vertices", c.fixed())
    print("fixed edges   ", b.fixed())
    print("fixed faces   ", a.fixed())

    faces = [tuple(o) for o in a.orbits()]

    flookup = {}
    for face in faces:
        for i in face:
            flookup[i] = face
    for face in faces:
        print(faces.index(face), ":", end=" ")
        for i in face:
            print(faces.index(flookup[b(i)]), end=" ")
        print()

    
    twists = [edge[0] for edge in edges if len(edge)==1]
    H = numpy.empty((len(faces), len(vertices)), dtype=object)
    H[:] = '.'
    for i in range(len(faces)):
        lhs = set(faces[i])
        for j in range(len(vertices)):
            rhs = set(vertices[j])
            meet = lhs.intersection(rhs)
            if meet:
                if meet.intersection(twists):
                    H[i, j] = 'Y'
                else:
                    H[i, j] = 'I'
    print('\n'.join(''.join(row) for row in H))

    labels = {}
    I = G.identity
    words = list(graph.get_words())
    for word in words:
        g = reduce(mul, [G.gens[w] for w in word], I)
        #print(word, g(0))
        label = str(g(0))
        labels[word] = label
#        if H:
#            labels[word+(2,)] = label
        #w = word + rel
        #labels[w] = label
#        for wrd in words:
#            iwrd = tuple({0:1, 1:0, 2:3, 3:2, 4:5, 5:4}[w] for w in reversed(wrd))
#            assert graph.follow_path(0, wrd+iwrd) == 0
#            #w = word + iwrd+rel+wrd
#            #labels[w] = label
#            w = word + iwrd+rel+wrd+(2,)
#            labels[w] = label

#    from bruhat.disc import render_group
#    render_group(5, 2, 4, labels=labels, name="output", maxsize=1000)
#    render_group(5, 2, 4, labels=labels, name="bring", maxsize=1000)
#    return


def make_bring():

    # ------------------------------------------------

    # Bring's curve rotation group
    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), (a,)*5, (b,)*2]
    rels += [ (a,a,a,b)*3 ]
    graph = Schreier(2*ngens, rels)
    #graph.DEBUG = True
    graph.build()
    assert len(graph) == 60 # == 12 * 5

    # --------------- Group ------------------
    G = graph.get_group()
    #burnside(G)
    assert len(G.gens) == 4
    ai, a, bi, b = G.gens
    assert a==~ai
    assert b==~bi

    assert a.order() == 5
    assert b.order() == 2
    H = Group.generate([a])

    faces = G.left_cosets(H)
    assert len(faces) == 12
    #hom = G.left_action(faces)

    ab = a*b
    K = Group.generate([ab])
    vertices = G.left_cosets(K)
    assert len(vertices) == 12

    J = Group.generate([b])
    edges = G.left_cosets(J)
    assert len(edges) == 30

    from bruhat.solve import shortstr, zeros2, dot2
    from numpy import alltrue, zeros, dot

    def get_adj(left, right):
        A = zeros2((len(left), len(right)))
        for i, l in enumerate(left):
          for j, r in enumerate(right):
            lr = l.intersect(r)
            A[i, j] = len(lr)
        return A

    Hz = get_adj(faces, edges)
    Hxt = get_adj(edges, vertices)
    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hxt))
    #print()

    assert alltrue(dot2(Hz, Hxt)==0)

    # ------------ _oriented version ---------------
    
    # Bring's curve reflection group
    ngens = 3
    a, ai, b, bi, c, ci = range(2*ngens)
    rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
    ]
    a1 = (b,a)
    b1 = (a,c)
    rels += [ (3*a1+b1)*3 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 120 # == 12 * 10

    # --------------- Group ------------------
    G = graph.get_group()
    assert len(G) == 120

    a, ai, b, bi, c, ci = G.gens
    a1 = b*a
    b1 = a*c

    assert a1.order() == 5
    assert b1.order() == 2

    bc = b*c
    ba = b*a
    assert bc.order() == 5
    assert ba.order() == 5
    L = Group.generate([a1, b1])
    assert len(L) == 60, len(L)
    orients = G.left_cosets(L)
    assert len(orients) == 2
    #L, gL = orients
    orients = list(orients)

    H = Group.generate([a, b])
    assert len([g for g in H if g in L]) == 5
    faces = G.left_cosets(H)
    assert len(faces) == 12 # unoriented faces
    #hom = G.left_action(faces)
    o_faces = [face.intersect(orients[0]) for face in faces] # oriented faces
    r_faces = [face.intersect(orients[1]) for face in faces] # reverse oriented faces

    K = Group.generate([b, c])
    vertices = G.left_cosets(K)
    assert len(vertices) == 12 # unoriented vertices
    o_vertices = [vertex.intersect(orients[0]) for vertex in vertices] # oriented vertices

    J = Group.generate([a, c])
    u_edges = G.left_cosets(J)
    assert len(u_edges) == 30 # unoriented edges

    J = Group.generate([c])
    edges = G.left_cosets(J)
    assert len(edges) == 60, len(edges) # oriented edges ?

    # Here we choose an orientation on each edge:
    pairs = {}
    for e in u_edges:
        pairs[e] = []
        for e1 in edges:
            if e.intersect(e1):
                pairs[e].append(e1)
        assert len(pairs[e]) == 2

    def shortstr(A):
        s = str(A)
        s = s.replace("0", '.')
        return s

    Hz = zeros((len(o_faces), len(u_edges)), dtype=int)
    for i, l in enumerate(o_faces):
      for j, e in enumerate(u_edges):
        le = l.intersect(e)
        if not le:
            continue
        e1, e2 = pairs[e]
        if e1.intersect(le):
            Hz[i, j] = 1
        elif e2.intersect(le):
            Hz[i, j] = -1
        else:
            assert 0

    Hxt = zeros((len(u_edges), len(o_vertices)), dtype=int)
    for i, e in enumerate(u_edges):
      for j, v in enumerate(o_vertices):
        ev = e.intersect(v)
        if not ev:
            continue
        e1, e2 = pairs[e]
        if e1.intersect(ev):
            Hxt[i, j] = 1
        elif e2.intersect(ev):
            Hxt[i, j] = -1
        else:
            assert 0

    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hxt))
    #print()

    assert alltrue(dot(Hz, Hxt)==0)


def make_hyperbolic_group(idx):
    # start with hyperbolic Coxeter reflection group: a--5--b--5--c
    # Then we add another generator "d" that halves these triangles.
    ngens = 8
    a, ai, b, bi, c, ci, d, di = range(ngens)
    rels_552 = [
        (ai, a), (bi, b), (ci, c), (di, d),
        (a,)*2, (b,)*2, (c,)*2, (d,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
        (ci,d,a,d),
        (a,d)*4, (b,d)*2, (d,c)*4,
    ]

    order = argv.get("order")

    # Bring's curve & some random generators found from running make_random_55.
    rels = [
#        (3*(b,a)+(a,c))*3,                 # 120  [[30,8,3]]
        (a,b,c,b)*3,
        (0, 2, 0, 2, 2, 4, 2, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 4,
            0, 4, 0, 2, 2, 4, 0, 2, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2,
            0, 4, 0, 2),                   # 160  [[40,10,4]]
        (2, 4, 0, 2, 0, 2, 0, 2, 0, 2, 2, 4, 0, 2, 0, 4,
            0, 2, 2, 4, 2, 4, 2, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 4,
            0, 4, 0, 2, 0, 4, 0, 2, 2, 4, 0, 2, 0, 2, 0, 4, 0, 2,
            0, 4, 0, 4, 2, 4, 0, 4),       # 320  [[80,18,5]]
        (0, 2, 0, 2, 0, 4, 0, 2, 2, 4, 0, 2, 0, 4, 0, 2, 0, 4,
            2, 4, 2, 4, 0, 2, 0, 4, 2, 4, 0, 4, 0, 2, 2, 4, 0, 4,
            2, 4, 0, 4, 2, 4, 0, 2, 0, 4, 0, 2, 0, 4, 0, 4, 0, 4,
            2, 4, 0, 4, 2, 4),             # 600 [[150,32,6]]
        (0, 2, 0, 4, 0, 2, 0, 2, 0, 4, 0, 2, 2, 4, 0, 2, 2, 4,
            0, 2, 2, 4, 0, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2, 0, 2,
            2, 4, 0, 4, 2, 4, 0, 4, 2, 4, 0, 4, 0, 2, 0, 2, 2, 4,
            2, 4, 2, 4, 2, 4),             # 720  [[180,38,4]]
        (2, 4, 0, 4, 0, 4, 0, 4, 2, 4, 0, 2, 2, 4, 2, 4, 0, 2,
            0, 2, 2, 4, 0, 2, 2, 4, 0, 4, 0, 4, 0, 2, 2, 4, 0, 4,
            0, 2, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2, 2, 4, 2, 4, 2, 4,
            0, 2, 0, 2, 0, 2, 2, 4, 2, 4, 0, 2, 0, 2, 0, 2, 2, 4,
            0, 2, 2, 4, 0, 2, 0, 4, 0, 2, 2, 4, 0, 4, 0, 4, 2, 4,
            0, 2, 0, 2, 2, 4, 0, 4, 0, 2), # 1320  [[330,68,6]]
        (2, 4, 2, 4, 0, 2, 0, 2, 0, 2, 0, 2, 2, 4, 2, 4, 0, 2,
            2, 4, 0, 4, 0, 4, 0, 2, 0, 4, 0, 4, 0, 4, 0, 4, 2, 4,
            0, 4, 0, 2, 2, 4, 0, 2, 0, 2, 0, 4, 0, 2, 2, 4, 2, 4,
            2, 4, 0, 4, 2, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 2, 0, 2,
            0, 2, 0, 4, 2, 4, 0, 4, 2, 4, 0, 4, 0, 4, 0, 4, 0, 4,
            0, 2, 0, 2, 0, 4, 0, 4, 0, 2), # 1920  [[480,98,5 or 6?]]
    ]

    # reduced words
    rels = [
        (a,b,c,b)*3,
        (0, 2, 0, 4, 2, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 4, 0, 2, 0, 4, 0, 2),
        (2, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 2, 0, 2, 0, 2, 0, 4),
        (0, 2, 0, 2, 0, 2, 0, 4, 0, 2, 0, 4, 2, 4, 2, 4, 0, 2,
          0, 4, 2, 4, 0, 4, 2, 4, 0, 4, 2, 4, 0, 2, 0, 4, 0, 2, 0, 4, 2, 4, 0, 4, 2, 4),
        (0, 2, 0, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 4, 2, 4, 2, 4, 2, 4),
        (2, 4, 0, 4, 2, 4, 0, 4, 2, 4, 0, 2, 0, 2, 0, 4, 2, 4,
          0, 2, 0, 2, 0, 4, 2, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2),
        (2, 4, 2, 4, 0, 2, 0, 2, 0, 2, 0, 4, 2, 4, 0, 2, 0, 4,
          2, 4, 0, 4, 2, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 2, 0, 2,
          0, 2, 0, 4, 2, 4, 0, 4, 2, 4, 0, 2, 0, 2, 0, 2),
    ]
    for rel in rels:
        assert len(rel)%2 == 0

    #for rel in rels:
    #    rel = reduce_word(rels_552, rel)
    #    print(rel)

    rel = rels[idx]
    graph = Schreier(ngens, rels_552 + [rel])
    graph.build()
    G = graph.get_group()

    return G


def test_cover():
    print("test_cover")

    G2 = make_hyperbolic_group(2)
    G1 = make_hyperbolic_group(1)

    assert len(G2) == 640
    assert len(G1) == 320

    hom = {}
    for g,h in zip(G2.gens, G1.gens):
        hom[g] = h

    hom = close_hom(hom)
    assert hom is not None

    assert len(hom) == len(G2)
    kern = [g for (g,h) in hom.items() if h.is_identity()]
    assert len(kern) * len(G1) == len(G2)

    swaps = [g for g in kern if not g.is_identity()]
    swap = swaps[0]
    assert (swap*swap).is_identity()

    s1 = Surface(G1)
    s2 = Surface(G2)

    # Check we have an unramified degree 2 covering
    #s1.faces # list of Coset's
    fwd = {}
    counts = {}
    for f2 in s2.faces:
        f1 = Coset([hom[g] for g in f2], G1.items)
        counts[f1] = counts.get(f1, 0) + 1
        fwd[f2] = f1
    
    print(list(counts.values()))
    
    for f2 in s2.faces:
        swapf2 = Coset([swap*g for g in f2], G2.items)
        assert(fwd[f2]==fwd[swapf2])


def make_codes_54():
    idx = argv.get("idx", 0)

    G = make_hyperbolic_group(idx)

    if argv.oriented:
        make_surface_54_oriented(G)
    else:
        make_surface_54(G)



class Surface(object):
    def __init__(self, G_0):
        a, ai, b, bi, c, ci, d, di = G_0.gens
    
        # 5,5 coxeter subgroup
        G_1 = Group.generate([a, b, c])
        assert G_0.is_subgroup(G_1)
    
        # orientation subgroup
        L_0 = Group.generate([a*b, a*d])
        assert G_0.is_subgroup(L_0)
        orients = G_0.left_cosets(L_0)
        assert len(orients) == 2
    
        H_1 = Group.generate([a, b])
        self.faces = G_1.left_cosets(H_1)
        self.act_faces = G_1.action_subgroup(H_1)

        print("faces:", len(self.faces))

        return
    
        K_1 = Group.generate([b, c])
        vertices = G_1.left_cosets(K_1)
    
        J_1 = Group.generate([a, c])
        edges = G_1.left_cosets(J_1)
    
        print("faces: %d, edges: %d, vertices: %d" % (
            len(faces), len(edges), len(vertices)))

        self.faces = faces
        self.edges = edges
        self.vertices = vertices
        self.act_faces = act_faces


def make_surface_54(G_0):
    "find Z/2 chain complex"

    from bruhat.solve import shortstr, zeros2, dot2
    from numpy import alltrue, zeros, dot

    print("|G_0| =", len(G_0))

    a, ai, b, bi, c, ci, d, di = G_0.gens

    # 5,5 coxeter subgroup
    G_1 = Group.generate([a, b, c])
    assert G_0.is_subgroup(G_1)

    # orientation subgroup
    L_0 = Group.generate([a*b, a*d])
    assert G_0.is_subgroup(L_0)
    orients = G_0.left_cosets(L_0)
    assert len(orients) == 2

    H_1 = Group.generate([a, b])
    faces = G_1.left_cosets(H_1)
    #face_vertices = G_0.left_cosets(H_1) # G_0 swaps faces & vertices

    K_1 = Group.generate([b, c])
    vertices = G_1.left_cosets(K_1)

    J_1 = Group.generate([a, c])
    edges_1 = G_1.left_cosets(J_1)

    print("faces: %d, edges: %d, vertices: %d" % (
        len(faces), len(edges_1), len(vertices)))

    J_0 = Group.generate([a, d])
    assert len(J_0) == 8

    edges_0 = G_0.left_cosets(J_0)
    print("edges_0:", len(edges_0))

    assert G_0.is_subgroup(J_0)
    act_edges_0 = G_0.action_subgroup(J_0)
    #print("here")
    #while 1:
    #    sleep(1)

    from bruhat.solve import shortstr, zeros2, dot2, array2, solve
    from numpy import alltrue, zeros, dot

    def get_adj(left, right):
        A = zeros2((len(left), len(right)))
        for i, l in enumerate(left):
          for j, r in enumerate(right):
            lr = l.intersect(r)
            A[i, j] = len(lr)>0
        return A

    Hz = get_adj(faces, edges_0)
    Hxt = get_adj(edges_0, vertices)
    Hx = Hxt.transpose()

    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hxt))
    #print()

    assert alltrue(dot2(Hz, Hxt)==0)

#    act_fold_faces = G_0.action_subgroup(H_1)

    dualities = []
    idx = 0
    for i, g in enumerate(G_0):
        assert g in G_0
        h_edge_0 = act_edges_0[g]
        n_edge = len(h_edge_0.fixed())

#        h_fold_faces = act_fold_faces[g]
#        n_fold_faces = len(h_fold_faces.fixed())
        count = 0
        #perms = [] # extract action on the fixed edges
        for e0 in edges_0:
            e1 = g*e0
            if e0 != e1:
                continue
            count += 1
        #    perm = e0.left_mul_perm(g)
        #    #print(perm, end=" ")
        #    perms.append(perm)
        assert count == n_edge

        if g in G_1:
            continue
        elif g.order() != 2:
            continue
        #elif g in L_0:
        #    continue
        #elif g.order() == 2 and g not in G_1 and g not in L_0:

        perm = h_edge_0.get_idxs()
        dualities.append(perm)
        #if g not in L_0:
        #    check_432(Hz, Hx, perm)
        #    return
        result = is_422_cat_selfdual(Hz, Hx, perm)

        print("%d: [|g|=%s,%s.%s.%s.%s]"%(
            len(dualities),
            g.order(),
            n_edge or ' ',  # num fixed edges
#            n_fold_faces or ' ',  # num fixed faces+vertexes
            ["    ", "pres"][int(g in G_1)], # preserves face/vertex distinction
            ["refl", "rot "][int(g in L_0)], # oriented or un-oriented
            ["        ","selfdual"][result],
        ), end=" ", flush=True)
        print()

    print()

    #dualities = Group.generate(dualities)
    #print(dualities)

    #check_dualities(Hz, Hxt, dualities)


def is_422_cat_selfdual(Hz, Hx, perm):
    "concatenate with the [[4,2,2]] code and see if we get a self-dual code"
    from numpy import alltrue, zeros, dot
    import qupy.ldpc.solve
    import bruhat.solve
    qupy.ldpc.solve.int_scalar = bruhat.solve.int_scalar
    from bruhat.solve import shortstrx, zeros2, dot2, array2, solve
    from qupy.ldpc.css import CSSCode
    Cout = CSSCode(Hz=Hz, Hx=Hx)
    #print(Cout)
    #Cin = CSSCode(Hz=array2([[1,1,1,1]]), Hx=array2([[1,1,1,1]]))
    #print(Cin)

    pairs = []
    singles = []
    for i in range(Cout.n):
        j = perm[i]
        if j < i:
            continue
        if i==j:
            singles.append(i)
        else:
            pairs.append((i, j))
    #print(singles, pairs)
    M = len(singles) + 4*len(pairs)

    # encoding matrices
    enc_z = zeros2(M, Cout.n)
    enc_x = zeros2(M, Cout.n)
    
    row = 0
    for col in singles:
        enc_z[row, col] = 1
        enc_x[row, col] = 1
        row += 1
    H = []
    for (i,j) in pairs:
        enc_z[row, i] = 1   # 1010
        enc_z[row+2, i] = 1
        enc_z[row, j] = 1   # 1100
        enc_z[row+1, j] = 1

        enc_x[row, i] = 1   # 1100
        enc_x[row+1, i] = 1
        enc_x[row, j] = 1   # 1010
        enc_x[row+2, j] = 1
        h = array2([0]*M)
        h[row:row+4] = 1
        H.append(h)

        row += 4
    assert row == M
    
    #print(shortstrx(enc_z, enc_x))

    Hz = dot2(enc_z, Cout.Hz.transpose()).transpose()
    Hx = dot2(enc_x, Cout.Hx.transpose()).transpose()
    assert alltrue(dot2(Hz, Hx.transpose()) == 0)

    Hz = numpy.concatenate((Hz, H))
    Hx = numpy.concatenate((Hx, H))
    assert alltrue(dot2(Hz, Hx.transpose()) == 0)

    C = CSSCode(Hz=Hz, Hx=Hx)
    assert C.k == Cout.k
    #print(C)

    lhs = (solve(Hz.transpose(), Hx.transpose()) is not None)
    rhs = (solve(Hx.transpose(), Hz.transpose()) is not None)
    return lhs and rhs



toric = None
def check_toric():
    global toric # arff !
    import qupy.ldpc.solve
    import bruhat.solve
    qupy.ldpc.solve.int_scalar = bruhat.solve.int_scalar
    from bruhat.solve import shortstr, zeros2, dot2, array2, solve
    from numpy import alltrue, zeros, dot
    l = argv.get("l", 2)
    from qupy.ldpc.toric import Toric2D
    toric = Toric2D(l)
    Hx, Hz = toric.Hx, toric.Hz
    assert alltrue(dot2(Hx, Hz.transpose()) == 0)

    from qupy.condmat.isomorph import Tanner, search
    src = Tanner.build2(Hx, Hz)
    #tgt = Tanner.build2(Hx, Hz)
    tgt = Tanner.build2(Hz, Hx) # weak duality

    mx, n = Hx.shape
    mz, n = Hz.shape

    fns = []
    perms = []
    for fn in search(src, tgt):
        assert len(fn) == mx+mz+n
        bitmap = []
        for i in range(n):
            bitmap.append( fn[i+mx+mz]-mx-mz )
        perm = tuple(bitmap)
        #print(bitmap)
        fixed = [i for i in range(n) if bitmap[i]==i]
        print("perm:", perm)
        print("fixed:", fixed)

        g = Perm(perm, list(range(n)))
        assert g.order() == 2

        perms.append(perm)

    for hx in Hx:
        print(toric.strop(hx))
        print("--->")
        hx = array2([hx[i] for i in perm])
        print(toric.strop(hx))
        print("--->")
        hx = array2([hx[i] for i in perm])
        print(toric.strop(hx))
        print()

    check_dualities(Hz, Hx.transpose(), perms)

def check_dualities(Hz, Hxt, dualities):
    from bruhat.solve import shortstr, zeros2, dot2, array2, solve, span
    from numpy import alltrue, zeros, dot

    import qupy.ldpc.solve
    import bruhat.solve
    qupy.ldpc.solve.int_scalar = bruhat.solve.int_scalar
    from qupy.ldpc.css import CSSCode
    Hz = Hz % 2
    Hx = Hxt.transpose() % 2
    code = CSSCode(Hz=Hz, Hx=Hx)
    print(code)
    n = code.n

    Lx, Lz = code.Lx, code.Lz

    # check we really do have weak dualities here:
    for perm in dualities:
        Hz1 = Hz[:, perm]
        Hxt1 = Hxt[perm, :]
        assert solve(Hxt, Hz1.transpose()) is not None
        assert solve(Hz1.transpose(), Hxt) is not None

        Lz1 = Lz[:, perm]
        Lx1 = Lx[:, perm]
        find_xz = solve(concatenate((Lx, Hx)).transpose(), Lz1.transpose()) is not None
        find_zx = solve(concatenate((Lz1, Hz1)).transpose(), Lx.transpose()) is not None
        #print(find_xz, find_zx)
        assert find_xz 
        assert find_zx 

        # the fixed points live simultaneously in the homology & the cohomology
        fixed = [i for i in range(n) if perm[i]==i]
        if len(fixed) == 0:
            continue
        v = array2([0]*n)
        v[fixed] = 1
        v.shape = (n, 1)
        find_xz = solve(concatenate((Lx, Hx)).transpose(), v) is not None
        find_zx = solve(concatenate((Lz, Hz)).transpose(), v) is not None
        #print(find_xz, find_zx)
        assert find_xz 
        assert find_zx 

    from qupy.ldpc.asymplectic import Stim as Clifford
    s_gates = []
    h_gates = []
    s_names = []
    for idx, swap in enumerate(dualities):

        fixed = [i for i in range(n) if swap[i] == i]
        print(idx, fixed)
        for signs in cross([(-1, 1)]*len(fixed)): # <------- does not scale !!! XXX
            v = [0]*n
            for i, sign in enumerate(signs):
                v[fixed[i]] = sign
            ux = numpy.dot(Hx, v)
            uz = numpy.dot(Hz, v)
            if numpy.alltrue(ux==0) and numpy.alltrue(uz==0):
                #print("*", end=" ")
                break
        #else:
        #    assert 0
        #print(v)
        #print()

        # transversal S/CZ
        g = Clifford.identity(n)
        name = []
        for i in range(n):
            j = swap[i]
            if j < i:
                continue
            if j==i:
                assert v[i] in [1, -1]
                if v[i] == 1:
                    op = Clifford.s_gate(n, i)
                    name.append("S_%d"%(i,))
                else:
                    op = Clifford.s_gate(n, i).inverse()
                    name.append("Si_%d"%(i,))
            else:
                op = Clifford.cz_gate(n, i, j)
                name.append("CZ_%d_%d"%(i,j))
            g = op*g
        name = "*".join(reversed(name))
        s_names.append(name)
        #print(g)
        #print()
        #assert g.is_transversal(code)

        s_gates.append(g)

        h = Clifford.identity(n)
        for i in range(n):
            h = h * Clifford.h_gate(n, i)

        for i in range(n):
            j = swap[i]
            if j <= i:
                continue
            h = h * Clifford.swap_gate(n, i, j)
        #print(g)
        #print()
        #assert h.is_transversal(code)
        h_gates.append(h)

    if 0:
        print()
        print("CZ:")
        CZ = Clifford.cz_gate(2, 0, 1)
        op = (1., [0, 0], [1,1])
        for i in range(4):
            print(op)
            op = CZ(*op)
        return

    for idx, sop in enumerate(s_gates):
        print("idx =", idx)
        #for hx in Hx:
        perm = dualities[idx]
        #for hx in span(Hx):
        for hx in Hx:
            #print("hx =", hx)
            #print(s_names[idx])
            phase, zop, xop = sop(1., None, hx)
            assert numpy.alltrue(xop==hx) # fixes x component
            print(phase, zop, xop, dot2(zop, xop))
            for (i, j) in enumerate(perm):
                if xop[i] and xop[j] and i < j:
                    print("pair", (i, j))
            if toric is None:
                continue
            print("xop =")
            print(toric.strop(xop))
            print("zop =")
            print(toric.strop(zop))
            print()
#            if phase != 1:
#                #for jdx, xx in enumerate(xop):
#                #    if xx:
#                #        print(jdx, end=" ")
#                #print("\nFAIL: phase =", phase)
#                #return
#                print("FAIL")
#                break
#            assert phase == 1, phase
#            print(phase)
#        else:
#            print("OK")




def make_surface_54_oriented(G0):
    " find integer chain complex "
    from bruhat.solve import shortstr, zeros2, dot2
    from numpy import alltrue, zeros, dot

    print()
    print("|G0| =", len(G0))

    a, ai, b, bi, c, ci, d, di = G0.gens

    # 5,5 coxeter subgroup
    G_55 = Group.generate([a, b, c])
    assert G0.is_subgroup(G_55)

    # orientation subgroup
    L = Group.generate([a*b, a*d])
    #assert len(L) == 60, len(L)
    print("L:", len(L))
    orients = G0.left_cosets(L)
    assert len(orients) == 2
    orients = list(orients)

    H = Group.generate([a, b])
    faces = G_55.left_cosets(H)
    fold_faces = G0.left_cosets(H)
    print("faces:", len(faces))
    print("fold_faces:", len(fold_faces))
    #hom = G.left_action(faces)
    o_faces = [face.intersect(orients[0]) for face in faces] # oriented faces
    r_faces = [face.intersect(orients[1]) for face in faces] # reverse oriented faces

    K = Group.generate([b, c])
    vertices = G_55.left_cosets(K)
    #assert len(vertices) == 12 # unoriented vertices
    o_vertices = [vertex.intersect(orients[0]) for vertex in vertices] # oriented vertices
    print("o_vertices:", len(o_vertices))

    J0 = Group.generate([a, d])
    assert len(J0) == 8

    assert G0.is_subgroup(J0)
    act_edges = G0.action_subgroup(J0)
    act_fold_faces = G0.action_subgroup(H)
    #print("stab:", len(act_edges.get_stabilizer()))
    print("|G0| =", len(G0))
    print("edges:", len(G0) // len(J0))
    edges = G0.left_cosets(J0)
    #for e in edges:
    #    k = choice(e.perms)
    #    e1 = e.left_mul(~k)
    #    assert set(e1.perms) == set(J0.perms)
    fold_perms = []
    for i, g in enumerate(G0):
        assert g in G0
        #h_edge = act_edges.tgt[act_edges.send_perms[i]]
        #h_fold_faces = act_fold_faces.tgt[act_fold_faces.send_perms[i]]
        h_edge = act_edges[g]
        h_fold_faces = act_fold_faces[g]
        n_edge = len(h_edge.fixed())
        n_fold_faces = len(h_fold_faces.fixed())
        count = 0
        #perms = [] # extract action on the fixed edges
        #for e0 in edges:
        #    #e1 = e0.left_mul(g)
        #    e1 = g*e0
        #    if e0 != e1:
        #        continue
        #    perm = e0.left_mul_perm(g)
        #    #print(perm, end=" ")
        #    perms.append(perm)
        #    count += 1
        #assert count == n_edge
        if g.order() == 2 and g not in G_55 and g not in L:
            fold_perms.append(g)
        else:
            continue
        #print([p.order() for p in perms] or '', end="")
        print("[|g|=%s,%s.%s.%s.%s]"%(
            g.order(),
            n_edge or ' ',  # num fixed edges
            n_fold_faces or ' ',  # num fixed faces+vertexes
            ["    ", "pres"][int(g in G_55)], # preserves face/vertex distinction
            ["refl", "rot "][int(g in L)], # oriented or un-oriented
        ), end=" ", flush=True)
        print()
    print()

    J = Group.generate([a, c])
    u_edges = G_55.left_cosets(J)
    print("u_edges:", len(u_edges))

    J = Group.generate([c])
    edges = G_55.left_cosets(J)
    print("edges:", len(edges))

    # Here we choose an orientation on each edge:
    pairs = {}
    for e in u_edges:
        pairs[e] = []
        for e1 in edges:
            if e.intersect(e1):
                pairs[e].append(e1)
        assert len(pairs[e]) == 2, len(pairs[e])

    def shortstr(A):
        s = str(A)
        s = s.replace("0", '.')
        return s

    Hz = zeros((len(o_faces), len(u_edges)), dtype=int)
    for i, l in enumerate(o_faces):
      for j, e in enumerate(u_edges):
        le = l.intersect(e)
        if not le:
            continue
        e1, e2 = pairs[e]
        if e1.intersect(le):
            Hz[i, j] = 1
        elif e2.intersect(le):
            Hz[i, j] = -1
        else:
            assert 0

    Hxt = zeros((len(u_edges), len(o_vertices)), dtype=int)
    for i, e in enumerate(u_edges):
      for j, v in enumerate(o_vertices):
        ev = e.intersect(v)
        if not ev:
            continue
        e1, e2 = pairs[e]
        if e1.intersect(ev):
            Hxt[i, j] = 1
        elif e2.intersect(ev):
            Hxt[i, j] = -1
        else:
            assert 0

    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hxt))
    #print()

    assert alltrue(dot(Hz, Hxt)==0)

    for perm in fold_perms:
        pass

    import qupy.ldpc.solve
    import bruhat.solve
    qupy.ldpc.solve.int_scalar = bruhat.solve.int_scalar
    from qupy.ldpc.css import CSSCode
    Hz = Hz % 2
    Hx = Hxt.transpose() % 2
    code = CSSCode(Hz=Hz, Hx=Hx)
    print(code)



if __name__ == "__main__":

    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    print("OK\n")

        
