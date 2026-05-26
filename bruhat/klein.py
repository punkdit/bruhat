#!/usr/bin/env python3

"""
Klein Quartic

"""

import sys
from random import randint, seed, choice
from functools import reduce
from operator import mul

#from bruhat.gset import Perm, Group, Coset, mulclose # FAIL
#from bruhat.action import Perm, Group, Coset, mulclose, close_hom
from bruhat.action import mulclose_hom, mulclose_names
from bruhat.util import cross
from bruhat.smap import SMap
from bruhat.argv import argv

from bruhat.todd_coxeter import Schreier
from bruhat.disc import mktriangle, Disc, Geodesic, Mobius, mulclose

from huygens.namespace import *


from huygens import config
config(text="pdflatex", latex_header=r"""
\usepackage{amsmath}
\usepackage{amssymb}
""")


def make_gap():
    """
    gap> LoadPackage("recog");
    true
    gap> LoadPackage("LINS");
    gap> F := FreeGroup("a","b","c");;
    gap> AssignGeneratorVariables(F);;
    #I  Assigned the global variables [ a, b, c ]
    gap> G := F/[a^8,b^2,c^3,(a*c)^2,c*a*b];;
    gap> gr := LowIndexNormalSubgroupsSearchForAll(G, 1000);;
    gap> L := List(gr);;
    gap> Length(L);
    gap> StructureDescription(FactorGroup(G,Grp(L[7])));
    "((C4 x C2) : C4) : S3"
    """



I = Mobius()

def render(l, m, n, c_words, f_words, e_words, v_words, h_words=[]):

    # build the rotation group generators
    ga, gb = [g.todisc() for g in mktriangle(l, m, n)]
    gc = ~(ga*gb)

    assert ga.order() == 2*l
    assert gb.order() == 2*m
    assert gc.order() == n, c.order()
    assert gc*ga*gb == I

    #parse = lambda w : reduce(mul, [{"a":ga, "b":gb, "c":gc}[wi] for wi in w], I)
    gens = [ga, gb, gc]
    parse = lambda w : reduce(mul, [gens[wi] for wi in w], I)
    act = lambda w,z: parse(w)(z)
    iact = lambda w,z: (~parse(w))(z)

    disc = Disc(scale=10)

    z_face = 0j
    z_vert = gc.inner_fixed()

    gamma = Geodesic.construct(z_vert, (~ga)(z_vert))
    z_edge = gamma.z2 # midpoint
    g_face = gamma.get_refl()

    gamma = Geodesic.construct(z_face, z_vert)
    g_edge = gamma.get_refl()
    g_vert = Mobius.conjugate()

    z_wert = g_vert(z_vert)

    z0 = (z_face + 2*z_edge) / 3

    # ------------------------------------------

    a, b, c = 0, 1, 2
    w = (b, a, a, b)
    z1 = act(w, z_face)
    zs = [z1]
    for i in range(7):
        zs.append(ga(zs[-1]))
    for i in range(8):
        z1 = zs[i]
        z2 = zs[(i+1)%8]
        disc.show_geodesic(z1, z2, attrs=st_thin+[orange.alpha(0.5)]+st_round)

    st = st_round+[grey.alpha(0.4)]
    for word in c_words:
        disc.show_geodesic(act(word, z_face), act(word, z_vert), attrs=[grey.alpha(0.5)])
        disc.show_geodesic(act(word, z_edge), act(word, z_vert), attrs=[0.5*grey])

    for (words,z,cl) in [
        (f_words, z_face, green), 
        (e_words, z_edge, blue), 
        (v_words, z_vert, red)]:
        for w in words:
            g = parse(w)
            disc.show_point(g(z), [cl.alpha(1.0)], radius=0.02)

    scale = 0.1

    #for w in h_words:
    #    z = act(w, z0)
    #    disc.show_label(z, r"$\star$", scale)

    for w in h_words:
        za = act(w, z_face)
        zb = act(w, z_vert)
        zc = act(w, z_wert)
        disc.show_polygon([za, zb, zc], None, [orange.alpha(0.2)])

    disc.fini()
    return disc.cvs



def test_render():

    maxsize = 10000

    l, m, n = 7, 2, 3
    l, m, n = 8, 2, 3

    ngens = 3
    a, b, c = range(ngens)
    rels = [ (a,)*l, (b,)*m, (c,)*n, (a, c)*2, (c,a,b) ]
    #klein = [(c,a,a,a,a,a)*4] # Klein quartic
    bolza = [((c,c)+(a,)*3)*2]
    print(bolza)

    graph = Schreier(ngens, rels+bolza)
    graph.build(maxsize=maxsize)
    print("graph:", len(graph))
    G = graph.get_gset()
    assert len(graph) == len(G)
    #print(G.structure_description()) # GL(2,3) != Octahedral group
    #return
    ga, gb, gc = G.gens
    names = mulclose_names([ga, gb, gc], "abc")
    assert len(names) == len(G)
    lookup = {g:tuple("abc".index(char) for char in names[g]) for g in G}

    graph = Schreier(ngens, rels)
    graph.build(maxsize=maxsize)
    faces = Schreier(ngens, rels)
    faces.build([(a,)], maxsize=maxsize)
    edges = Schreier(ngens, rels)
    edges.build([(b,)], maxsize=maxsize)
    verts = Schreier(ngens, rels)
    verts.build([(c,)], maxsize=maxsize)

    d = 16
    get = lambda graph : [word for word in graph.search_words() if len(word) < d]
    c_words = get(graph)
    f_words = get(faces)
    e_words = get(edges)
    v_words = get(verts)

    Hs = [H for H in G.subgroups()]
    Hs.sort(key = len)
    for H in Hs:
        print(len(H), H.structure_description())
    print()

    cvs = Canvas()
    for H in Hs:
        #print(len(H), end=" ", flush=True)
        if len(H) != 16:
            continue

        h_words = [lookup[h] for h in H]
        print(h_words)
        h_words = graph.search_words(h_words)

        fg = render(l, m, n, c_words, f_words, e_words, v_words, h_words)
        cvs.append(fg)
        cvs.show_page()

    name = "images/klein.pdf"
    print("writePDFfile", name)
    cvs.writePDFfile(name)


def main_render(words, name=None):

    (l, m, n, maxsize) = (7,2,3,800)

    # build the rotation group generators
    ga, gb = [g.todisc() for g in mktriangle(l, m, n)]
    gc = ~(ga*gb)

    assert ga.order() == 14
    assert gb.order() == 4
    assert gc.order() == 3, c.order()
    assert gc*ga*gb == I

    parse = lambda w : reduce(mul, [{"a":ga, "b":gb, "c":gc}[wi] for wi in w], I)

    cvs = Canvas()
    cvs.append(Scale(2.))
    disc = Disc(cvs)

    z_face = 0j
    z_vert = (ga*gb).inner_fixed()

    gamma = Geodesic.construct(z_vert, (~ga)(z_vert))
    z_edge = gamma.z2 # midpoint
    g_face = gamma.get_refl()
    gamma = Geodesic.construct(z_face, z_vert)
    g_edge = gamma.get_refl()
    g_vert = Mobius.conjugate()

    #z0 = (z_face + z_vert + z_edge) / 3
    z0 = (z_face + 2*z_edge) / 3

    gens = [g_face, g_edge, g_vert]
    gens = gens + [~g for g in gens]
    G = mulclose(gens, verbose=True, maxsize=maxsize)

    faces, edges, verts = [], [], []
    for g in G:
        faces.append(g(z_face))
        edges.append(g(z_edge))
        verts.append(g(z_vert))

    for g in G:
        disc.show_geodesic(g(z_face), g(z_vert), attrs=st_round+[grey.alpha(0.1)])
    #for g in G:
    #    disc.show_geodesic(g(z_face), g(z_edge), attrs=st_round+[grey])
    for g in G:
        disc.show_geodesic(g(z_vert), g(z_edge), attrs=st_round)


    for [cl, zs] in ([green, faces], [blue, edges], [red, verts]):
    #for [cl, zs] in ([red, verts],):
        for z in zs:
            disc.show_point(z, [cl], radius=0.02)

    scale = 0.1
    #disc.show_label(z0, r"$\star$", scale)

    for word in words:
        g = parse(word)
        #disc.show_label(g(z0), r"$%s$"%(word or r"\star"), scale)
        g = ~g
        disc.show_point(g(z0), [black], radius=0.02)

    disc.fini()
    if name is None:
        name = ("poincare-rotation-%d%d%d"%(l,m,n))
    #disc.save(name)
    print("writePDFfile", name)
    disc.cvs.writePDFfile("i"+"mages/"+name+".pdf")


def main():

    ngens = 3
    a, b, c = range(ngens)
    rels = [ (a,)*7, (b,)*2, (c,)*3, (a, c)*2, (c,a,b) ]
    rels += [ (c, a, a, a, a, a)*4 ]
    graph = Schreier(ngens, rels)
    graph.build()
    assert len(graph) == 168, len(graph)

    for gen in graph.get_gens():
        p = [gen[i] for i in range(len(gen))]
        #print(p)

    words = []
    for w in graph.get_words():
        word = ''.join('abc'[i] for i in w)
        #if len(word) < 5:
        words.append(word)
        print(w, word)

    G = graph.get_gset()
    print(G)
    ga, gb, gc = G.gens
    names = mulclose_names([ga, gb, gc], "abc")
    assert len(names) == len(G)
    print(G.structure_description())
    #for g in G:
        #print(names[g])
    #return

    count = 0
    for H in G.subgroups():
        print("[%d]"%(len(H),), end='', flush=True)
        if len(H) == 8:
            print()
            words = [names[h] for h in H]
            main_render(words, "klein-G-%d"%count)
            count += 1
    print()

    
    

if __name__ == "__main__":

    from time import sleep, time
    start_time = time()
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

    else:
        test()
        test_coxeter()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))

        
