#!/usr/bin/env python
import sys, os
from math import *

from pyx import canvas, path, deco, trafo, style, text, color, deformer
from pyx.color import rgb

text.set(mode="latex") 
text.set(docopt="12pt")
text.preamble(r"\usepackage{amsmath,amsfonts,amssymb}")



north = [text.halign.boxcenter, text.valign.top]
northeast = [text.halign.boxright, text.valign.top]
northwest = [text.halign.boxleft, text.valign.top]
south = [text.halign.boxcenter, text.valign.bottom]
southeast = [text.halign.boxright, text.valign.bottom]
southwest = [text.halign.boxleft, text.valign.bottom]
east = [text.halign.boxright, text.valign.middle]
west = [text.halign.boxleft, text.valign.middle]
center = [text.halign.boxcenter, text.valign.middle]

black = rgb(0., 0., 0.) 
blue = rgb(0.1, 0.1, 0.9)
lred = rgb(0.8, 0.4, 0.2)
green = rgb(0.0, 0.6, 0.0)
white = rgb(1., 1., 1.) 
#shade = rgb(0.75, 0.55, 0)
grey = rgb(0.75, 0.75, 0.75)
yellow = rgb(1., 1., 0.)

st_dashed = [style.linestyle.dashed]
st_dotted = [style.linestyle.dotted]
st_round = [style.linecap.round]

st_thick = [style.linewidth.thick]
st_Thick = [style.linewidth.Thick]
st_THick = [style.linewidth.THick]
st_THICK = [style.linewidth.THICK]
#print dir(style.linewidth)


N = 100
d = {} # map (a, b) -> a^2 + b^2
for a in range(N):
  for b in range(a, N):
    k = a**2 + b**2
    d[k] = d.get(k, ()) + ((a, b),)


ks = []
for k in d.keys():
  if len(d[k])>1:
    #print k, d[k]
    ks.append(k)

ks.sort()
#for k in ks:
#    print k, d[k]


c = canvas.canvas()

m = 0.1
c.stroke(path.line(-m, 0., +m, 0.))
c.stroke(path.line(0., -m, 0., +m))

for k in ks:

    r = k**0.5
    #c.stroke(path.circle(0., 0., r))
    c.stroke(path.path(path.arc(0., 0., r, 45., 90)), [grey])


for k in ks:

    r = 0.1*len(d[k]) - 0.1
    for x, y in d[k]:
        c.fill(path.circle(x, y, r), [green])


c.writePDFfile("pic-rational.pdf")






