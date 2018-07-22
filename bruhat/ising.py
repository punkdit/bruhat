#!/usr/bin/env python3

import sys, os

from math import exp
import random

import numpy
scalar = numpy.float64
bit = numpy.int32


qrnd = random.random


class IsingSquare(object):

    def __init__(self, n, T):
        #state = numpy.zeros((n, n), dtype=bit)
        state = numpy.random.multinomial(1, [0.5, 0.5], (n, n))[:, :, 0]
        state = 2*state - 1
        state[:] = -1
        #print state
        self.n = n
        self.state = state
        self.T = T # tempurature

    def update(self):
        n = self.n
        state = self.state
        beta = 1./self.T
        for count in range(100):
            i = random.randint(0, n-1)
            j = random.randint(0, n-1)
            b = state[(i-1)%n, j%n] +\
                state[(i+1)%n, j%n] +\
                state[i%n, (j-1)%n] +\
                state[i%n, (j+1)%n]
            p = 1. / (1+exp(-2*beta*b))
            if random.random() <= p:
                state[i, j] = +1
            else:
                state[i, j] = -1
#        state[n/2:n/2+5, n/2:n/2+5] = +1
#        state[n/2, 0:n] = +1

    def render(self, ctx, width, height):
        self.update()
        n = self.n
        state = self.state
        dw = 1.*width/n
        dh = 1.*height/n

        ctx.set_source_rgb(1, 1, 1)
        ctx.rectangle(0, 0, width, height)
        ctx.fill()

        ctx.set_source_rgb(0, 0, 0)
        for i in range(n):
            for j in range(n):
                if state[i, j]>0:
                    ctx.rectangle(i*dw, j*dh, dw-1, dh-1)
                    ctx.fill()


class IsingTree(object):

    def __init__(self, w, T):
        self.w = w # tree width (diameter)
        self.T = T
        even = []
        odd = []
        n = 1
        for i in range(1, w):
#            _row = []
#            assert len(row)==i
#            for j in range(2**i):
            if i%2:
                even += range(n, n+2**i)
            else:
                odd += range(n, n+2**i)
            n += 2**i
        self.sequence = even + odd
        state = numpy.random.multinomial(1, [0.5, 0.5], (n,))[:, 0]
        state = 2*state - 1
        state[:] = -1
        state[0] = 0 # invalid index
        #print state
        self.state = state
        self.n = n
        self.count = 0
        self.h = 0.

    def update(self):
        n = self.n
        h = self.h
        state = self.state
        beta = 1./self.T
        for i in self.sequence:
            j = (i-1)//2 or 3-i
            if 2*i+1 < n:
                b = state[2*i+1] + state[2*i+2] + state[j]
            else:
                b = state[j]
            b += h
            p = 1. / (1+exp(-2*beta*b))
            #if random.random() <= p:
            if qrnd() <= p:
                state[i] = +1
            else:
                state[i] = -1
        self.count += 1

    def render(self, ctx, width, height):
        for i in range(1):
            self.update()
        n = self.n
        w = self.w
        state = self.state

        ctx.set_source_rgb(0, 0, 0)
        ctx.rectangle(0, 0, width, height)
        ctx.fill()

        sh = 40.
        dh = 1.*(height-sh)/(w-1)

        j0 = 1
        for i in range(1, w):
            dw = 1.*width/(2**i)
            for j in range(2**i):
#                print j0+j, state[j0+j],
                if state[j0+j]>0:
                    ctx.set_source_rgb(1, 1, 1)
                    ctx.rectangle(j*dw, (i-1)*dh, max(1, dw-1.), dh-1)
                    ctx.fill()
                #else:
                #    ctx.set_source_rgb(1, 1, 1)
                #    ctx.rectangle(j*dw, (i-1)*dh, max(1, dw-1.), dh)
                #    ctx.fill() # white outline ?
            j0 += 2**i
#        print
        r = self.state.sum()
        ctx.set_source_rgb(1, 1, 1)
        ctx.set_font_size(20)
        ctx.move_to(0, height-5.)
        ctx.show_text("order = %.4f"%(1.*r/(n-1)))
        ctx.move_to(0.7*width, height-5.)
        ctx.show_text("count = %d"%self.count)


def build_fractal(depth):
    # first work with the directed graph, so we don't double count.
    nodes = [[1], []]
    depth -= 1

    # we replace each edge with four edges and two nodes:
    #
    #  ---
    #   |
    #   v
    #  ---
    #
    # goes to:
    #
    #  -------
    #   |   |
    #   v   v
    #  --- ---
    #   |   |
    #   v   v
    #  -------
    #


    while depth:
        n = len(nodes)
        top = len(nodes)
        next = len(nodes)
        for i in range(n):
            links = nodes[i]
            m = len(links)
            newlinks = []
            for j in range(m):
                # i --> links[j]
                k = links[j]
                newlinks.append(next)
                newlinks.append(next+1)
                nodes.append([k])
                nodes.append([k])
                next += 2
            nodes[i] = newlinks
            assert next == len(nodes)
        depth -= 1

    #print nodes

    # Now we create backlinks for the undirected graph
    newnodes = [[] for i in range(len(nodes))]
    for i in range(len(nodes)):
        links = nodes[i]
        for j in links:
            newnodes[i].append(j)
            newnodes[j].append(i)
    return newnodes

#print build_fractal(2)
nodes = build_fractal(3)
assert nodes == [
    [4, 5, 6, 7], [8, 9, 10, 11], [8, 9, 4, 5], [10, 11, 6, 7], 
    [0, 2], [0, 2], [0, 3], [0, 3], [2, 1], [2, 1], [3, 1], [3, 1]]


class IsingFractal(object):
    def __init__(self, depth, T):

        dy = 1./(1+2.**(depth-1))
        rects = [(0., 0., 1., dy), (0., 1.-dy, 1., dy)]
        for i in range(2, depth+1):
            dx = 2.**(1-i)
            #print "i = ", i, "dx=", dx
            x = 0.
            nodes = build_fractal(i)
            for j in range(len(rects), len(nodes)):
                links = nodes[j]
                y = sum(rects[k][1] for k in links) / len(links)
                rects.append((x, y, dx, dy))
                x += dx
                if x >= 1-1e-8:
                    x = 0.
            #print zip(xs, ys)
        #print
        self.rects = rects
        #print rects

        nodes = build_fractal(depth)
        n = len(nodes)
        state = numpy.random.multinomial(1, [0.5, 0.5], (n,))[:, 0]
        state = 2*state - 1
        state[:] = -1

        self.count = 0
        self.T = T
        self.n = n
        self.state = state
        self.nodes = nodes

    def update(self):
        n = self.n
        nodes = self.nodes
        state = self.state
        beta = 1./self.T
        idxs = list(range(n))
        random.shuffle(idxs)
        for i in idxs:
            b = 0.
            for j in nodes[i]:
                b += state[j]
            p = 1. / (1+exp(-2*beta*b))
            if random.random() <= p:
                state[i] = +1
            else:
                state[i] = -1
        self.count += 1

    def render(self, ctx, width, height):
        self.update()
        #print self.state

        n = self.n
        state = self.state
        rects = self.rects

        ctx.set_source_rgb(0, 0, 0)
        ctx.rectangle(0, 0, width, height)
        ctx.fill()

        for i in range(n):
            if state[i]<0:
                ctx.set_source_rgb(*OFF)
            else:
                ctx.set_source_rgb(*ON)

            x, y, dx, dy = rects[i]
            ctx.rectangle(x*width, y*height, dx*width-1, dy*height)
            ctx.fill()


def mkcolour(r, g, b):
    return (r/255., g/255., b/255.)

ON = mkcolour(252, 233, 79) # yellow
ON = mkcolour(196, 160, 0) # brown
ON = mkcolour(207, 122, 83) # brown
ON = mkcolour(239, 41, 41) # red

OFF = mkcolour(138, 226, 52) # green
OFF = mkcolour(115, 210, 22) # green
OFF = mkcolour(78, 154, 6)
OFF = mkcolour(93, 148, 165)


def save():

    #ising = IsingTree(10, 0.5)
    ising = IsingFractal(9, 4.0)

    import cairo

    width, height = 600, 600
    surface = cairo.SVGSurface("ising.svg", width, height)
    ctx = cairo.Context(surface)
    #print ctx.get_antialias()
    ctx.set_antialias(cairo.ANTIALIAS_GRAY)
    for i in range(1000):
        ising.update()
    ising.T = 2.0
    for i in range(1000):
        ising.update()
    ising.render(ctx, width, height)

    return



def test():

    import render

    if "square" in sys.argv:
        ising = IsingSquare(64, 2.0)
        render.display(ising.render, 400, 400)
    
    elif "tree" in sys.argv:
        ising = IsingTree(5, 0.5)
        render.display(ising.render, 800, 250)
    
    elif "fractal" in sys.argv:
        ising = IsingFractal(7, 0.5)
        render.display(ising.render, 400, 400)




if __name__ == "__main__":
    #test()

    save()




