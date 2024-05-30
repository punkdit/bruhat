#!/usr/bin/env python

from math import cos, sin, pi
import os

#from mpmath import kleinj, mpc # sympy
import numpy
from PIL import Image

from numba import jit

from bruhat.argv import argv

EPSILON = 1e-6


def povray_sphere(png_name, xcam=-5, ycam=3, zcam=2, radius=1.0,
        light_source = "light_source { <-30, 0, 30> White }",
        background = "background { color red 0.97 green 0.96 blue 0.98}",
    ):
    src = """
     #include "colors.inc"
      camera {
        location <%s, %s, %s>
        look_at 0
        angle 30
        right <1,0,0>*image_width/image_height
      }
      %s
      %s
      sphere {
        <0, 0, 0>, %s
        pigment {
        image_map {
        png "%s"
        map_type 1 // spherical
        interpolate 2
        } }
      }
    """ % (xcam, ycam, zcam, light_source, background, radius, png_name)
    return src+"\n"




@jit
def weirstrass(z, tau=1j, N=5):
    z = complex(z)
    u = 1 / (z**2)
    for i in range(-N, N+1):
      for j in range(-N, N+1):
        if i==j==0:
            continue
        lmda = i + tau*j
        #if abs(z-lmda) < EPSILON:
        #    assert 0, cstr(z)
        #assert abs(lmda) > EPSILON, (i,j)
        u += 1 / ((z-lmda)**2) - 1/(lmda**2)
    return u


def cstr(z):
    s = "%.4f+%.4f"%(z.real, z.imag)
    s = s.replace("+-", "-")
    if s.startswith("-"):
        pass
    else:
        s = " "+s
    return s


def render(N, tau, W=1):
    A = numpy.zeros((N, N), dtype=complex)
    for i in range(1, N):
      for j in range(1, N):
        z = i + j*1j + (0.098 + 0.12j)
        z = W*z/N
        u = weirstrass(z, tau, N=5)
        #u1 = weirstrass(z, N=10)
        #print(cstr(z), cstr(u), cstr(u1), cstr(weirstrass(z,20)))
        A[i, j] = u
      #print('.', end='', flush=True)
    #print()

    return make_image(A)


def make_image(A):
    N = len(A)
    A0 = A.real
    #A0 -= A0.min() - 0.1
    #A0 /= A0.max()

    A1 = A.imag
    #A1 -= A1.min() - 0.1
    #A1[A1>100] = 100
    #A1 /= A1.max()

    #print()
    #print(A0.min(), A0.max())
    #print(A1.min(), A1.max())
    A = numpy.zeros((N, N, 3), dtype=numpy.uint8)
    A[:,:,0] = A0
    A[:,:,1] = A1
    A[:,:,2] = 0.2*(A1 + A0)
    im = Image.fromarray(A)
    return im


def stereographic(N):
    # https://en.wikipedia.org/wiki/Stereographic_projection
    A = numpy.zeros((N, N), dtype=complex)
    for i in range(N):
        theta = 2*pi*i/(N-1)
        for j in range(N):
            phi = pi*j/(N-1)
            bot = 1-cos(phi)
            if abs(bot) < EPSILON:
                continue
            else:
                R = sin(phi) / (1-cos(phi))
            x, y = R*cos(theta), R*sin(theta)
            z = x + 1j*y
            A[j,i] = 100*z
    im = make_image(A)
    return im


def main():

    theta = 2*pi/3
    tau = cos(theta) + 1j*sin(theta)
    N = 128*8

    #im = render(N, 1j, 2)
    #im.save("weierstrass_gauss.png")

    #im = render(N, tau, 2)
    #im.save("weierstrass_eisenstein.png")

    im = stereographic(N)
    name = "stereographic"
    im.save(name+".png")

    s = povray_sphere(name+".png", background=" ")
    f = open(name+".pov", "w")
    print(s, file=f)
    f.close()

    #cmd = "povray %s +O%s -D +W%s +H%s +WT%s" % (pov_name, output_name, width, height, threads)
    cmd = "povray stereographic.pov -D +W1000 +H800, +Oriemann_sphere.png"
    print(cmd)
    os.system(cmd)


    name = "worldmap"
    s = povray_sphere(name+".png", background=" ")
    f = open(name+".pov", "w")
    print(s, file=f)
    f.close()

    cmd = "povray worldmap.pov -D +W1000 +H800, +Oworldmap_sphere.png"
    print(cmd)
    os.system(cmd)





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
        main()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))



