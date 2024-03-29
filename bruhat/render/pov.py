#!/usr/bin/env python3

"""
Duplicate OpenGL coordinate system...

See:
https://gamedev.stackexchange.com/questions/153078/what-can-i-do-with-the-4th-component-of-gl-position
"""

import sys
from math import sin, cos, pi, sqrt

import numpy

scalar = numpy.float64
EPSILON = 1e-6

class Mat(object):
    def __init__(self, cs):
        A = numpy.array(cs, dtype=scalar)
        if len(A.shape)==1:
            m = len(A)
            A.shape = (m, 1) # col vector
        assert len(A.shape)==2, A.shape
        self.A = A
        self.shape = A.shape

    def strvec(self):
        v = self.A[:, 0]
        s = str(list(v))
        return "Mat(%s)"%(s,)

    def __str__(self):
        if self.shape[1] == 1:
            return self.strvec()
        A = self.A
        rows = [', '.join(["%.6f"%x for x in row]) for row in A]
        rows = ["[%s]"%row for row in rows]
        rows = "[%s]"%("\n".join(rows),)
        rows = rows.replace(".000000", ".      ")
        return rows
    __repr__ = __str__

    def __eq__(self, other):
        other = Mat.promote(other)
        assert self.shape == other.shape
        err = numpy.abs(self.A - other.A).sum()
        return err < EPSILON

    def __ne__(self, other):
        return not self.__eq__(other)

    def copy(self):
        A = self.A.copy()
        return Mat(A)

    @classmethod
    def promote(cls, item):
        if isinstance(item, Mat):
            return item
        m = Mat(item)
        return m

    @classmethod
    def identity(cls, n):
        A = numpy.identity(n)
        return cls(A)

    def __add__(self, other):
        other = Mat.promote(other)
        assert self.shape == other.shape
        A = self.A + other.A
        return Mat(A)

    def __sub__(self, other):
        other = Mat.promote(other)
        assert self.shape == other.shape
        A = self.A - other.A
        return Mat(A)

    def __neg__(self):
        A = -self.A
        return Mat(A)

    def __mul__(self, other):
        other = Mat.promote(other)
        assert self.shape[1] == other.shape[0]
        A = numpy.dot(self.A, other.A)
        return Mat(A)

    def __rmul__(self, r):
        A = r*self.A
        return Mat(A)

    def __getitem__(self, idx):
        if type(idx) is tuple:
            return self.A[idx] # <------ return
        elif type(idx) is slice:
            A = self.A[idx]
            return Mat(A) # <----- return
        if self.shape[1] == 1:
            idx = (idx, 0)
        return self.A[idx]

    def __setitem__(self, idx, value):
        if type(idx) is tuple:
            pass
        elif self.shape[1] == 1:
            idx = (idx, 0)
        self.A[idx] = value

    @classmethod
    def frustum(cls, left, right, bottom, top, nearval, farval):
        # mesa/src/mesa/math/m_matrix.c
        """
       GLfloat x, y, a, b, c, d;
       GLfloat m[16];
    
       x = (2.0F*nearval) / (right-left);
       y = (2.0F*nearval) / (top-bottom);
       a = (right+left) / (right-left);
       b = (top+bottom) / (top-bottom);
       c = -(farval+nearval) / ( farval-nearval);
       d = -(2.0F*farval*nearval) / (farval-nearval);  /* error? */
    
    #define M(row,col)  m[col*4+row]
       M(0,0) = x;     M(0,1) = 0.0F;  M(0,2) = a;      M(0,3) = 0.0F;
       M(1,0) = 0.0F;  M(1,1) = y;     M(1,2) = b;      M(1,3) = 0.0F;
       M(2,0) = 0.0F;  M(2,1) = 0.0F;  M(2,2) = c;      M(2,3) = d;
       M(3,0) = 0.0F;  M(3,1) = 0.0F;  M(3,2) = -1.0F;  M(3,3) = 0.0F;
    #undef M
    
       matrix_multf( mat, m, MAT_FLAG_PERSPECTIVE );
        """
        pass # TODO

    @classmethod
    def rotate(cls, angle, x, y, z):
        # angle in degrees
        s = sin(angle * pi / 180.0)
        c = cos(angle * pi / 180.0)
        M = cls.identity(4)
        r = sqrt(x*x + y*y + z*z)
        if r < EPSILON:
            return
        x /= r
        y /= r
        z /= r
        xx = x * x
        yy = y * y
        zz = z * z
        xy = x * y
        yz = y * z
        zx = z * x
        xs = x * s
        ys = y * s
        zs = z * s
        one_c = 1.0 - c
        
        M[0,0] = (one_c * xx) + c
        M[0,1] = (one_c * xy) - zs
        M[0,2] = (one_c * zx) + ys
        M[0,3] = 0.0
        
        M[1,0] = (one_c * xy) + zs
        M[1,1] = (one_c * yy) + c
        M[1,2] = (one_c * yz) - xs
        M[1,3] = 0.0
        
        M[2,0] = (one_c * zx) - ys
        M[2,1] = (one_c * yz) + xs
        M[2,2] = (one_c * zz) + c
        M[2,3] = 0.0
        
        
        M[3,0] = 0.0
        M[3,1] = 0.0
        M[3,2] = 0.0
        M[3,3] = 1.0
        return M

    @classmethod
    def translate(cls, *args):
        "modelled after glTranslate"
        n = len(args)+1
        A = numpy.identity(n)
        for i,val in enumerate(args):
            A[i, n-1] = val
        M = cls(A)
        return M

    @classmethod
    def scale(cls, sx, sy, sz):
        "modelled after glScale"
        A = numpy.identity(4)
        A[0, 0] = sx
        A[1, 1] = sy
        A[2, 2] = sz
        M = cls(A)
        return M

    @classmethod
    def perspective(cls, fovy, aspect, z_near, z_far):
        "modelled after gluPerspective"
        theta = fovy / 2 * pi / 180
    
        delta_z = z_far - z_near
        sine = sin(theta)
        if (delta_z == 0) or (sine == 0) or (aspect == 0):
            return
        cotangent = cos(theta) / sine
    
        A = numpy.identity(4)
        A[0,0] = cotangent / aspect
        A[1,1] = cotangent
        A[2,2] = -(z_far + z_near) / delta_z
        #A[2,3] = -1
        A[3,2] = -1
        #A[3,2] = -2 * z_near * z_far / delta_z
        A[2,3] = -2 * z_near * z_far / delta_z
        A[3,3] = 0

        M = Mat(A)
        return M

    def norm(self):
        A = self.A
        r = (A*A).sum()**0.5
        return r

    def normalized(self):
        r = self.norm()
        assert r>EPSILON
        A = self.A / r
        return Mat(A)

    def cross(self, other):
        assert self.shape == (3, 1)
        assert other.shape == (3, 1)
        cs = [
            self[1]*other[2] - self[2]*other[1],
            self[2]*other[0] - self[0]*other[2],
            self[0]*other[1] - self[1]*other[0]]
        return Mat(cs)

    def dot(self, other):
        assert self.shape == (3, 1)
        assert other.shape == (3, 1)
        r = (self.A*other.A).sum()
        return r

    @classmethod
    def lookat(cls, eye, center, up):
        "modelled after gluLookAt"
        eye = cls.promote(eye)
        center = cls.promote(center)
        up = cls.promote(up)

        forward = center - eye
        forward = forward.normalized()
    
        side = forward.cross(up)
        side = side.normalized()
    
        up = side.cross(forward)
    
        #M = cls.identity(4)
        M = numpy.identity(4)
        M[0,0] = side[0]
        M[1,0] = side[1]
        M[2,0] = side[2]
    
        M[0,1] = up[0]
        M[1,1] = up[1]
        M[2,1] = up[2]
    
        M[0,2] = -forward[0]
        M[1,2] = -forward[1]
        M[2,2] = -forward[2]
        M = M.transpose()
        M = Mat(M)
    
        M1 = cls.translate(-eye[0], -eye[1], -eye[2])
        M = M*M1
        return M


# ----------------------------------------------------------------------

def test_perspective():
    width, height = 640, 480
    proj = Mat.identity(4)
    
    M = Mat.perspective(45., width/height, 0.1, 100.)
    proj = M * proj
    
    assert proj == Mat([
        [ 1.8106601,  0.,         0.,         0.,       ],
        [ 0.,         2.4142137,  0.,         0.,       ],
        [ 0.,         0.,        -1.002002,  -0.2002002,  ],
        [ 0.,         0.,        -1., 0.,       ]])

test_perspective()


# ----------------------------------------------------------------------


from bruhat.render.base import SCALE_CM_TO_POINT
from bruhat.render.front import *


def mkpath(pts, closepath=True):
    pts = [path.moveto(*pts[0])]+[path.lineto(*p) for p in pts[1:]]
    if closepath:
        pts.append(path.closepath())
    p = path.path(*pts)
    return p


class GItem(object):
    def __init__(self, verts, epsilon=1e-4):
        assert len(verts) >= 3
        v0 = verts[0]
        for v in verts[1:]:
            v0 = v0 + v
        center = (1./len(verts))*v0

        if epsilon is not None and len(verts)>1:
            # try to cover up the seams. 
            # does not look good with alpha blending
            verts = [p + epsilon*(p-center).normalized() for p in verts]
        self.verts = verts
        self.center = center

    def render(self, cvs):
        pass


class GPoly(GItem):
    def __init__(self, verts, fill=None, stroke=None, epsilon=1e-2):
        GItem.__init__(self, verts, epsilon)
        self.fill = fill
        self.stroke = stroke

        v0, v1, v2 = verts[:3]
        a = v1-v0
        b = v2-v0
        ab = a.cross(b)
        assert ab.norm() > EPSILON
        self.normal = ab.normalized()

    def render(self, view, cvs):
        GItem.render(self, cvs)
        #fill, stroke = view.illuminate(self)
        fill = self.fill
        stroke = self.stroke
        verts = [view.trafo_canvas(v) for v in self.verts]
        v = self.center
        n = self.normal
        if fill is not None:
            fill = view.illuminate(v, n, fill)
        if stroke is not None:
            stroke = view.illuminate(v, n, stroke)
        cvs.append(Polygon(verts, fill, stroke))
        

class GMesh(GItem):
    def __init__(self, verts, normals, fill, epsilon=1e-2):
        GItem.__init__(self, verts, epsilon)

        assert len(verts) >= 3
        assert len(verts) == len(normals)
        v0, v1, v2 = verts[:3]
        a = v1-v0
        b = v2-v0
        ab = a.cross(b)
        normal = ab.normalized()
        self.normals = normals
        for n in normals:
            assert normal.dot(n) > 0.
        self.fill = fill

    def render(self, view, cvs):
        GItem.render(self, cvs)
        verts = [view.trafo_canvas(v) for v in self.verts]
        fill = self.fill
        fills = [view.illuminate(v, n, fill) 
            for (v,n) in zip(self.verts, self.normals)]
        cvs.append(Polymesh(verts, fills))
        

#class GBall(GItem):
#    def __init__(self, point, radius):
#        GItem.__init__(self, [point])
#        self.radius = radius
#        
#    def render(self, cvs):
#        GItem.render(self, cvs)


class Light(object):
    def __init__(self, position, color):
        assert position.shape == (3, 1)
        self.position = position
        self.color = color


class View(object):
    def __init__(self, _width=640, _height=480):
        #global width, height, viewport, proj, model
        scale = 0.05 # XXX 1./SCALE_CM_TO_POINT
        width, height = scale*_width, scale*_height
        self.viewport = (0., 0., width, height)
        self.proj = Mat.identity(4) # Projection matrix
        self.model = Mat.identity(4) # ModelView matrix
        self.stack = []
        self.gitems = []
        self.lights = []
    
    def perspective(self):
        #global proj
        width, height = self.viewport[2:]
        M = Mat.perspective(45., width/height, 0.1, 100.)
        self.proj = M * self.proj
    
    def translate(self, x, y, z): # XXX use Mat
        #global model
        M = Mat.translate(x, y, z)
        self.model = self.model*M
    
    def lookat(self, eye, center, up):
        #global model
        M = Mat.lookat(eye, center, up)
        self.model = self.model*M

    def rotate(self, angle, x, y, z):
        M = Mat.rotate(angle, x, y, z)
        self.model = self.model*M

    def scale(self, sx, sy, sz):
        M = Mat.scale(sx, sy, sz)
        self.model = self.model*M

    def save(self):
        model = self.model
        self.stack.append(model.copy())

    def restore(self):
        if not self.stack:
            return
        self.model = self.stack.pop()

    # ------------------------------------------------
    
    def trafo_view(self, point):
        assert isinstance(point, Mat), type(point)
        assert point.shape == (3, 1), repr(point)
        x, y, z = point
        v = [x, y, z, 1.]
        v = self.model * v
        assert abs(v[3]-1.) < EPSILON, "model matrix should not do this.."
        v = v[:3]
        return v

    def trafo_view_distance(self, point):
        assert isinstance(point, Mat), type(point)
        assert point.shape == (3, 1), repr(point)
        x, y, z = point
        v = [x, y, z, 0.]
        v = self.model * v
        v = v[:3]
        return v

    def trafo_camera(self, point):
        assert point.shape == (3, 1)
        x, y, z = point
        v = self.proj * [x, y, z, 1.]
        x, y, z, w = v
        #return x/w, y/w, z/w
        return x, y, z

    def trafo_canvas(self, point):
        x, y, z = point
        v = self.proj * [x, y, z, 1.]
        x, y, z, w = v
        x, y = x/w, y/w
    
        x0, y0, width, height = self.viewport
        w2, h2 = width/2, height/2
        x = x0 + w2 + x*w2
        y = y0 + h2 + y*h2
        return x, y
    
    def get_depth(self, gitem):
        v = gitem.center
        x, y, z = self.trafo_camera(v)
        return -z

    # -----------------------------------------
    # class Scene ?

    def add_gitem(self, face):
        self.gitems.append(face)

    def add_poly(self, verts, *args, **kw):
        verts = [self.trafo_view(v) for v in verts]
        gitem = GPoly(verts, *args, **kw)
        self.add_gitem(gitem)

    def add_mesh(self, verts, normals, *args, **kw):
        verts = [self.trafo_view(v) for v in verts]
        normals = [self.trafo_view_distance(n) for n in normals]
        gitem = GMesh(verts, normals, *args, **kw)
        self.add_gitem(gitem)

    def add_ball(self, point, radius, *args, **kw):
        point = self.trafo_view(point)
        gitem = Ball(point, radius, *args, **kw)
        self.add_gitem(gitem)

    def add_light(self, position, color):
        position = Mat.promote(position)
        position = self.trafo_view(position)
        light = Light(position, color)
        self.lights.append(light)

    def prepare_canvas(self, bg=color.rgb.black, clip=True):
        cvs = canvas.canvas()
        cvs.append(style.linewidth.THick)
    
        x0, y0, width, height = self.viewport
        p = mkpath([(x0, y0), (x0+width, y0), (x0+width, y0+height), (x0, y0+height)])
        cvs.fill(p, [bg])
        if clip:
            cvs.clip(p)

        cvs.append(style.linejoin.bevel)
        return cvs

#    def XXXilluminate(self, gitem):
#        light = self.lights[0]
#        v = (light.position - gitem.center).normalized()
#        #v = Mat([0., 0., 1.]).normalized()
#        x = v.dot(gitem.normal)
#        #print(x)
#        assert x <= 1+EPSILON, x
#        x = max(0.3, x)
#        fill = gitem.fill
#        stroke = gitem.stroke
#        if fill is not None:
#            r, g, b, a = fill
#            fill = (x*r, x*g, x*b, a)
#        if stroke is not None:
#            r, g, b, a = stroke
#            stroke = (x*r, x*g, x*b, a)
#        return fill, stroke

    def illuminate(self, vert, normal, color):
        light = self.lights[0]
        v = (light.position - vert).normalized()
        x = v.dot(normal)
        assert x <= 1+EPSILON, x
        x = max(0.3, x) # XXX diffuse
        r, g, b, a = color
        color = (x*r, x*g, x*b, a)
        return color

    def render(self, *args, **kw):
        cvs = self.prepare_canvas(*args, **kw)

        gitems = list(self.gitems)

        # XXX sorting by depth does not always work...
        # XXX try subdividing your GItem's ?
        gitems.sort(key = self.get_depth)

        for gitem in gitems:
            gitem.render(self, cvs)

        return cvs


# ----------------------------------------------------------------------


def make_sphere(view, radius, slices=8, stacks=8, fill=color.rgb.white):
    z = -radius
    dz = 2*radius / stacks
    dtheta = 2*pi/slices
    dphi = pi/stacks

    for i in range(stacks):
        phi0 = dphi*i
        phi1 = dphi*(i+1)
        r0 = radius*sin(phi0)
        r1 = radius*sin(phi1)
        z0 = -radius*cos(phi0)
        z1 = -radius*cos(phi1)
        for j in range(slices):
            theta0 = j * dtheta
            theta1 = (j+1) * dtheta
            x0, y0 = r0*cos(theta0), r0*sin(theta0)
            x1, y1 = r1*cos(theta0), r1*sin(theta0)
            x2, y2 = r1*cos(theta1), r1*sin(theta1)
            x3, y3 = r0*cos(theta1), r0*sin(theta1)
            verts = [
                Mat([x3, y3, z0]),
                Mat([x2, y2, z1]),
                Mat([x1, y1, z1]),
                Mat([x0, y0, z0]),
            ]
            if i==0:
                verts.pop(0)
            elif i==stacks-1:
                verts.pop(1)
            normals = [v.normalized() for v in verts]
            view.add_mesh(verts, normals, fill)


def make_cylinder(view, radius0, radius1, height, slices=8, fill=color.rgb.white):

    assert radius0 > EPSILON
    assert radius1 > EPSILON
    assert height > EPSILON

    dtheta = 2*pi/slices

    r0 = radius0
    r1 = radius1
    z = Mat([0, 0, height])
    for j in range(slices):
        theta0 = j * dtheta
        theta1 = (j+1) * dtheta
        v0 = Mat([cos(theta0), sin(theta0), 0.])
        v1 = Mat([cos(theta1), sin(theta1), 0.])
        dv0 = Mat([-sin(theta0), cos(theta0), 0.])
        dv1 = Mat([-sin(theta1), cos(theta1), 0.])

        verts = [r0*v1, r1*v1 + z, r1*v0 + z, r0*v0]
        n0 = dv1.cross(verts[1]-verts[0]).normalized()
        n1 = n0
        n2 = dv0.cross(verts[2]-verts[3]).normalized()
        n3 = n2
        normals = [n0, n1, n2, n3]
        view.add_mesh(verts, normals, fill)


def make_cone(view, radius, height, slices=8, fill=color.rgb.white):

    assert radius > EPSILON
    assert height > EPSILON

    dtheta = 2*pi/slices

    r = radius
    z = Mat([0, 0, height])
    for j in range(slices):
        theta0 = j * dtheta
        theta1 = (j+1) * dtheta
        v0 = Mat([cos(theta0), sin(theta0), 0.])
        v1 = Mat([cos(theta1), sin(theta1), 0.])
        dv0 = Mat([-sin(theta0), cos(theta0), 0.])
        dv1 = Mat([-sin(theta1), cos(theta1), 0.])

        verts = [r*v1, z, r*v0]
        n0 = dv1.cross(verts[1]-verts[0]).normalized()
        n2 = dv0.cross(verts[1]-verts[2]).normalized()
        n1 = (0.5*(n0+n2)).normalized()
        normals = [n0, n1, n2]
        view.add_mesh(verts, normals, fill)


def make_torus(view, inner, outer, slices=16, stacks=16, fill=color.rgb.white):

    dphi = 2*pi/stacks
    dtheta = 2*pi/slices

    for i in range(stacks):
        phi0 = dphi*i
        phi1 = dphi*(i+1)

        u0 = Mat([cos(phi0), sin(phi0), 0.])
        v0 = Mat([0., 0., 1.])
        u1 = Mat([cos(phi1), sin(phi1), 0.])
        v1 = Mat([0., 0., 1.])

        for j in range(slices):
            theta0 = dtheta*j
            theta1 = dtheta*(j+1)

            n0 = sin(theta0)*u0 + cos(theta0)*v0
            n1 = sin(theta0)*u1 + cos(theta0)*v1
            n2 = sin(theta1)*u1 + cos(theta1)*v1
            n3 = sin(theta1)*u0 + cos(theta1)*v0
            x0 = outer*u0 + inner*n0
            x1 = outer*u1 + inner*n1
            x2 = outer*u1 + inner*n2
            x3 = outer*u0 + inner*n3
            verts = [x3, x2, x1, x0]
            normals = [n3, n2, n1, n0]
            view.add_mesh(verts, normals, fill)


# ----------------------------------------------------------------------


def main():

    global cvs

    from bruhat import platonic
    polytopes = [
        platonic.make_tetrahedron(),
        platonic.make_cube(),
        platonic.make_octahedron(),
        platonic.make_dodecahedron(),
        platonic.make_icosahedron()]

    polytopes = [
        [[Mat(list(v)) for v in face] for face in polygon]
        for polygon in polytopes]


    from bruhat.argv import argv
    frames = argv.get("frames", 1)

    R = 6.0
    y = 2.
    theta = 0.
    for frame in range(frames):

        view = View()
        view.perspective()
        theta += 0.004*pi
        x, z = R*sin(theta), R*cos(theta)
        #R -= 0.01
        #y += 0.01
        view.lookat([x, y, z], [0., 0, 0], [0, 1, 0]) # eye, center, up
        #view.lookat([x, 0., z], [0., 0, 0], [0, 1, 0]) # eye, center, up

        point = [x, y, z]
        view.add_light(point, (1., 1., 1., 1.))

        #stroke = (0.4, 0.4, 0.4, 1.)
        stroke = None

        fills = [
            (0.9, 0.8, 0., 1.0),
            (0.9, 0.8, 0., 0.8),
            (0.6, 0.2, 0.4, 0.8),
            (0.2, 0.2, 0.4, 0.8),
            (0.8, 0.6, 0.4, 0.8),
            (0.0, 0.6, 0.4, 0.8),
        ]

        if 1:
            #view.rotate(frame, 1, 0, 0)
            #view.rotate(0.1*frame, 0, 1, 0)
            fill = fills[0]
            #make_torus(view, 0.5, 2., 32, 32, fill)
            #view.scale(2, 1, 1)
            #make_sphere(view, 1., 16, 12, fill)
            view.translate(0, -1, 0)
            view.rotate(-90, 1, 0, 0)
            #make_cylinder(view, 0.5, 1.0, 2., 16, fill)
            make_cone(view, 0.5, 1., 16, fill)

        elif 1:
            view.translate(-4., 0, 2)
    
            for idx, polygon in enumerate(polytopes):
                fill = fills[idx]
    
                view.save()
                view.rotate(-frame*(idx+1), 1, 1, 0)
                for verts in polygon:
                    view.add_poly(verts, fill, stroke)
                view.restore()
    
                view.translate(+5, 0, 0)
        
                view.rotate(360./5, 0, 1, 0)
    
        bg = color.rgb(0.2, 0.2, 0.2, 1.0)
        cvs = view.render(bg=bg)
        cvs.writePNGfile("frames/%.4d.png"%frame)
        if frame == 0:
            cvs.writePDFfile("frames/%.4d.pdf"%frame)

        print(".", end="", flush=True)


    print("OK")





if __name__ == "__main__":

    main()


