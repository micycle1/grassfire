
from operator import sub as _sub, mul as _mul, truediv as _div, add as _add
from math import fsum, sqrt
from grassfire.calc import near_zero

import logging

def dot(v1, v2):
    """Returns dot product of v1 and v2 """
    assert len(v1) == len(v2), 'Vector dimensions should be equal'
    return fsum(p * q for p, q in zip(v1, v2))


def add(a, b):
    """Add a vector b to a, or add a scalar"""
    if hasattr(b, '__iter__'):
        if len(a) != len(b):
            raise ValueError('Vector dimensions should be equal')
        return tuple(map(_add, a, b))
    else:
        return tuple(ai + b for ai in a)


def mul(a, b):
    """Multiply a vector either element-wise with another vector, or with a
    scalar."""
    if hasattr(b, '__iter__'):
        assert len(a) == len(b), 'Vector dimensions should be equal'
        return tuple(map(_mul, a, b))
    else:
        return tuple(ai * b for ai in a)


def sub(a, b):
    """Subtract a vector b from a, or subtract a scalar"""
    if hasattr(b, '__iter__'):
        assert len(a) == len(b), 'Vector dimensions should be equal'
        return tuple(map(_sub, a, b))
    else:
        return tuple(ai - b for ai in a)


def div(a, b):
    """Element-wise division with another vector, or with a scalar."""
    if hasattr(b, '__iter__'):
        assert len(a) == len(b), 'Vector dimensions should be equal'
        return tuple(map(_div, a, b))
    else:
        return tuple(ai / b for ai in a)


def norm2(v):
    """Returns the norm of v, *squared*."""
    return dot(v, v)


def norm(a):
    """L2 norm ('distance')"""
    return norm2(a) ** 0.5


def unit(v):
    """Returns the unit vector in the direction of v."""
    return div(v, norm(v))


def make_vector(end, start):
    """Creates a vector from the start to the end.
    Vector is made based on two points: end -(minus) start.
    """
    return sub(end, start)


def dist(start, end):
    """Distance between two positons"""
    return norm(make_vector(end, start))


def coefficients_from_points(p, q):
    """ given 2 points p and q, give the coefficients (a,b,c) for the line
    that runs from p to q: ax + by + c = 0

    when a = 0, then the line is horizontal
    when b = 0, then the line is vertical

    see also:
    https://www.mathcentre.ac.uk/resources/uploaded/mc-ty-strtlines-2009-1.pdf
    """
    (px, py) = p
    (qx, qy) = q
    if py == qy:
        a = 0.
        if qx > px:
            b = 1.
            c = -py
        elif qx == px: # q and p at same x, direction unknown
            b = 0.
            c = 0.
        else:
            b = -1.
            c = py
    elif qx == px:
        b = 0.
        if qy > py: # q above p
            a = -1.
            c = px
        elif qy == py: # q and p at same y, direction unknown
            a = 0.
            c = 0.
        else:
            a = 1.
            c = -px
    else:
        a = py - qy
        b = qx - px
        c = -px * a - py * b
    return tuple(map(float, (a, b, c)))


def coefficients_perpendicular_through_point(la, lb, px, py):
    a = -lb
    b = la
    c = lb * px - la * py
    return tuple(map(float, (a, b, c)))

def coefficients_bisector_of_lines(pa, pb, pc, qa, qb, qc):
    n1 = sqrt(pa * pa + pb * pb)
    n2 = sqrt(qa * qa + qb * qb)
    a = n2 * pa + n1 * qa
    b = n2 * pb + n1 * qb
    c = n2 * pc + n1 * qc
    if (a == 0 and b == 0):
        a = n2 * pa - n1 * qa
        b = n2 * pb - n1 * qb
        c = n2 * pc - n1 * qc
    return tuple(map(float, (a, b, c)))


class Line2:
    def __init__(self, w, b, normalize=True):
        w = tuple(map(float, w))
        b = float(b)
        self.w = w
        self.b = b
        if normalize:
            self._normalize()

    def __eq__(self, other):
        return self.w == other.w and self.b == other.b

    def _normalize(self):
        nrm = norm(self.w)
        self.b /= nrm
        self.w = unit(self.w)

    def __repr__(self):
        return "Line2(w={}, b={})".format(self.w, self.b)

    def translated(self, v):
        """ translated line in direction of vector *v*"""
        d = dot(self.w, v)
        return Line2(self.w, self.b - d, normalize=(d != 0.0))

    def perpendicular(self, through):
        """ perp line """
        coeff = coefficients_perpendicular_through_point(self.w[0], self.w[1],
                                                         through[0], through[1])
        return Line2(coeff[:2], coeff[2])

    def bisector(self, other):
        """ bisector defined by self and other Line2 """
        a, b, c = coefficients_bisector_of_lines(self.w[0], self.w[1], self.b,
                                                 other.w[0], other.w[1], other.b)
        return Line2([a, b], c)

    def signed_distance(self, pt):
        """ signed distance of line *self* to point *pt*

        + = left
        - = right
        0 = on
        """
        return dot(self.w, pt) + self.b

    @property
    def through(self):
        """
        Returns nD-point through which this hyperplane passes
        """
        return mul(self.w, -self.b)

    @classmethod
    def from_points(cls, start, end): # !end! - start
        coeff = coefficients_from_points(start, end)
        ln = cls(coeff[:2], coeff[2])
        assert end != start
        dist_end = ln.signed_distance(end)
        dist_start = ln.signed_distance(start)
        return ln

    def y_at_x(self, x):
        (a, b), c = self.w, self.b
        return (-a*x-c) / b;

    def ends(self):
        ccw = rotate90ccw(self.w)
        cw = rotate90cw(self.w)
        end = add(mul(cw, 1000.), self.through)
        start = add(mul(ccw, 1000.), self.through)
        return (start, end)

    def visualize(self):
        return "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format(*self.ends())

    def at_time(self, now):
        """ Return Line2 at specific time """
        if now == 0.0:
            return Line2(self.w, self.b, normalize=False)
        else:
            logging.debug(' (constructing new line at t= {})'.format(now))
            return self.translated(mul(self.w, now)) # at_time(0) / # at_time(1)
class WaveFront:
    """ A line that remembers start and end point from which it is constructed """
    def __init__(self, start, end, line=None):
        if line is None:
            self.line = Line2.from_points(start, end)
        else:
            assert isinstance(line, Line2) == True
            self.line = line
        self.start = tuple(map(float, start))
        self.end = tuple(map(float, end))

    def __str__(self):
        return "WF start: {0} end: {1} line: {2[0]}·x + {2[1]}·y + {3} = 0.".format(self.start, self.end, self.line.w, self.line.b)


class WaveFrontIntersector:
    """ """
    def __init__(self, wf_left, wf_right):
        self.left = wf_left
        self.right = wf_right
        logging.debug(self.left)
        logging.debug(self.right)

    def get_bisector(self):
        intersector = LineLineIntersector(self.left.line, self.right.line)
        res = intersector.intersection_type()
        if res == LineLineIntersectionResult.LINE:
            bi = add(mul(self.left.line.w, 0.5), mul(self.right.line.w, 0.5))
        elif res == LineLineIntersectionResult.POINT:
            left_translated = self.left.line.translated(self.left.line.w)
            right_translated = self.right.line.translated(self.right.line.w)

            intersector_inner = LineLineIntersector(left_translated, right_translated)
            inner_res = intersector_inner.intersection_type()
            assert inner_res == LineLineIntersectionResult.POINT
            bi = make_vector(end=intersector_inner.result, start=intersector.result)
        elif res == LineLineIntersectionResult.NO_INTERSECTION:
            added = add(self.left.line.w, self.right.line.w)
            bi = added
        else:
            raise RuntimeError(f"Unknown intersection type: {res}")
        magn = norm(bi)
        logging.debug("magnitude of bisector: {}".format(magn))
        return bi

    def get_intersection_at_t(self, t):
        intersector = LineLineIntersector(self.left.line.at_time(t), self.right.line.at_time(t))
        if intersector.intersection_type() == LineLineIntersectionResult.POINT:
            return intersector.result
        else:
            raise ValueError('parallel lines, can not compute point of intersection')


class LineLineIntersectionResult:
    NO_INTERSECTION = 0
    POINT = 1
    LINE = 2


class LineLineIntersector:
    def __init__(self, one, other):
        self.one = one
        self.other = other
        self.result = None

    def intersection_type(self):
        one, other = self.one, self.other
        assert len(one.w) == 2
        assert len(other.w) == 2
        (a1, b1), c1 = one.w, one.b
        (a2, b2), c2 = other.w, other.b
        denom = a1 * b2 - a2 * b1
        if near_zero(denom) == True:
            x1 = a1 * c2 - a2 * c1
            x2 = b1 * c2 - b2 * c1
            if near_zero(x1) == True and near_zero(x2) == True:
                self.result = self.one
                return LineLineIntersectionResult.LINE
            else:
                self.result = None
                return LineLineIntersectionResult.NO_INTERSECTION
        else:
            num1 = b1 * c2 - b2 * c1
            num2 = a2 * c1 - a1 * c2
            x, y, w = num1, num2, denom
            xw = x / w
            yw = y / w
            self.result = (xw, yw)
            return LineLineIntersectionResult.POINT


def rotate90ccw(v):
    """Rotate 2d vector 90 degrees counter clockwise

    (x, y) -> (-y, x)
    """
    return (-(v[1]), v[0])


def rotate90cw(v):
    """Rotate 2d vector 90 degrees clockwise

    (x, y) -> (y, -x)
    """
    return (v[1], -v[0])