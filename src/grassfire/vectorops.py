"""Operations that allow tuples/lists (or any type that implements __getitem__
and __iter__) to be used as vectors"""

import logging
import math
from operator import sub as _sub, mul as _mul, truediv as _div, add as _add
from grassfire.calc import near_zero


def sub(a, b):
    """Subtract a vector b from a, or subtract a scalar"""
    if hasattr(b, '__iter__'):
        if len(a) != len(b):
            raise ValueError('Vector dimensions should be equal')
        return tuple(map(_sub, a, b))
    else:
        return tuple(ai - b for ai in a)


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
        if len(a) != len(b):
            raise ValueError('Vector dimensions should be equal')
        return tuple(map(_mul, a, b))
    else:
        return tuple(ai * b for ai in a)


def div(a, b):
    """Element-wise division with another vector, or with a scalar."""
    if hasattr(b, '__iter__'):
        if len(a) != len(b):
            raise ValueError('Vector dimensions should be equal')
        return tuple(map(_div, a, b))
    else:
        return tuple(ai / b for ai in a)


def make_vector(end, start):
    """Creates a vector from the start to the end.

    Vector is made based on two points: end -(minus) start.
    """
    return sub(end, start)


def dot(v1, v2):
    """Returns dot product of v1 and v2 """
    if len(v1) != len(v2):
        raise ValueError('Vector dimensions should be equal')
    return sum(p * q for p, q in zip(v1, v2))


def norm2(v):
    """Returns the norm of v, *squared*."""
    return dot(v, v)


def norm(a):
    """L2 norm"""
    return math.sqrt(norm2(a))

def dist(start, end):
    """Distance between two positons"""
    return norm(make_vector(end, start))


def unit(v):
    """Returns the unit vector in the direction of v."""
    return div(v, norm(v))


def cross(a, b):
    """Cross product between a 3-vector or a 2-vector"""
    if len(a) != len(b):
        raise ValueError('Vector dimensions should be equal')
    if len(a) == 3:
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])
    elif len(a) == 2:
        return a[0] * b[1] - a[1] * b[0]
    else:
        raise ValueError('Vectors must be 2D or 3D')


def angle(v1, v2):
    """angle between 2 vectors"""
    return math.acos(dot(v1, v2) / (norm(v1) * norm(v2)))


def angle_unit(v1, v2):
    """angle between 2 *unit* vectors

    does not compute the norm(v1)*norm(v2), as it is assumed to be 1
    """
    d = dot(v1, v2)
    d = max(-1.0, min(1.0, d))
    acos_d = math.acos(d)

    return d, acos_d


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