"""Tests for vectorops module"""

import math
from grassfire.vectorops import (
    sub, add, div, mul, angle, make_vector, unit, rotate90cw
)


def test():
    v1, v2 = map(float, range(0, 3)), map(float, range(4, 7))

    assert sub(v2, v1) == (4, 4, 4), sub(v2, v1)
    assert sub(v1, v2) == (-4, -4, -4), sub(v1, v2)
    assert sub(v1, 5) == (-5, -4, -3), sub(v1, 5)

    assert add(v1, v2) == (4, 6, 8)
    assert add(v1, 10) == (10, 11, 12)

    assert div((4, 4, 4), 4) == (1, 1, 1)
    assert div(v1, v2) == (0, 1. / 5., 2. / 6.), div(v1, v2)

    assert mul(v1, v2) == (0, 1 * 5, 2 * 6), mul(v1, v2)
    assert mul(v2, 10) == (40, 50, 60)

    assert angle((1, 1), (-1, 1)) == math.pi * 0.5


def test_bisector():

    v1 = make_vector((0, 0), (1, -1))
    u1 = unit(rotate90cw(v1))

    v2 = make_vector((-1, -1), (0, 0))
    u2 = unit(rotate90cw(v2))
    v1 = make_vector((0, 0), (0.01, -10))
    v2 = make_vector((-0.01, -10), (0, 0))

    u1 = unit(rotate90cw(v1))
    u2 = unit(rotate90cw(v2))


if __name__ == "__main__":
    test()
    test_bisector()
