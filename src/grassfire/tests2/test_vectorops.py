"""Tests for vectorops module"""

import math
from grassfire.vectorops import (
    sub, add, div, mul, angle, make_vector, unit, rotate90cw
)


def test():
    v1, v2 = map(float, range(0, 3)), map(float, range(4, 7))
    # v1 = [0,1,2]
    # v2 = [4,5,6]

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

    # test almost equal
#    print unit((1, 1)) == (0.5 * math.sqrt(2.), 0.5 * math.sqrt(2.))

#     assert distance((0, 0), (3, 4)) == 5
#     assert distance2((0, 0), (3, 4)) == 25


def test_bisector():
#    print "Bisector", bisector((1., 0.), (0., 1.))

    v1 = make_vector((0, 0), (1, -1))
    u1 = unit(rotate90cw(v1))

    v2 = make_vector((-1, -1), (0, 0))
    u2 = unit(rotate90cw(v2))

    # test almost equal
#    print bisector(u1, u2)
#    print (0., math.sqrt(2.))
#
    v1 = make_vector((0, 0), (0.01, -10))
    v2 = make_vector((-0.01, -10), (0, 0))

    u1 = unit(rotate90cw(v1))
    u2 = unit(rotate90cw(v2))

#    print u1
#    print u2
#    print bisector(u1, u2)

#    print "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format((0, 0),
#                                                            (0.01, -1000))
#    print "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format((-0.01, -1000),
#                                                            (0, 0))
#    print "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format((0, 0),
#                                                            bisector(u1, u2))
#     print "angle between vectors", angle(u1, u2), math.degrees(angle(u1, u2))
#     alpha = 0.5 * math.pi + 0.5 * angle(u1, u2)
#     print math.degrees(alpha)
#     print 1. / math.sin(alpha)
#
#     for i in range(1, 100, 10):
#         print (1, -i)
#         v1 = make_vector((0, 0), (1, -i))
#         u1 = unit(rotate90cw(v1))
#         v2 = make_vector((-1, -i), (0, 0))
#         u2 = unit(rotate90cw(v2))
#
#         print "u1,u2", u1, u2
#         print unit(add(u1, u2))
#         # test almost equal
#         print "bisector found", bisector(u1, u2)


if __name__ == "__main__":
    test()
    test_bisector()
