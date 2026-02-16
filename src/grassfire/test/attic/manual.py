"""Manual tests, for checking output manually in QGIS
"""
from grassfire import ToPointsAndSegments, calc_skel, init_event_list, event_loop
from grassfire.collapse import compute_collapse_time

from tri.delaunay import Edge, output_edges
import time
from qgis._core import QGis


if __name__ == "__main__":
    import logging
    import sys

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    root.addHandler(ch)


def test_poly():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(0, 0), (9, 0), (11, -.1), (11.1,0), (22,0), (14,10), (2,8), (0, 5), (0,0)]])
    skel = calc_skel(conv, output=True, pause=True)

def test_cocircular():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(0,1), (1,0), (2,0), (3,1), (3,2), (2,3), (1,3), (0,2), (0,1)]])
    calc_skel(conv, output=True, pause=True)

def test_cocircular1():
    ok = (3.8,0.8) # this works
    fail = (4,1) # substitute with this and we get a lot of simultaneous events!
    conv = ToPointsAndSegments()
    conv.add_polygon([[(0,1), (1,0), (3,0), ok, (4,3), (3,4), (1,4), (0,3), (0,1)]])
    calc_skel(conv, output=True, pause=True)

def test_diamant():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(-1,0), (0,-1), (1,0), (0,5), (-1,0)]])
    skel = calc_skel(conv, pause=True, output=True)
    assert len(skel.segments()) == 8

def test_diamantlike():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(-15,0), (-1,0), (0,-1), (1,0), (15,0), (0,15), (-15,0)]])
    skel = calc_skel(conv, pause=True, output=True)
    assert len(skel.segments()) == (7+6)

def test_parallellogram():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(-15,0), (0,0), (15,25), (0, 25), (-15,0)]])
    calc_skel(conv)

def test_simple_poly():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(0, 0), (22,0), (14,10), (2,8), (0, 6.5), (0,0)]])
    skel = calc_skel(conv, output=True, pause = True)
    assert len(skel.segments()) == 12


def test_single_line():
    conv = ToPointsAndSegments()
    conv.add_point((0, 0))
    conv.add_point((10, 0))
    conv.add_segment((0, 0), (10,0))
    calc_skel(conv)

def test_three_lines():
    conv = ToPointsAndSegments()
    conv.add_point((0, 0))
    conv.add_point((10, 0))
    conv.add_point((-2, 8))
    conv.add_point((-2, -8))
    conv.add_segment((0, 0), (10,0))
    conv.add_segment((0, 0), (-2,8))
    conv.add_segment((0, 0), (-2,-8))
    skel = calc_skel(conv, output=True, pause = True)

def test_arrow_four_lines():
    conv = ToPointsAndSegments()

    conv.add_point((0, 0))
    conv.add_point((10, 0))
    conv.add_point((12,0.5))
    conv.add_point((8,5))
    conv.add_point((8,-5))

    conv.add_segment((0,0), (10,0))
    conv.add_segment((8,5), (12,0.5))
    conv.add_segment((12,0.5), (8,-5))

    calc_skel(conv)


def test_single_point():
    conv = ToPointsAndSegments()
    conv.add_point((0, 0))
    skel = calc_skel(conv)


def test_triangle():
    conv = ToPointsAndSegments()
    conv.add_point((10,0))
    conv.add_point((-2,8))
    conv.add_point((-2,-8))
    conv.add_segment((10,0), (-2,8))
    conv.add_segment((-2,8), (-2,-8))
    conv.add_segment((-2,-8), (10,0))
    skel = calc_skel(conv)
    assert len(skel.segments()) == 6

def test_quad():
    conv = ToPointsAndSegments()
    conv.add_point((4,5))
    conv.add_point((-2,8))
    conv.add_point((-2,-8))
    conv.add_point((14,10))
    conv.add_segment((14,10), (-2,8))
    conv.add_segment((-2,8), (-2,-8))
    conv.add_segment((4,5), (14,10))
    conv.add_segment((-2,-8), (4,5))
    skel = calc_skel(conv)
    assert len(skel.segments()) == 9


def test_tri_intermediate_pt():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(-1, 0), (0,0), (1,0), (0,14), (-1,0)]])
    skel = calc_skel(conv, output=True)
    assert len(skel.segments()) == 8


def test_tri_2intermediate_pts():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(-1, 0), (-0.1,-1), (0.1,-1), (1,0), (0,14), (-1,0)]])
    skel = calc_skel(conv, pause=True, output=True)
    assert len(skel.segments()) == 8

def test_tri_intermediate_pt_sharp():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(-1, 0), (0, 13.5), (1,0), (0,14), (-1,0)]])
    skel = calc_skel(conv, output=True)
    assert len(skel.segments()) == 8

def test_pointy_star():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(-5,10), (-.1,0), (-5, -10), (0,-9), (5,-10), (.1,0), (5,10), (0,9), (-5,10)]])
    skel = calc_skel(conv, pause=True, output=True)
    assert len(skel.segments()) == 8

def test_2triangle_1side_collapse():
    conv = ToPointsAndSegments()
    conv.add_polygon([[(0.1,0), (10,0), (9, 10), (2,9.5),(0,0.1), (0.1, 0)]])
    skel = calc_skel(conv, pause=True, output=True)
    assert len(skel.segments()) == (7+5)

def test_two_lines_par():
    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((10,0))
    conv.add_point((12,1))
    conv.add_point((22,1))

    conv.add_segment((0,0), (10,0))
    conv.add_segment((12,1), (22,1))

    calc_skel(conv)

def test_circular():
    from math import pi, cos, sin
    ring = []
    pi2 = 2 * pi
    ct = 10
    alpha = pi2 / ct
    for i in range(ct):
        ring.append( (cos(i* alpha), sin(i* alpha)))
    ring.append(ring[0])
    print ring
    conv = ToPointsAndSegments()
    conv.add_polygon([ring])
    calc_skel(conv, output=True, pause=True)


def test_cross():
    ring = [(0,0), (10, 0), (10,-10), (15, -10), (15,0), (25,0), (25,5), (15,5), (15,15), (10,15), (10,5), (0,5), (0,0)]
    conv = ToPointsAndSegments()
    conv.add_polygon([ring])
    skel = calc_skel(conv, output=True, pause=True)


def test_rocket():
    """Contains zero triangle to flip ...
    """
    ring = [(0,0), (10, 0), (15,5), (10,9), (1,7), (6,4), (0,0)]
    conv = ToPointsAndSegments()
    conv.add_polygon([ring])
    skel = calc_skel(conv, output=True, pause=True)
    print "DONE"


def test_bottom_circle():
    from math import pi, cos, sin, degrees
    ring = []
    pi2 = 2 * pi
    ct = 8
    alpha = pi / ct
    print degrees(alpha)
    for i in range(ct+1):
        ring.append( (cos(pi+i*alpha), sin(pi+i*alpha)))
    ring.append(ring[0])
    print ring
    conv = ToPointsAndSegments()
    conv.add_polygon([ring])
    calc_skel(conv, output=True, pause=True)

def test_bottom_circle_top_square():
    from math import pi, cos, sin, degrees
    ring = []
    pi2 = 2 * pi
    ct = 6
    alpha = pi / ct
    print degrees(alpha)
    for i in range(ct+1):
        ring.append( (cos(pi+i*alpha), sin(pi+i*alpha)))
    ring.extend([(1, 10), (-1,10)])
    ring.append(ring[0])
    conv = ToPointsAndSegments()
    conv.add_polygon([ring])
    calc_skel(conv, pause=True, output=True)

def test_polyline():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((10,-1))
    conv.add_point((22,1))
    conv.add_point((30,-5))

    conv.add_segment((0,0), (10,-1))
    conv.add_segment((10,-1), (22,1))
    conv.add_segment((22,1), (30,-5))

    calc_skel(conv)

def test_1_segment():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((10,0))

    conv.add_segment((0,0), (10,0))

    calc_skel(conv)

def test_2_segments():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((10,0))
    conv.add_point((22,0))
    conv.add_point((30,0))

    conv.add_segment((0,0), (10,0))
    conv.add_segment((22,0), (30,0))

    calc_skel(conv)


def test_2_perp_segments():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((10,0))
    conv.add_point((12,2))
    conv.add_point((12,10))

    conv.add_segment((0,0), (10,0))
    conv.add_segment((12,2), (12,10))

    calc_skel(conv)

def test_45_deg_segments():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((10,0))
    conv.add_point((12,2))
    conv.add_point((14,4))

    conv.add_segment((0,0), (10,0))
    conv.add_segment((12,2), (14,4))

    calc_skel(conv)

def test_30_deg_segments():

    conv = ToPointsAndSegments()

    conv.add_point((0,5))
    conv.add_point((9,0.5))
    conv.add_point((12,2))
    conv.add_point((14,4))

    conv.add_segment((0,5), (9,0.5))
    conv.add_segment((12,2), (14,4))

    calc_skel(conv)


def test_4_segments():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((10,0))
    conv.add_point((22,0))
    conv.add_point((30,0))

    conv.add_point((16,-3))
    conv.add_point((16,-6))

    conv.add_point((16,2))
    conv.add_point((16,6))

    conv.add_segment((0,0), (10,0))
    conv.add_segment((22,0), (30,0))
    conv.add_segment((16,-3), (16,-6))
    conv.add_segment((16,2), (16,6))

    conv.add_segment((0,0), (16,-6))
    conv.add_segment((16,-6), (30,0))
    conv.add_segment((30,0), (16,6))
    conv.add_segment((16,6), (0,0))

    calc_skel(conv)

def test_cocircular_segments():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((1,1))

    conv.add_point((3,0))
    conv.add_point((2,1))

    conv.add_point((0,3))
    conv.add_point((1,2))

    conv.add_point((3,3))
    conv.add_point((2,2))

    conv.add_segment((0,0), (1,1))
    conv.add_segment((3,0), (2,1))
    conv.add_segment((0,3), (1,2))
    conv.add_segment((3,3), (2,2))

    calc_skel(conv)

def test_parallel_movement():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((1,0))
    conv.add_point((2,0))
    conv.add_point((3,0))

    conv.add_segment((0,0), (1,0))
    conv.add_segment((1,0), (2,0))
    conv.add_segment((2,0), (3,0))

    calc_skel(conv)


def test_crash_vertex():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((1,0))


    conv.add_point((0,2))
    conv.add_point((0.5,1.5))
    conv.add_point((1,2))

    conv.add_segment((0,0), (1,0))
    conv.add_segment((0,2), (0.5,1.5))
    conv.add_segment((1,2), (0.5,1.5))

    calc_skel(conv)


def test_2triangle():
    from math import sqrt
    conv = ToPointsAndSegments()
    polygon = [[(1,0),
                (1, -10), (11,-10),
                (11, 10), (1,10), (1,2), (1-(sqrt(3)),1), (1,0)]]
    conv.add_polygon(polygon)
    calc_skel(conv)


def test_2triangle_eq_sides():
    from math import sqrt
    conv = ToPointsAndSegments()
    polygon = [[(1,0),
                (-5, 0), (-4.25, -6.65),# (-1,-0.9),
                (1, -10), (11,-10),
                (11, 10), (1,10), (1,2), (1-(sqrt(3)),1), (1,0)]]
    conv.add_polygon(polygon)
    calc_skel(conv)

def test4_3_3():
    tri = KineticTriangle()

    a = KineticVertex()
    a.origin = (0,0)
    a.velocity = (-sqrt(2),sqrt(2))

    b = KineticVertex()
    b.origin = (1,0)
    b.velocity = (sqrt(2),sqrt(2))

    c = KineticVertex()
    c.origin = (0.5, 1.5)
    c.velocity = (0,-1)

    tri.vertices = [a, b, c]

    print tri
    Mv = tuple(map(sub, c.origin, a.origin))
    print Mv

    m =  map(sub, b.origin, a.origin)
    m = norm(m)
    print m
    n = perp(m)
    print n
    distance_v_e = dot(Mv, n)
    print "dist", distance_v_e
    s = c.velocity
    crash_time = distance_v_e / (1 - dot(s, n))
    print "time vertex crashes on edge:", crash_time
    coeff = area_collapse_time_coeff(a, b, c)
    print solve_quadratic(coeff[0], coeff[1], coeff[2])
    import matplotlib.pyplot as plt
    print "roots", numpy.roots(coeff)
    x = range(-40, 100)
    y = [quadratic(y, coeff[0], coeff[1], coeff[2]) for y in x]
    plt.plot(x,y)
    plt.show()

    with open("/tmp/ktris.wkt", "w") as fh:
        output_triangles([tri], fh)

def test_compute_0():
    tri = KineticTriangle()
    a = KineticVertex()
    a.origin = (0,0)
    a.velocity = (sqrt(2),sqrt(2))
    b = KineticVertex()
    b.origin = (2,0)
    b.velocity = (-sqrt(2),sqrt(2))
    c = KineticVertex()
    c.origin = (1, 3)
    c.velocity = (0,-1)
    tri.neighbours = [True, True, True] # make them not None for the test
    tri.vertices = [a, b, c]
    print compute_collapse_time(tri)
    with open("/tmp/ktris.wkt", "w") as fh:
        output_triangles([tri], fh)

def test_compute_1():
    tri = KineticTriangle()
    a = KineticVertex()
    a.origin = (0,0)
    a.velocity = (sqrt(2),sqrt(2))
    b = KineticVertex()
    b.origin = (2,0)
    b.velocity = (-sqrt(2),sqrt(2))
    c = KineticVertex()
    c.origin = (1, 3)
    c.velocity = (0,-1)
    tri.neighbours = [True, True, None] # make them not None for the test
    tri.vertices = [a, b, c]
    print compute_collapse_time(tri)
    with open("/tmp/ktris.wkt", "w") as fh:
        output_triangles([tri], fh)

def test_compute_2():
    tri = KineticTriangle()
    a = KineticVertex()
    a.origin = (0,0)
    a.velocity = (sqrt(2),sqrt(2))
    b = KineticVertex()
    b.origin = (2,0)
    b.velocity = (-sqrt(2),sqrt(2))
    c = KineticVertex()
    c.origin = (1, 3)
    c.velocity = (0,-1)
    tri.neighbours = [None, True, None] # make them None for the test
    tri.vertices = [a, b, c]
    print compute_collapse_time(tri)
    with open("/tmp/ktris.wkt", "w") as fh:
        output_triangles([tri], fh)

def test_compute_3():
    tri = KineticTriangle()
    a = KineticVertex()
    a.origin = (0,0)
    a.velocity = (sqrt(2),sqrt(2))
    b = KineticVertex()
    b.origin = (1,0)
    b.velocity = (-sqrt(2),sqrt(2))
    c = KineticVertex()
    c.origin = (0.5, 1)
    c.velocity = (0,-1)
    tri.neighbours = [None, None, None] # make them None for the test
    tri.vertices = [a, b, c]
    print compute_collapse_time(tri)
    with open("/tmp/ktris.wkt", "w") as fh:
        output_triangles([tri], fh)


def test_flip():
    """Flip 2 triangles
    """
    tri0 = KineticTriangle()
    tri1 = KineticTriangle()
    tri2 = KineticTriangle()
    tri2.vertices = [None, "a", "c"]
    tri3 = KineticTriangle()
    tri3.vertices = [None, "b", "a"]
    tri4 = KineticTriangle()
    tri4.vertices = [None, "d", "b"]
    tri5 = KineticTriangle()
    tri5.vertices = [None, "c", "d"]

    tri2.neighbours[0] = tri0
    tri3.neighbours[0] = tri0
    tri4.neighbours[0] = tri1
    tri5.neighbours[0] = tri1

    tri0.vertices = ["a","b","c"]
    tri1.vertices = ["c","b","d"]

    tri0.neighbours = [tri1, tri2, tri3]
    tri1.neighbours = [tri4, tri5, tri0]

    assert tri1 in tri0.neighbours
    assert tri2 in tri0.neighbours
    assert tri3 in tri0.neighbours

    flip(tri0, 0, tri1, 2)

    assert tri0 in tri1.neighbours
    assert tri2 in tri1.neighbours
    assert tri5 in tri1.neighbours

    assert tri1 in tri0.neighbours
    assert tri3 in tri0.neighbours
    assert tri4 in tri0.neighbours

    assert "a" in tri0.vertices
    assert "b" in tri0.vertices
    assert "d" in tri0.vertices

    assert "a" in tri1.vertices
    assert "c" in tri1.vertices
    assert "d" in tri1.vertices

def test_split():
    conv = ToPointsAndSegments()
    conv.add_point((0, 0))
    conv.add_point((10, 0))
    conv.add_point((10, 20))
    close = (5, 4)
    conv.add_point(close)
    conv.add_point((0, 20))
    conv.add_segment((0,0), (10,0))
    conv.add_segment((10,0), (10,20))
    conv.add_segment((10,20), close)
    conv.add_segment(close, (0,20))
    conv.add_segment((0,20), (0,0))

    skel = calc_skel(conv, pause = True, output=True)
    assert len(skel.segments()) == 12


def test_left_right_for_vertex():
    kva = KineticVertex()
    kvb = KineticVertex()
    kvc = KineticVertex()
    kvb_ = KineticVertex()
    kvc_ = KineticVertex()

    kva.left = kvb, 0
    kva.right = kvc, 0
    print kva._left
    print kva._right

    kva.left = kvb_, 10
    kva.right = kvc_, 10
    print kva._left
    print kva._right

    assert kva.left_at(5) is kvb
    assert kva.left_at(15) is kvb_

    kva.left_at(-1)

def test_collinear_bisectors():
    from calc import bisector
    a = (10,0)
    b = (0,0)
    c = (10.,0)
    d = (-10, 0)
    bi = bisector(a, b, c)
    print "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format(b, bi)
    bi = bisector(c, b, a)
    print "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format(b, bi)
    bi = bisector(d, b, a)
    print bi
    bi = bisector(a, b, d)
    print bi
    p1, p2, p3 = (-9.171572875253812, 2.4142135623730936), (0.0, 2.414213562373094), (9.171572875253812, 2.4142135623730936)
    print bisector(p1, p2, p3)
    print bisector(p3, p2, p1)

def test_ordering():
    tri1 = ("1", )
    tri2 = ("2", )
    tri3 = ("3", )
    ea = Event(when = 1, tri = tri1)
    eb = Event(when = 1, tri = tri2)
    ec = Event(when = 2, tri = tri3)
    L = [ec, eb, ea]
    queue = OrderedSequence(cmp=compare_event_by_time)
    for e in L: queue.add(e)
    queue.remove(ec)
    queue.remove(Event(when=1, tri=tri1))
    while queue:
        e = queue.popleft()
        assert e.triangle is tri2
    assert len(queue) == 0

def helper_make_test_collapse_time():
    from math import sqrt
    from tri import triangulate, ToPointsAndSegments
    from grassfire import init_skeleton

    conv = ToPointsAndSegments()
    conv.add_polygon([[(-1,0), (0,-1), (1,0), (0,5), (-1,0)]])

    conv = ToPointsAndSegments()
    ring = [(0,0), (3,0), (3.8,2), (4,0), (6,0), (6.3, 2), (7,0), (10,0), (10,5), (7,5), (6.5, 3), (6,5), (4,5), (3.5,3), (3,5), (0,5), (0,0)]
    conv.add_polygon([ring])




    conv = ToPointsAndSegments()
    ring = [(0,0), (3,0), (3.8,2), (4,0), (6,0), (6.3, 2), (7,0), (10,0), (10,5), (7,5), (6.5, 3), (6,5), (4,5), (3.5,3), (3,5), (0,5), (0,0)]
    conv.add_polygon([ring])

    conv = ToPointsAndSegments()
    ring = [(0,0), (3,0), (1.5, sqrt(3)), (0,0)]
    conv.add_polygon([ring])
    dt = triangulate(conv.points, None, conv.segments)
    skel = init_skeleton(dt)
    print "triangles = {}"
    for t in skel.triangles:
        print "###", id(t)
        print "k = KineticTriangle()"
        print "V = []"
        for v in t.vertices:
            print "v = KineticVertex()"
            print "v.origin =", v.origin
            print "v.velocity =", v.velocity
            print "V.append(v)"
        print "k.vertices = V"
        print "triangles[",id(t),"] = k"
    print ""
    print "### neighbour relationships"
    for t in skel.triangles:
        print "n = [", ", ".join(["triangles[{0}]".format(id(n)) if n is not None else "None" for n in t.neighbours]), "]"
        print "triangles[",id(t),"].neighbours = n"


def test_flip_loop():
    """This has an infinite sequence in there...
    """
    conv = ToPointsAndSegments()
    ring = [ (0,0), (3,0), (3.8,2), (4,0), (6,0), (6.3, 2), (7,0), (10,0), (13,4), (10,5), (7,5), (6.5, 3), (6,5), (4,5), (3.5,3), (3,5), (0,5), (-2,2), (0,0)]
    conv.add_polygon([ring])
    skel = calc_skel(conv, pause=True, output=True)

def test_flip_loop2():
    """This makes an infinite event loop with infinite triangles (because of flipping)!
    """
    conv = ToPointsAndSegments()
    ring = [ (0,0), (3,0), (3.8,2), (5,0), (6.3, 2), (7,0), (10,0), (13,4), (10,5), (7,5), (6.5, 3), (5,5), (3.5,3), (3,5), (0,5), (-2,2), (0,0)]
    conv.add_polygon([ring])
    skel = calc_skel(conv, pause=True, output=True)

def test_inf_quad():
    conv = ToPointsAndSegments()
    ring = [(0,0), (4.2,5.4), (6.9,0.05), (10,2), (5,10), (0,0)]
    conv.add_polygon([ring])
    skel = calc_skel(conv, pause=True, output=True)

test_inf_quad()

def test_inf_flat():
    conv = ToPointsAndSegments()
    ring = [(0,0), (10,0), (50,1), (90,0), (100,0), (49,2), (0,0)]
    conv.add_polygon([ring])
    skel = calc_skel(conv, pause=True, output=True)


def test_inf_teeth():
    conv = ToPointsAndSegments()
    ring = [(0,7), (10,0), (8,8), (16,8), (17,0), (17.1,10.5), (0,10), (0,7)]
    conv.add_polygon([ring])
    skel = calc_skel(conv, pause=True, output=True)


def test_appendix_a4():
    """Contains zero triangle to flip ...
    """
    ring = [(0,0), (1,1), (2,0), (6,0), (6.1,3), (2.5, 5.5), (-0.5, 3), (0,0) ]
    conv = ToPointsAndSegments()
    conv.add_polygon([ring])
    skel = calc_skel(conv, output=True, pause=True)
    print "DONE"

def test_infinite():

    conv = ToPointsAndSegments()

    conv.add_point((0,0))
    conv.add_point((1,0))

    conv.add_segment((0,0), (1,0))

    calc_skel(conv, pause=True, output=True)

def output_edges_at_T(edges, T, fh):
    fh.write("id;side;wkt\n")
    for e in edges:
        segment = e.segment
        s = segment[0].position_at(T), segment[1].position_at(T)
        fh.write("{0};{1};LINESTRING({2[0][0]} {2[0][1]}, {2[1][0]} {2[1][1]})\n".format(id(e.triangle), e.side, s))

def output_triangles_at_T(tri, T, fh):
    """Output list of triangles as WKT to text file (for QGIS)"""
    fh.write("id;time;wkt;n0;n1;n2;v0;v1;v2;finite;info\n")
    for t in tri:
        if t is None:
            continue
        fh.write("{0};{6};{1};{2[0]};{2[1]};{2[2]};{3[0]};{3[1]};{3[2]};{4};{5}\n".format(id(t),
                                                                                      t.str_at(T),
                                                                                      [id(n) for n in t.neighbours],
                                                                                      [id(v) for v in t.vertices],
                                                                                      t.is_finite,
                                                                                      t.info,
                                                                                      T))


def test_infinite2():
    """3 segments with terminal vertices at convex hull
    """
    conv = ToPointsAndSegments()
    l0 = [ ( 0.032020441647887, 0.050549836508082), (0.556388841835153, 0.835771552524547) ]
    l1 = [ ( 0.597646254032629, 0.835771552524547), (1.133992612599807, 0.029255688277127) ]
    l2 = [ ( 1.118022001426591, -0.000023765540436), (0.065292548258754, -0.000023765540436) ]
    for line in l0, l1, l2:
        conv.add_point(line[0])
        conv.add_point(line[1])
        conv.add_segment(*line)
    skel = calc_skel(conv, pause = True, output = True)
    print skel.vertices
    print skel.triangles

    return




def test_infinite3():
    """6 segments with terminal vertices at convex hull
    """
    from math import sqrt
    conv = ToPointsAndSegments()
    l0 = [(0.0, 1.0), (1.0, 1.0)]
    l1 = [(1,1), (1,0)]
    l2 = [(5,0), (5,1,)]
    l3 = [(5,1), (6,1,)]
    l4 = [(2,2 + sqrt(3)/2.*4), (3,1 + sqrt(3)/2.*4,)]
    l5 = [(3,1 + sqrt(3)/2.*4), (4, 2 + sqrt(3)/2.*4)]
    for line in l0, l1, l2, l3, l4, l5:
        conv.add_point(line[0])
        conv.add_point(line[1])
        conv.add_segment(*line)
    skel = calc_skel(conv, pause = True, output = True)
    print skel.vertices
    print skel.triangles
    el = init_event_list(skel)
    event_loop(el, skel, pause=True)

    return

from PyQt4.QtCore import QFileSystemWatcher

isDone = True
ct = 0
def watch():
    global isDone, ct
    if isDone:
        isDone = False
        ct += 1
        iface.mapCanvas().refresh()
        iface.mapCanvas().saveAsImage("/tmp/image{0:04d}.png".format(ct))
        isDone = True

watcher = QFileSystemWatcher()
watcher.addPath( '/tmp/wavefront_edges_progress.wkt' )
watcher.fileChanged.connect( watch )






from PyQt4.QtCore import QFileSystemWatcher
watcher = QFileSystemWatcher()
watcher.addPath( '/tmp/rays.wkt' )
watcher.fileChanged.connect( iface.mapCanvas().refresh )


from PyQt4.QtCore import QFileSystemWatcher
watcher = QFileSystemWatcher()
iface.mapCanvas().setCachingEnabled(False)
isDone = True

def watch():
    global isDone
    print "watch() called"
    if isDone:
        print "refreshing"
        isDone = False
        iface.mapCanvas().refresh()
        isDone = True

watcher.addPath( '/tmp/signal' )
watcher.fileChanged.connect( watch )


def watch2():
    iface.mapCanvas().refresh()
