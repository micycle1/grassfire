"""Tests for line2d module"""

from grassfire.line2d import (
    Line2, LineLineIntersector, WaveFront, WaveFrontIntersector,
    rotate90cw
)
from grassfire.vectorops import mul, make_vector


def as_wkt(start, end):
    """Helper to format WKT linestring"""
    return "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format(start, end)


def main():

    print("")
    e0 = Line2.from_points( (-0.25, 0), (0, 1) )
    print(e0)
    e0t = e0.translated(mul(e0.w,50.))
    print(e0t)

    print("")
    e1 = Line2.from_points( (0, 1), (1,0) )
    print(e1)
    e1t = e1.translated(e1.w)
    print(e1t)

    print("")

    pts = []
    intersect = LineLineIntersector(e0, e1)

    if intersect.intersection_type() == 1:
        print(intersect.result)
        pts.append(intersect.result)

    intersect = LineLineIntersector(e0t, e1t)
    if intersect.intersection_type() == 1:
        print(intersect.result)
        pts.append(intersect.result)

    if len(pts) > 1:
        print(pts)
        velocity = make_vector(end=pts[1], start=pts[0])
        print(velocity)

    xaxis = Line2.from_points( (0, 0), (1, 0) )
    yaxis = Line2.from_points( (0, 0), (0, 1) )

    with open('/tmp/lines.wkt', 'w') as fh:
        fh.write("wkt")
        fh.write("\n")

        fh.write(as_wkt(*xaxis.visualize()))
        fh.write("\n")

        fh.write(as_wkt(*yaxis.visualize()))
        fh.write("\n")

        fh.write(as_wkt(*e0.visualize()))
        fh.write("\n")

        fh.write(as_wkt(*e0t.visualize()))
        fh.write("\n")

        fh.write(as_wkt(*(e0t.perpendicular( (10, -10) ).visualize())))
        fh.write("\n")


        fh.write(as_wkt(*e1.visualize()))
        fh.write("\n")
        fh.write(as_wkt(*e1t.visualize()))
        fh.write("\n")


    pos_xaxis = Line2.from_points( (0, 0), (-1, 0) )
    neg_xaxis_plus1 = Line2.from_points( (-1, -1), (0, 0) )
    bi_p = pos_xaxis.bisector(neg_xaxis_plus1).perpendicular((0, 0))
    bi_n = neg_xaxis_plus1.bisector(pos_xaxis)

    la = Line2.from_points( (-1, 0), (0, 0) )
    lb = Line2.from_points( (-1, -1), (0,0) )
    bi_ab  = la.bisector(lb)
    bi_ba  = lb.bisector(la)


    with open('/tmp/bisect.wkt', 'w') as fh:
        fh.write("wkt")
        fh.write("\n")

        fh.write(as_wkt(*la.visualize()))
        fh.write("\n")

        fh.write(as_wkt(*lb.visualize()))
        fh.write("\n")

        fh.write(as_wkt(*bi_ab.visualize()))
        fh.write("\n")

        fh.write(as_wkt(*bi_ba.visualize()))
        fh.write("\n")


def create_bevel_wavefront(first, last):
    """
    turning ccw around a vertex that needs to be made into a bevel
    vertex, we start at the first edge, and stop at the last edge
    """
    assert first.start == last.end
    perp = first.line.perpendicular(through=first.start)
    wf_bevel = WaveFront(first.start, first.start, line = perp)
    del perp
    print(wf_bevel)
    wf_intersect = WaveFrontIntersector(first, wf_bevel)
    bisector1 = wf_intersect.get_bisector()
    intersection = wf_intersect.get_intersection_at_t(0.0)
    print(bisector1)
    print(intersection)
    assert intersection == first.start == last.end
    wf_intersect = WaveFrontIntersector(wf_bevel, last)
    bisector2 = wf_intersect.get_bisector()
    intersection = wf_intersect.get_intersection_at_t(0.0)
    print(bisector2)
    print(intersection)
    assert intersection == first.start == last.end
    return wf_bevel, bisector1, bisector2


def test_wf():
    left = WaveFront(start=(-1,1), end= (0,0))
    right = WaveFront(start=(0,0), end=(1,1))
    print(left)
    print(right)
    wf_intersect = WaveFrontIntersector(left, right)
    bisector = wf_intersect.get_bisector()
    print(bisector)
    left = WaveFront(start=(0,0), end= (10,0))
    right = WaveFront(start=(10,0), end=(20,0))
    print(left)
    print(right)
    wf_intersect = WaveFrontIntersector(left, right)
    bisector = wf_intersect.get_bisector()
    print(bisector)
    left = WaveFront(start=(0,0), end=(0,10))
    right = WaveFront(start=(0,10), end=(0,0))
    print(left)
    print(right)
    wf_intersect = WaveFrontIntersector(left, right)
    bisector = wf_intersect.get_bisector()
    print(bisector)


def test_create_bevel_wavefront():
    first = WaveFront(start=(10,12.5), end=(0,0))
    last = WaveFront(start=(0,0), end=(10,12.5))
    bevel0, bi1, bi2 = create_bevel_wavefront(first, last)

    bevel1, bi1, bi2 = create_bevel_wavefront(last, first)

    with open('/tmp/lines.wkt', 'w') as fh:
        fh.write("WKT")
        fh.write("\n")

    with open('/tmp/lines.wkt', 'a') as fh:
        for t in range(1, 10):
            for item in [
                    first.line.at_time(t).visualize(),
                    last.line.at_time(t).visualize(),
                    bevel0.line.at_time(t).visualize(),
                    bevel1.line.at_time(t).visualize()]:
                fh.write(item)
                fh.write("\n")


if __name__ == "__main__":
    test_wf()
    test_create_bevel_wavefront()
