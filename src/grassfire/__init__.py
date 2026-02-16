from tri.delaunay.helpers import ToPointsAndSegments
from tri.delaunay.insert_kd import triangulate
from tri.delaunay.iter import FiniteEdgeIterator, TriangleIterator
from tri.delaunay.inout import output_triangles

from grassfire.inout import output_offsets, output_skel
from grassfire.initialize import init_skeleton, internal_only_skeleton
from grassfire.events import init_event_list, event_loop
from grassfire.transform import get_transform, get_box

__version__ = "0.1.dev0"
__license__ = "MIT License"
__author__ = "Martijn Meijers"
__all__ = ["calc_skel"]


def calc_skel(conv, pause=False, output=False, shrink=True, internal_only=False):
    """Perform the calculation of the skeleton, given points and segments

    Returns:
        skel -- skeleton structure
    """
    if shrink:
        box = get_box(conv.points)
        transform = get_transform(box)
        pts = list(map(transform.forward, conv.points))
    else:
        pts = conv.points
    dt = triangulate(pts, conv.infos, conv.segments, output)
    if output:
        with open("/tmpfast/edges.wkt", "w") as fh:
            fh.write("id;wkt\n")
            edgeit = FiniteEdgeIterator(dt, constraints_only=True)
            for j, edge in enumerate(edgeit):
                fh.write(
                    "{0};LINESTRING({1[0][0]} {1[0][1]}, {1[1][0]} {1[1][1]})\n".format(
                        j, edge.segment
                    )
                )
    skel = init_skeleton(dt)
    if internal_only:
        skel = internal_only_skeleton(skel)
    if shrink:
        skel.transform = transform

    for kv in skel.vertices:
        x, y = kv.start_node.pos
        assert -2.0 <= x <= 2.0, (x, "start")
        assert -2.0 <= y <= 2.0, (y, "start")
    el = init_event_list(skel)
    last_evt_time = event_loop(el, skel, pause)
    if output:
        output_offsets(skel, last_evt_time)
        output_skel(skel, last_evt_time + 10)
        from grassfire.inout import visualize

        visualize([], skel, last_evt_time + 10)
    return skel


def calc_offsets(skel, now, ct=100):
    inc = now / ct
    for t in range(ct):
        t *= inc
        for v in skel.vertices:
            if (v.starts_at <= t and v.stops_at is not None and v.stops_at > t) or (
                v.starts_at <= t and v.stops_at is None
            ):
                try:
                    yield (
                        v.position_at(t),
                        v.left_at(t).position_at(t),
                        t,
                        id(v),
                        id(v.left_at(t)),
                    )
                except AttributeError:
                    continue
