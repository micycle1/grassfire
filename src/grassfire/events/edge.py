import logging

from tri.delaunay.tds import cw, ccw, Edge

from grassfire.events.lib import stop_kvertices, compute_new_kvertex, \
    update_circ, replace_kvertex, schedule_immediately, near_zero
from grassfire.events.lib import get_fan, is_infinitely_fast
from grassfire.events.parallel import handle_parallel_fan
from grassfire.line2d import WaveFrontIntersector
def handle_edge_event(evt, step, skel, queue, immediate, pause):
    """Handles triangle collapse, where exactly 1 edge collapses"""
    t = evt.triangle

    assert len(evt.side) == 1, len(evt.side)
    e = evt.side[0]
    is_wavefront_collapse = t.neighbours[e] is None
    now = evt.time
    v1 = t.vertices[ccw(e)]
    v2 = t.vertices[cw(e)]
    if is_wavefront_collapse and not v1.is_stopped and not v2.is_stopped:
        assert v1.right is v2
        assert v2.left is v1
    a = v1.wfl
    b = v1.wfr
    if is_wavefront_collapse and not v1.is_stopped and not v2.is_stopped:
        assert v2.wfl is b
    c = v2.wfr
    intersector = WaveFrontIntersector(a, c)
    bi = intersector.get_bisector()
    pos_at_now = None
    try:
        intersector = WaveFrontIntersector(a, c)
        pos_at_now = intersector.get_intersection_at_t(now)
        logging.debug("POINT({0[0]} {0[1]});a;c".format(pos_at_now))
    except ValueError:
        pass

    sk_node, newly_made = stop_kvertices([v1, v2], step, now, pos=pos_at_now)
    if newly_made:
        skel.sk_nodes.append(sk_node)
    kv = compute_new_kvertex(v1.ul, v2.ur, now, sk_node, len(skel.vertices) + 1, v1.internal or v2.internal, pause)
    kv.wfl = v1.wfl
    kv.wfr = v2.wfr

    if v1.left:
        logging.debug(v1.left.position_at(now))
    else:
        logging.warning("no v1.left")
    if v2.right:
        logging.debug(v2.right.position_at(now))
    else:
        logging.warning("no v2.right")
    skel.vertices.append(kv)
    update_circ(v1.left, kv, now)
    update_circ(kv, v2.right, now)
    assert kv.wfl is kv.left.wfr
    assert kv.wfr is kv.right.wfl
    a = t.neighbours[ccw(e)]
    b = t.neighbours[cw(e)]
    n = t.neighbours[e]

    fan_a = []
    fan_b = []
    if a is not None:
        a_idx = a.neighbours.index(t)
        a.neighbours[a_idx] = b
        fan_a = replace_kvertex(a, v2, kv, now, cw, queue, immediate)

        if fan_a:
            e = Edge(fan_a[-1], cw(fan_a[-1].vertices.index(kv)))
            orig, dest = e.segment
            import math
            if(near_zero(math.sqrt(orig.distance2_at(dest, now)))):
                schedule_immediately(fan_a[-1], now, queue, immediate)
    if b is not None:
        b_idx = b.neighbours.index(t)
        b.neighbours[b_idx] = a
        fan_b = replace_kvertex(b, v1, kv, now, ccw, queue, immediate)
        if fan_b:
            e = Edge(fan_b[-1], ccw(fan_b[-1].vertices.index(kv)))
            orig, dest = e.segment
            import math
            if(near_zero(math.sqrt(orig.distance2_at(dest, now)))):
                schedule_immediately(fan_b[-1], now, queue, immediate)

    if n is not None:
        n.neighbours[n.neighbours.index(t)] = None
        if n.event is not None and n.stops_at is None:
            schedule_immediately(n, now, queue, immediate)
    t.stops_at = now
    if kv.inf_fast:
        if fan_a and fan_b:
            fan_a = list(fan_a)
            fan_a.reverse()
            fan_a.extend(fan_b)
            handle_parallel_fan(fan_a, kv, now, ccw, step, skel, queue, immediate, pause)
            return
        elif fan_a:
            handle_parallel_fan(fan_a, kv, now, cw, step, skel, queue, immediate, pause)
            return
        elif fan_b:
            handle_parallel_fan(fan_b, kv, now, ccw, step, skel, queue, immediate, pause)
            return


def handle_edge_event_3sides(evt, step, skel, queue, immediate):
    """Handle a collapse of a triangle with 3 sides collapsing.
    It does not matter whether the 3-triangle has wavefront edges or not.

    Important: The triangle vertices should collapse to 1 point.

    The following steps are performed:
    - stop the 3 kinetic vertices of the triangle
    - optionally make a new skeleton node
    - schedule all neighbours, if any, for immediate processing
      (these also collapse to the same point)
    """
    now = evt.time
    t = evt.triangle


    assert len(evt.side) == 3
    sk_node, newly_made = stop_kvertices(t.vertices, step, now)
    if newly_made:
        skel.sk_nodes.append(sk_node)
    for n in t.neighbours:
        if n is not None and n.event is not None and n.stops_at is None:
            n.neighbours[n.neighbours.index(t)] = None
            schedule_immediately(n, now, queue, immediate)
    t.stops_at = now


def handle_edge_event_1side(evt, step, skel, queue, immediate, pause):
    """Handle a collapse of a triangle with 1 side collapsing.

    Important: The triangle collapses to a line segment.
    """
    t = evt.triangle



    assert len(evt.side) == 1, len(evt.side)
    e = evt.side[0]
    now = evt.time
    v0 = t.vertices[e]
    v1 = t.vertices[ccw(e)]
    v2 = t.vertices[cw(e)]
    sk_node, newly_made = stop_kvertices([v1, v2], step, now)
    if newly_made:
        skel.sk_nodes.append(sk_node)
    kv = compute_new_kvertex(v1.ul, v2.ur, now, sk_node, len(skel.vertices) + 1, v1.internal or v2.internal, pause)
    skel.vertices.append(kv)
    sk_node, newly_made = stop_kvertices([v0, kv], step, now)
    if newly_made:
        skel.sk_nodes.append(sk_node)
    t.stops_at = now

