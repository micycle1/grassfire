import logging
from tri.delaunay.tds import cw, ccw
from grassfire.events.lib import stop_kvertices, compute_new_kvertex, update_circ, replace_kvertex
from grassfire.events.parallel import handle_parallel_fan
from grassfire.line2d import WaveFrontIntersector
def handle_split_event(evt, step, skel, queue, immediate, pause):
    """Handles a split event where a wavefront edge is hit on its interior
    This splits the wavefront in two pieces
    """
    t = evt.triangle


    assert len(evt.side) == 1
    e = evt.side[0]
    now = evt.time
    v = t.vertices[(e) % 3]
    n = t.neighbours[e]
    assert n is None
    v1 = t.vertices[(e + 1) % 3]
    v2 = t.vertices[(e + 2) % 3]



    assert v1.wfr is v2.wfl
    a = v.wfr
    b = v1.wfr
    assert v2.wfl is b
    c = v.wfl
    intersector = WaveFrontIntersector(a, b)           #
    bi0 = intersector.get_bisector()
    intersector = WaveFrontIntersector(a, c)           #
    bi1 = intersector.get_bisector()
    try:
        intersector = WaveFrontIntersector(a, b)
        pos_at_now = intersector.get_intersection_at_t(now)
    except ValueError:
        pass

    try:
        intersector = WaveFrontIntersector(a, c)
        pos_at_now = intersector.get_intersection_at_t(now)
    except ValueError:
        pass
    try:
        intersector = WaveFrontIntersector(b,c)
        pos_at_now = intersector.get_intersection_at_t(now)
    except ValueError:
        pass

    sk_node, newly_made = stop_kvertices([v], step, now)
    if newly_made:
        skel.sk_nodes.append(sk_node)

    assert v1.ur is v2.ul
    vb = compute_new_kvertex(v.ul, v2.ul, now, sk_node, len(skel.vertices) + 1, v.internal or v2.internal, pause)
    vb.wfl = v.wfl
    vb.wfr = v2.wfl
    skel.vertices.append(vb)

    if pause:
        from grassfire.inout import interactive_visualize
        interactive_visualize(queue, skel, step, now)
    va = compute_new_kvertex(v1.ur, v.ur, now, sk_node, len(skel.vertices) + 1, v.internal or v1.internal, pause)
    va.wfl = v1.wfr
    va.wfr = v.wfr

    skel.vertices.append(va)

    if pause:
        interactive_visualize(queue, skel, step, now)

    update_circ(v.left, vb, now)
    update_circ(vb, v2, now)
    assert vb.left.wfr is vb.wfl
    assert vb.right.wfl is vb.wfr

    update_circ(v1, va, now)
    update_circ(va, v.right, now)

    assert va.left.wfr is va.wfl
    assert va.right.wfl is va.wfr
    b = t.neighbours[(e + 1) % 3]
    assert b is not None
    b.neighbours[b.neighbours.index(t)] = None
    fan_b = replace_kvertex(b, v, vb, now, ccw, queue, immediate)

    if pause:
        from grassfire.inout import interactive_visualize
        interactive_visualize(queue, skel, step, now)
    a = t.neighbours[(e + 2) % 3]
    assert a is not None
    a.neighbours[a.neighbours.index(t)] = None
    fan_a = replace_kvertex(a, v, va, now, cw, queue, immediate)

    if pause:
        from grassfire.inout import interactive_visualize
        interactive_visualize(queue, skel, step, now)

    t.stops_at = now
    if va.inf_fast:
        handle_parallel_fan(fan_a, va, now, cw, step, skel, queue, immediate, pause)
    if vb.inf_fast:
        handle_parallel_fan(fan_b, vb, now, ccw, step, skel, queue, immediate, pause)