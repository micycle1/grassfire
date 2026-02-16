import logging
from tri.delaunay.tds import cw, ccw
from grassfire.events.lib import stop_kvertices, compute_new_kvertex, update_circ, replace_kvertex
from grassfire.events.parallel import handle_parallel_fan
from grassfire.line2d import WaveFrontIntersector


# ------------------------------------------------------------------------------
# Split event handler
def handle_split_event(evt, step, skel, queue, immediate, pause):
    """Handles a split event where a wavefront edge is hit on its interior
    This splits the wavefront in two pieces
    """
    t = evt.triangle

    logging.info("* split          :: tri>> #{} [{}]".format(id(t), t.info))

    logging.debug("{}".format(t.neighbours))
    assert len(evt.side) == 1
    e = evt.side[0]
    now = evt.time
    v = t.vertices[(e) % 3]
    n = t.neighbours[e]
    assert n is None
    v1 = t.vertices[(e + 1) % 3]
    v2 = t.vertices[(e + 2) % 3]

    logging.debug("v1 := {} [{}]".format(id(v1), v1.info))
    logging.debug("v2 := {} [{}]".format(id(v2), v2.info))


    assert v1.wfr is v2.wfl

    # ---- new use of wavefronts ------------------------------ #
    # split leads to 2 bisectors
    a = v.wfr
    b = v1.wfr
    assert v2.wfl is b
    c = v.wfl

    # TODO: are we having the correct bisectors here?
    intersector = WaveFrontIntersector(a, b)           #
    bi0 = intersector.get_bisector()
    # print(bi0)
    #
    intersector = WaveFrontIntersector(a, c)           #
    bi1 = intersector.get_bisector()
    # print(bi1)

    # the position of the node is witnessed by 3 pairs of wavefronts
    # NOTE uses try-except to determine parallel (bad approach)
    try:
        intersector = WaveFrontIntersector(a, b)
        pos_at_now = intersector.get_intersection_at_t(now)
    except ValueError:
        pass
    # print("POINT({0[0]} {0[1]})".format(pos_at_now))

    try:
        intersector = WaveFrontIntersector(a, c)
        pos_at_now = intersector.get_intersection_at_t(now)
        # print("POINT({0[0]} {0[1]})".format(pos_at_now))
    except ValueError:
        pass
    try:
        intersector = WaveFrontIntersector(b,c)
        pos_at_now = intersector.get_intersection_at_t(now)
        # print("POINT({0[0]} {0[1]})".format(pos_at_now))            #
    except ValueError:
        pass

    # ---- new use of wavefronts ------------------------------ #

    sk_node, newly_made = stop_kvertices([v], step, now)
    # add the skeleton node to the skeleton
    if newly_made:
        skel.sk_nodes.append(sk_node)
#     assert v1.right is v2
#     assert v2.left is v1

    assert v1.ur is v2.ul

    # the position of the stop_node can be computed by:
    # taking the original wavefronts and translate these lines to 'now'
    # then intersect them -> stop_node position
    # maybe this is more robust than relying on the direction vector and the velocity of the original point

    # a bisector based on the original line equations
    # BI = compute_crossing_bisector(v.ul, v2.ul, now)
    vb = compute_new_kvertex(v.ul, v2.ul, now, sk_node, len(skel.vertices) + 1, v.internal or v2.internal, pause)
    # FIXME: new wavefront
    vb.wfl = v.wfl
    vb.wfr = v2.wfl
#     logging.debug("""
# -- BISECTOR BI   {}
#             bi_b {}
#             vb.v {}""".format(BI, bi_b, vb.velocity)
#     )
    skel.vertices.append(vb)

    # if not vb.inf_fast:
    #     if sign(vb.velocity[0]) != sign(bi_b[0]) or sign(vb.velocity[1]) != sign(bi_b[1]):
    #         vb.velocity = mul(unit(bi_b), norm(vb.velocity))

    if pause:
        logging.debug('split l.194 -- computed new vertex B')
        from grassfire.inout import interactive_visualize
        interactive_visualize(queue, skel, step, now)

    # BI = compute_crossing_bisector(v1.ur, v.ur, now)
    va = compute_new_kvertex(v1.ur, v.ur, now, sk_node, len(skel.vertices) + 1, v.internal or v1.internal, pause)
    va.wfl = v1.wfr
    va.wfr = v.wfr

    skel.vertices.append(va)

    if pause:
        logging.debug('split l.211  -- computed new vertex A')
        interactive_visualize(queue, skel, step, now)

    # if not va.inf_fast:
    #     if sign(va.velocity[0]) != sign(bi_a[0]) or sign(va.velocity[1]) != sign(bi_a[1]):
    #         va.velocity = mul(unit(bi_a), norm(va.velocity))

    logging.debug("-- update circular list at B-side: {} [{}]".format(id(vb), vb.info))
    update_circ(v.left, vb, now)
    update_circ(vb, v2, now)

    # FIXME: why do these assertions not hold?
    assert vb.left.wfr is vb.wfl
    assert vb.right.wfl is vb.wfr

    logging.debug("-- update circular list at A-side: {} [{}]".format(id(va), va.info))
    update_circ(v1, va, now)
    update_circ(va, v.right, now)

    logging.debug("-- [{}]".format(va.info))
    logging.debug("   [{}]".format(va.right.info))
    logging.debug("   {}".format(va.right.wfl))
    logging.debug("   {}".format(va.wfr))

    # FIXME: why do these assertions not hold?
    assert va.left.wfr is va.wfl
    assert va.right.wfl is va.wfr

    # updates (triangle fan) at neighbour 1
    b = t.neighbours[(e + 1) % 3]
    assert b is not None
    b.neighbours[b.neighbours.index(t)] = None
    fan_b = replace_kvertex(b, v, vb, now, ccw, queue, immediate)

    if pause:
        logging.debug('split l.243 -- replaced vertex B')
        from grassfire.inout import interactive_visualize
        interactive_visualize(queue, skel, step, now)


    # updates (triangle fan) at neighbour 2
    a = t.neighbours[(e + 2) % 3]
    assert a is not None
    a.neighbours[a.neighbours.index(t)] = None
    fan_a = replace_kvertex(a, v, va, now, cw, queue, immediate)

    if pause:
        logging.debug('split l.255 -- replaced vertex A')
        from grassfire.inout import interactive_visualize
        interactive_visualize(queue, skel, step, now)

    t.stops_at = now

    # handle infinitely fast vertices
    if va.inf_fast:
        handle_parallel_fan(fan_a, va, now, cw, step, skel, queue, immediate, pause)
    if vb.inf_fast:
        handle_parallel_fan(fan_b, vb, now, ccw, step, skel, queue, immediate, pause)