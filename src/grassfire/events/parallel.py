import logging

from tri.delaunay.tds import ccw, cw, Edge
from tri.delaunay.tds import apex, orig, dest

from grassfire.events.lib import stop_kvertices, update_circ, \
    compute_new_kvertex, replace_kvertex, schedule_immediately, \
    is_infinitely_fast
from grassfire.events.lib import get_fan
from grassfire.inout import interactive_visualize
from grassfire.calc import near_zero
from grassfire.primitives import KineticVertex
from grassfire.vectorops import dist



def flip(t0, side0, t1, side1):
    """Performs a flip of triangle t0 and t1

    If t0 and t1 are two triangles sharing a common edge AB,
    the method replaces ABC and BAD triangles by DCA and DBC, respectively.

    Pre-condition: input triangles share a common edge and this edge is known.
    """
    apex0, orig0, dest0 = apex(side0), orig(side0), dest(side0)
    apex1, orig1, dest1 = apex(side1), orig(side1), dest(side1)
    assert t0.vertices[orig0] is t1.vertices[dest1]
    assert t0.vertices[dest0] is t1.vertices[orig1]
    assert t0.neighbours[apex0] is not None
    assert t1.neighbours[apex1] is not None
    A, B, C, D = t0.vertices[apex0], t0.vertices[
        orig0], t1.vertices[apex1], t0.vertices[dest0]
    AB, BC, CD, DA = t0.neighbours[dest0], t1.neighbours[
        orig1], t1.neighbours[dest1], t0.neighbours[orig0]
    apex_around = []
    for neighbour, corner in zip([AB, BC, CD, DA],
                                 [A, B, C, D]):
        if neighbour is None:
            apex_around.append(None)
        else:
            apex_around.append(ccw(neighbour.vertices.index(corner)))
    for neighbour, side, t in zip([AB, BC, CD, DA],
                                  apex_around,
                                  [t0, t0, t1, t1]):
        if neighbour is not None:
            neighbour.neighbours[side] = t
    t0.vertices = [A, B, C]
    t0.neighbours = [BC, t1, AB]
    t1.vertices = [C, D, A]
    t1.neighbours = [DA, t0, CD]
def handle_parallel_fan(fan, pivot, now, direction, step, skel, queue, immediate, pause):
    """Dispatches to correct function for handling parallel wavefronts

    fan: list of triangles, sorted (in *direction* order)
    pivot: the vertex that is going infinitely fast
    now: time at which parallel event happens
    direction: how the fan turns

    skel: skeleton structure
    queue: event queue
    immediate: queue of events that should be dealt with when finished handling
    the current event

    pause: whether we should stop for interactivity
    """
    if not fan:
        raise ValueError("we should receive a fan of triangles to handle them")
        return

    if pause:
        logging.debug(" -- {}".format(len(fan)))
        logging.debug("    triangles in the fan: {}".format([(id(_), _.info) for _ in fan]))
        logging.debug("    fan turns: {}".format(direction))
        logging.debug("    pivot: {} [{}]".format(id(pivot), pivot.info))
        interactive_visualize(queue, skel, step, now)

    assert pivot.inf_fast
    first_tri = fan[0]
    last_tri = fan[-1]
    if first_tri.neighbours.count(None) == 3:
        assert first_tri is last_tri #FIXME: is this true?
        dists = []
        for side in range(3):
            edge = Edge(first_tri, side)
            d = dist(*map(lambda x: x.position_at(now), edge.segment))
            dists.append(d)
        dists_sub_min = [near_zero(_ - min(dists)) for _ in dists]
        if near_zero(min(dists)) and dists_sub_min.count(True) == 1:
            logging.debug("Smallest edge collapses? {}".format(near_zero(min(dists))))
            assert dists_sub_min.count(True) == 1
            side = dists_sub_min.index(True)
            pivot = first_tri.vertices[dists_sub_min.index(True)]
            handle_parallel_edge_event_even_legs(first_tri, first_tri.vertices.index(pivot), pivot, now, step, skel, queue, immediate)
            return
        else:
            handle_parallel_edge_event_3tri(first_tri, first_tri.vertices.index(pivot), pivot, now, step, skel, queue, immediate)
            return

    if first_tri is last_tri:
        assert len(fan) == 1
    if direction is cw:
        left = fan[0]
        right = fan[-1]
    else:
        assert direction is ccw
        left = fan[-1]
        right = fan[0]

    left_leg_idx = ccw(left.vertices.index(pivot))
    left_leg = Edge(left, left_leg_idx)
    left_dist = dist(*map(lambda x: x.position_at(now), left_leg.segment))
    right_leg_idx = cw(right.vertices.index(pivot))
    right_leg = Edge(right, right_leg_idx)
    right_dist = dist(*map(lambda x: x.position_at(now), right_leg.segment))
    dists = [left_dist, right_dist]
    dists_sub_min = [near_zero(_ - min(dists)) for _ in dists]
    unique_dists = dists_sub_min.count(True)
    if unique_dists == 2:
        if len(fan) == 1:
            handle_parallel_edge_event_even_legs(first_tri, first_tri.vertices.index(pivot), pivot, now, step, skel, queue, immediate)
        elif len(fan) == 2:

            all_2 = True
            for t in fan:
                left_leg_idx = ccw(t.vertices.index(pivot))
                left_leg = Edge(t, left_leg_idx)
                left_dist = dist(*map(lambda x: x.position_at(now), left_leg.segment))

                right_leg_idx = cw(t.vertices.index(pivot))
                right_leg = Edge(t, right_leg_idx)
                right_dist = dist(*map(lambda x: x.position_at(now), right_leg.segment))

                dists = [left_dist, right_dist]

                dists_sub_min = [near_zero(_ - min(dists)) for _ in dists]
                logging.debug("  {}".format([left_dist, right_dist]))
                logging.debug("  {}".format(dists_sub_min))
                unique_dists = dists_sub_min.count(True)
                if unique_dists != 2:
                    all_2 = False
            if all_2 == True:
                for t in fan:
                    handle_parallel_edge_event_even_legs(t, t.vertices.index(pivot), pivot, now, step, skel, queue, immediate)
            else:
                t0 = fan[0]
                t1 = fan[1]
                side0 = t0.neighbours.index(t1)
                side1 = t1.neighbours.index(t0)
                flip(t0, side0, t1, side1)
                t0_has_inf_fast = [v.inf_fast for v in t0.vertices]
                t1_has_inf_fast = [v.inf_fast for v in t1.vertices]

                if True in t0_has_inf_fast:
                    handle_parallel_edge_event_even_legs(t0, t0.vertices.index(pivot), pivot, now, step, skel, queue, immediate)

                if True in t1_has_inf_fast:
                    handle_parallel_edge_event_even_legs(t1, t1.vertices.index(pivot), pivot, now, step, skel, queue, immediate)

            if pause:
                interactive_visualize(queue, skel, step, now)

        else:
            raise NotImplementedError('More than 2 triangles in parallel fan with 2 equal sized legs on the outside -- does this exist?')
        if len(fan) == 1 or (len(fan) == 2 and all_2 == True):
            for t in fan:
                assert t.stops_at is not None
    else:
        shortest_idx = dists_sub_min.index(True)
        if shortest_idx == 1: # right is shortest, left is longest
            handle_parallel_edge_event_shorter_leg(right_leg.triangle, right_leg.side, pivot, now, step, skel, queue, immediate, pause)
        elif shortest_idx == 0: # left is shortest, right is longest
            handle_parallel_edge_event_shorter_leg(left_leg.triangle, left_leg.side, pivot, now, step, skel, queue, immediate, pause)


def handle_parallel_edge_event_shorter_leg(t, e, pivot, now, step, skel, queue, immediate, pause):
    """Handles triangle collapse, where exactly 1 edge collapses

    One of the vertices of the triangle moves *infinitely* fast.

    There are 2 cases handled in this function

    a. triangle with long left leg, short right leg
    b. triangle with long right leg, short left leg

    Arguments:
    t -- triangle that collapses
    e -- short side over which pivot moves inf fast

    """


    assert pivot.inf_fast

    v1 = t.vertices[ccw(e)]
    v2 = t.vertices[cw(e)]
    v3 = t.vertices[e]
    assert pivot is v1 or pivot is v2

    to_stop = []
    for v in [v1, v2]:
        if not v.inf_fast:
            to_stop.append(v)
    sk_node, newly_made = stop_kvertices(to_stop, step, now)
    if newly_made:
        skel.sk_nodes.append(sk_node)
    if pivot.stop_node is None:
        assert pivot.stop_node is None
        assert pivot.stops_at is None
        pivot.stop_node = sk_node
        pivot.stops_at = now
    t.stops_at = now
    assert t.vertices.index(pivot) != e
    kv = compute_new_kvertex(v1.ul, v2.ur, now, sk_node, len(skel.vertices) + 1, v1.internal or v2.internal, pause)
    kv.wfl = v1.left.wfr
    kv.wfr = v2.right.wfl

    a = t.neighbours[ccw(e)]
    b = t.neighbours[cw(e)]
    n = t.neighbours[e]
    skel.vertices.append(kv)
    update_circ(v1.left, kv, now)
    update_circ(kv, v2.right, now)
    fan_a = []
    fan_b = []
    if a is not None:
        logging.debug("- replacing vertex for neighbours at side A {} [{}]".format(id(a), a.info))
        a_idx = a.neighbours.index(t)
        a.neighbours[a_idx] = b
        fan_a = replace_kvertex(a, v2, kv, now, cw, queue, immediate)
        if pause:
            interactive_visualize(queue, skel, step, now)

    if b is not None:
        logging.debug("- replacing vertex for neighbours at side B {} [{}]".format(id(b), b.info))
        b_idx = b.neighbours.index(t)
        b.neighbours[b_idx] = a
        fan_b = replace_kvertex(b, v1, kv, now, ccw, queue, immediate)
        if pause:
            interactive_visualize(queue, skel, step, now)
    if n is not None:
        n.neighbours[n.neighbours.index(t)] = None
        if n.event is not None and n.stops_at is None:
            schedule_immediately(n, now, queue, immediate)
    if kv and kv.inf_fast:

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
def handle_parallel_edge_event_even_legs(t, e, pivot, now, step, skel, queue, immediate):


    assert t.vertices.index(pivot) == e
    assert t.vertices[e] is pivot
    v1 = t.vertices[ccw(e)]
    v2 = t.vertices[cw(e)]
    sk_node, newly_made = stop_kvertices([v1,v2], step, now)
    if newly_made:
        skel.sk_nodes.append(sk_node)
    pivot.stop_node = sk_node
    pivot.stops_at = now
    t.stops_at = now

    n = t.neighbours[e]
    msg = "schedule adjacent neighbour for *IMMEDIATE* processing" if n is not None else "no neighbour to collapse simultaneously"
    if n is not None:
        n.neighbours[n.neighbours.index(t)] = None
        if n.event is not None and n.stops_at is None:
            schedule_immediately(n, now, queue, immediate)


def handle_parallel_edge_event_3tri(t, e, pivot, now, step, skel, queue, immediate):

    assert t.vertices.index(pivot) == e
    assert t.vertices[e] is pivot

    left_leg_idx = ccw(t.vertices.index(pivot))
    left_leg = Edge(t, left_leg_idx)
    left_dist = dist(*map(lambda x: x.position_at(now), left_leg.segment))
    v1 = t.vertices[ccw(e)]

    right_leg_idx = cw(t.vertices.index(pivot))
    right_leg = Edge(t, right_leg_idx)
    right_dist = dist(*map(lambda x: x.position_at(now), right_leg.segment))
    v2 = t.vertices[cw(e)]

    assert v1 is not pivot
    assert v2 is not pivot

    assert pivot in t.vertices
    assert v1 in t.vertices
    assert v2 in t.vertices

    from grassfire.vectorops import dot, norm
    magn_v1 = norm(v1.velocity)
    magn_v2 = norm(v2.velocity)


    dists = [left_dist, right_dist]
    dists_sub_min = [near_zero(_ - min(dists)) for _ in dists]
    if magn_v2 < magn_v1:
        sk_node, newly_made = stop_kvertices([v2], step, now)
        if newly_made:
            skel.sk_nodes.append(sk_node)
        v1.stop_node = sk_node
        v1.stops_at = now
    else:
        sk_node, newly_made = stop_kvertices([v1], step, now)
        if newly_made:
            skel.sk_nodes.append(sk_node)
        v2.stop_node = sk_node
        v2.stops_at = now
    pivot.stop_node = sk_node
    pivot.stops_at = now
    t.stops_at = now

    for kv in t.vertices:
        assert kv.stops_at is not None
