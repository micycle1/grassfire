import logging

from grassfire.calc import is_close, near_zero
from grassfire.collapse import (compute_collapse_time,
                                compute_new_edge_collapse_event)
from grassfire.primitives import KineticVertex, SkeletonNode
from grassfire.vectorops import mul


def is_infinitely_fast(fan, now):
    """Determine whether all triangles in the fan collapse
    at the same time, if so, the vertex needs to be infinitely fast"""
    times = [tri.event.time if tri.event is not None else -1 for tri in fan]
    is_inf_fast = all(map(near_zero, [time - now for time in times]))
    if fan and is_inf_fast:
        return True
    else:
        return False


def stop_kvertices(V, step, now, pos=None):
    """ Stop a list of kinetic vertices *V* at time *now*, creating a new node.

    If one of the vertices was already stopped before, at a node, use that
    skeleton node

    Returns tuple of (new node, False) in case all vertices are stopped for the
    first time, otherwise it returns (node, True) to indicate that were already
    stopped once.
    """
    sk_node = None

    logging.debug("stopping kinetic vertices, @t:={}".format(now))
    for v in V:
        logging.debug(" - kv #{} [{}] inf-fast:={}".format(id(v), v.info, v.inf_fast))

    for v in V:
        stopped = v.stops_at is not None
        time_close = near_zero(v.starts_at - now)
        logging.debug("Vertex starts at same time as now: {}".format(time_close))
        logging.debug("Kinetic vertex is not stopped: {}".format(stopped))
        if stopped:
            logging.debug("Stop_node of vertex")
            sk_node = v.stop_node
        elif time_close:
            logging.debug("Time close")
            assert not stopped
            sk_node = v.start_node
        else:
            v.stops_at = now
    if sk_node is not None:
        logging.debug("Skeleton node already there")
        for v in V:
            v.stop_node = sk_node
            v.stops_at = now
        is_new_node = False
    else:
        if pos is None:
            logging.debug("Make new skeleton node")
            l = [v.position_at(now) for v in V]
            ct = len(l)
            sumx, sumy = 0., 0.
            for x, y in l:
                sumx += x
                sumy += y
            pos = sumx / ct, sumy / ct
        else:
            logging.debug("Make new skeleton node - using external position: {}".format(pos))
        sk_node = SkeletonNode(pos, step)
        for v in V:
            v.stop_node = sk_node
        is_new_node = True
    for v in V:
        assert v.stop_node is not None
        assert v.is_stopped == True
        assert v.stops_at == now

    logging.debug("POINT({0[0]} {0[1]});sknode_new_pos".format(sk_node.pos))
    return sk_node, is_new_node


def compute_new_kvertex(ul, ur, now, sk_node, info, internal, pause=False):
    """Based on the two wavefront directions and time t=now, compute the
    velocity and position at t=0 and return a new kinetic vertex

    Returns: KineticVertex
    """
    kv = KineticVertex()
    kv.info = info
    kv.starts_at = now
    kv.start_node = sk_node
    kv.internal = internal

    logging.debug('== New vertex: {} [{}] =='.format(id(kv), kv.info))
    logging.debug('bisector calc')
    logging.debug(' ul: {}'.format(ul))
    logging.debug(' ur: {}'.format(ur))
    logging.debug(' sk_node.pos: {}'.format(sk_node.pos))
    import math

    from grassfire.vectorops import add, angle_unit
    logging.debug(' >>> {}'.format(angle_unit(ul.w, ur.w)))
    u1, u2 = ul.w, ur.w
    direction = add(u1, u2)
    logging.debug(" direction: {}".format(direction))
    d, acos_d = angle_unit(u1, u2)
    if pause:
        with open('/tmpfast/support_lines.wkt', 'w') as fh:
            from grassfire.inout import interactive_visualize
            fh.write("wkt\tlr\toriginal")
            fh.write("\n")
            b = ul.translated(mul(ul.w,now)).bisector( ur.translated(mul(ur.w,now)) )
            bperp = b.perpendicular(sk_node.pos)

            for l, lr, original in zip([ul, ur, ul.translated(mul(ul.w,now)), ur.translated(mul(ur.w,now)), b, bperp], ["l", "r", "l", "r", "b", "b"], [True, True, False, False, None, None]):
                fh.write(l.at_time(0).visualize())
                fh.write("\t")
                fh.write(lr)
                fh.write("\t")
                fh.write(str(original))
                fh.write("\n")
    if all(map(near_zero, direction)) or near_zero(acos_d - math.pi) or d < math.cos(math.radians(179.999999)):
        logging.debug(" OVERRULED - vectors cancel each other out / angle ~180° -> parallel wavefront!")
        bi = (0, 0)
    else:
        from grassfire.line2d import (LineLineIntersectionResult,
                                      LineLineIntersector, make_vector)
        intersect = LineLineIntersector(ul, ur)
        tp = intersect.intersection_type()
        if tp == LineLineIntersectionResult.NO_INTERSECTION:
            bi = (0, 0)
        elif tp == LineLineIntersectionResult.POINT:
            pos_at_t0 = intersect.result
            ul_t = ul.translated(ul.w)
            ur_t = ur.translated(ur.w)
            intersect_t = LineLineIntersector(ul_t, ur_t)
            intersect_result_t = intersect_t.intersection_type()
            assert intersect_result_t == LineLineIntersectionResult.POINT
            bi = make_vector(end=intersect_t.result, start=pos_at_t0)
        elif tp == LineLineIntersectionResult.LINE:
            bi = tuple(ul.w[:])
            neg_velo = mul(mul(bi, -1.0), now)
            pos_at_t0 = add(sk_node.pos, neg_velo)
        else: # Handle unexpected case to avoid unbound variable
             raise RuntimeError(f"Unknown intersection type: {tp}")

    kv.velocity = bi #was: bisector(ul, ur)
    logging.debug(' kv.velocity: {}'.format(kv.velocity))
    from grassfire.vectorops import norm
    magn_v = norm(kv.velocity)
    logging.debug(' magnitude of velocity: {}'.format(magn_v))
    logging.debug('== New vertex ==')
    if kv.velocity == (0, 0): ### or abs(kv.velocity[0]) > 100 or abs(kv.velocity[1]) > 100:
        kv.inf_fast = True
        kv.origin = sk_node.pos
    else:
        kv.origin = pos_at_t0
    kv.ul = ul
    kv.ur = ur
    return kv


def get_fan(t, v, direction):
    """Gets a list of triangles that are the fan of
    vertex *v*, while turning *direction*, starting at triangle *t*

    This function assumes that the fan is finite (i.e. passes
    a triangle that has as neighbour = None (wavefront))
    """
    fan = []
    start = t
    while t is not None:
        side = t.vertices.index(v)
        fan.append(t)
        t = t.neighbours[direction(side)]
        assert t is not start # prevent infinite loops
    return fan


def replace_kvertex(t, v, newv, now, direction, queue, immediate):
    """Replace kinetic vertex at incident triangles

    Returns fan of triangles that were replaced
    """
    logging.debug("replace_kvertex, start at: {0} [{1}] dir: {2}".format(id(t), t.info, direction))
    fan = []
    first = True
    while t is not None:
        logging.debug(" @ {} [{}]".format(id(t), t.info))
        logging.debug(t.event)
        if t.event is not None and near_zero(now - t.event.time):
            logging.debug(near_zero(now - t.event.time))
            logging.debug("""

            SAME SAME TIME... ARE WE PARALLEL?

            """)
            if t.event.tp == 'flip':
                logging.debug(t.neighbours[t.event.side[0]]) # -- can have become None
                logging.error('Error with current event -- as we do not handle flip now, we run the risk of inconsistency -- in fan: {0} $'.format(t.event))

        side = t.vertices.index(v)
        fan.append(t)
        t.vertices[side] = newv
        logging.debug(
            "Placed vertex #{} [{}] (inf fast? {}) at side {} of triangle {} [{}]".format(
                id(newv),
                newv.info,
                newv.inf_fast,
                side,
                id(t),
                t.info
                ))
        if newv.inf_fast and t.event is not None:  # infinitely fast
            queue.discard(t.event)
            if t.event in immediate:
                immediate.remove(t.event)
        else:  # vertex moves with normal speed
            replace_in_queue(t, now, queue, immediate)
        t = t.neighbours[direction(side)]
    return tuple(fan)


def replace_in_queue(t, now, queue, immediate):
    """Replace event for a triangle in the queue """
    if t.event is not None:
        queue.discard(t.event)
        if t.event in immediate:
            immediate.remove(t.event)
    else:
        logging.debug(
            "triangle #{0} without event not removed from queue".format(
                id(t)))
    e = compute_collapse_time(t, now)
    if e is not None:
        logging.debug("new event in queue {}".format(e))
        queue.add(e)
    else:
        logging.debug("no new events".format(e))
        return


def update_circ(v_left, v_right, now):
    """Update neighbour list of kinetic vertices (this is a circular list going
    around the wavefront edges

    Note that for a vertex often 2 sides need to be updated.
    """
    if v_left is not None:
        logging.debug("update_circ at right of #{} [{}] lies #{} [{}]".format(
                                                  id(v_left),
                                                  v_left.info,
                                                  id(v_right),
                                                  v_right.info if v_right is not None else ""))
        v_left.right = v_right, now
    if v_right is not None:
        logging.debug("update_circ at left  of #{} [{}] lies #{} [{}]".format(
                                                  id(v_right),
                                                  v_right.info,
                                                  id(v_left),
                                                  v_left.info if v_left is not None else ""))
        v_right.left = v_left, now


def schedule_immediately(tri, now, queue, immediate):
    """Schedule a triangle for immediate processing

    Computes a new event for the triangle, where we only check
    how the triangle collapses (look at side lengths at time=now);
    it is assumed that the triangle collapses (neighbor relations say so)!

    The original event is removed from the event queue and the new
    event is added to the immediate queue.

    """
    logging.debug("Scheduling triangle [{}] for direct collapse".format(tri.info))
    queue.discard(tri.event)
    if tri.event in immediate:
        immediate.remove(tri.event)
    E = compute_new_edge_collapse_event(tri, now)
    tri.event = E
    if tri.neighbours.count(None) == 3:
        tri.event.side = list(range(3))
    assert len(tri.event.side) > 0
    immediate.append(tri.event)
