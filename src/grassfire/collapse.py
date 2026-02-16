import bisect
import logging
import math

from tri.delaunay.tds import cw, ccw, Edge

from grassfire.calc import get_unique_times, near_zero
from grassfire.inout import output_edges_at_T, output_triangles_at_T, output_vertices_at_T
from grassfire.primitives import Event, InfiniteVertex
from grassfire.vectorops import add, dot, mul, sub, norm
from predicates import orient2d_xy as orient2d


def find_gt(a, x):
    """Find value greater than x, ignores None values."""
    a = [x for x in a if x is not None]
    a = [v for v in a if not near_zero(v - x)]
    L = sorted(a)
    i = bisect.bisect_right(L, x)
    if i != len(a):
        return L[i]
    return None


def find_gte(a, x):
    """Find greater than or equal to x, ignores None values."""
    L = sorted([x for x in a if x is not None])
    i = bisect.bisect_left(L, x)
    if i != len(L):
        return L[i]
    return None


def vertex_crash_time(org, dst, apx):
    """Returns time when vertex crashes on edge.

    This method assumes that velocity of wavefront is unit speed.

    Input:
        org, dst: kinetic vertices incident with wavefront edge
        apx: kinetic vertex opposite of wavefront edge
    """
    Mv = tuple(sub(apx.origin, org.origin))
    assert org.ur is not None
    assert org.ur == dst.ul, "#{} #{} :: {} vs {}".format(id(org), id(dst), org.ur, dst.ul)

    n = tuple(org.ur.w)  # was: org.ur

    s = apx.velocity
    dist_v_e = dot(Mv, n)

    s_proj = dot(s, n)

    denom = 1.0 - s_proj
    if not near_zero(denom):
        t_v = dist_v_e / denom
        return t_v
    return None


def area_collapse_times(o, d, a):
    coeff = area_collapse_time_coeff(o, d, a)
    solution = solve_quadratic(coeff[0], coeff[1], coeff[2])
    solution.sort()
    return solution


def compute_event_0triangle(tri, now, sieve):
    assert tri.neighbours.count(None) == 0
    o, d, a = tri.vertices

    times_area_collapse = area_collapse_times(o, d, a)
    for time in times_area_collapse:
        if near_zero(abs(time - now)):
            dists = [
                d.distance2_at(a, now),
                a.distance2_at(o, now),
                o.distance2_at(d, now),
            ]

            indices = []
            for i, _ in enumerate([math.sqrt(_) for _ in dists]):
                if near_zero(_):
                    indices.append(i)
            if len(indices) == 1:
                side = indices[0]
                return Event(when=now, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)
            elif len(indices) == 3:
                raise ValueError("0-triangle collapsing to point")
            else:
                largest_dist = max(dists)
                side = dists.index(largest_dist)
                return Event(when=now, tri=tri, side=(side,), tp="flip", tri_tp=tri.type)

    times_edge_collapse = [
        collapse_time_edge(o, d),
        collapse_time_edge(d, a),
        collapse_time_edge(a, o),
    ]

    dists = [
        o.distance2_at(d, times_edge_collapse[0]),
        d.distance2_at(a, times_edge_collapse[1]),
        a.distance2_at(o, times_edge_collapse[2]),
    ]

    indices = []
    for i, _ in enumerate([math.sqrt(_) for _ in dists]):
        if near_zero(_):
            indices.append(i)
    t_e_c = []
    for i in indices:
        t_e_c.append(times_edge_collapse[i])
    time_edge_collapse = sieve(t_e_c, now)
    time_area_collapse = sieve(times_area_collapse, now)

    if time_edge_collapse is None and time_area_collapse is None:
        return None
    elif time_edge_collapse is not None and time_area_collapse is not None:
        if near_zero(abs(time_area_collapse - time_edge_collapse)):
            time = time_edge_collapse
            dists = [
                d.distance2_at(a, time),
                a.distance2_at(o, time),
                o.distance2_at(d, time),
            ]
            zeros = [near_zero(d - min(dists)) for d in dists]
            sides_collapse = zeros.count(True)
            if sides_collapse == 3:
                return Event(when=time, tri=tri, side=(0, 1, 2), tp="edge", tri_tp=tri.type)
            elif sides_collapse == 1:
                side = zeros.index(True)
                return Event(when=time, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)
            else:
                time = time_area_collapse
                dists = [
                    d.distance2_at(a, time),
                    a.distance2_at(o, time),
                    o.distance2_at(d, time),
                ]
                largest_dist = max(dists)
                side = dists.index(largest_dist)
                return Event(when=time, tri=tri, side=(side,), tp="flip", tri_tp=tri.type)
        elif time_area_collapse < time_edge_collapse:
            time = time_area_collapse
            dists = [
                d.distance2_at(a, time),
                a.distance2_at(o, time),
                o.distance2_at(d, time),
            ]
            largest_dist = max(dists)
            side = dists.index(largest_dist)
            return Event(when=time, tri=tri, side=(side,), tp="flip", tri_tp=tri.type)
        elif time_edge_collapse is not None:
            time = time_edge_collapse
            dists = [
                d.distance2_at(a, time),
                a.distance2_at(o, time),
                o.distance2_at(d, time),
            ]
            zeros = [near_zero(_) for _ in dists]
            sides_collapse = zeros.count(True)
            if sides_collapse == 3:
                return Event(when=time, tri=tri, side=(0, 1, 2), tp="edge", tri_tp=tri.type)
            elif sides_collapse == 1:
                side = zeros.index(True)
                return Event(when=time, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)
            else:
                raise ValueError("can this happen?")

    else:
        if time_edge_collapse is not None:
            time = time_edge_collapse
            dists = [
                d.distance2_at(a, time),
                a.distance2_at(o, time),
                o.distance2_at(d, time),
            ]
            zeros = [near_zero(_) for _ in dists]
            sides_collapse = zeros.count(True)
            if sides_collapse == 3:
                return Event(when=time, tri=tri, side=(0, 1, 2), tp="edge", tri_tp=tri.type)
            elif sides_collapse == 1:
                side = zeros.index(True)
                return Event(when=time, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)
            else:
                raise ValueError(
                    "0 triangle with 2 or 0 side collapse,"
                    "while edge collapse time computed?"
                )
        elif time_area_collapse is not None:
            time = time_area_collapse
            dists = [
                d.distance2_at(a, time),
                a.distance2_at(o, time),
                o.distance2_at(d, time),
            ]
            largest_dist = max(dists)
            side = dists.index(largest_dist)
            return Event(when=time, tri=tri, side=(side,), tp="flip", tri_tp=tri.type)
        else:
            raise ValueError("problem!!!")


def compute_event_1triangle(tri, now, sieve):
    assert tri.neighbours.count(None) == 1
    wavefront_side = tri.neighbours.index(None)

    o, d, a = tri.vertices
    ow, dw, aw = [
        tri.vertices[ccw(wavefront_side)],
        tri.vertices[cw(wavefront_side)],
        tri.vertices[wavefront_side],
    ]

    times_vertex_crash = [vertex_crash_time(ow, dw, aw)]
    for time in times_vertex_crash:
        if time is None:
            continue
        if near_zero(abs(time - now)):
            time = now

            dists = [
                d.distance2_at(a, now),
                a.distance2_at(o, now),
                o.distance2_at(d, now),
            ]
            indices = []
            for i, _ in enumerate([math.sqrt(_) for _ in dists]):
                if near_zero(_):
                    indices.append(i)
            if len(indices) == 1:
                side = indices[0]
                return Event(when=now, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)

            dists = [
                math.sqrt(d.distance2_at(a, time)),
                math.sqrt(a.distance2_at(o, time)),
                math.sqrt(o.distance2_at(d, time)),
            ]
            longest_side = dists.index(max(dists))
            tp = "split" if longest_side == wavefront_side else "flip"
            return Event(when=time, tri=tri, side=(longest_side,), tp=tp, tri_tp=tri.type)

    time_vertex_crash = sieve(times_vertex_crash, now)

    time_area_collapse = sieve(area_collapse_times(o, d, a), now)

    time_edge_collapse = sieve([collapse_time_edge(ow, dw)], now)

    if time_edge_collapse is None and time_vertex_crash is None:
        time = sieve(
            solve_quadratic(*area_collapse_time_coeff(*tri.vertices)),
            now,
        )
        if time is None:
            return None
        elif near_zero(time - now) is True:
            return Event(
                when=now,
                tri=tri,
                side=(wavefront_side,),
                tp="split",
                tri_tp=tri.type,
            )
        else:
            dists = [
                d.distance2_at(a, time) if tri.neighbours[0] is not None else -1,
                a.distance2_at(o, time) if tri.neighbours[1] is not None else -1,
                o.distance2_at(d, time) if tri.neighbours[2] is not None else -1,
            ]
            logging.debug(" {}".format(dists))
            sides = (dists.index(max(dists)),)
            return Event(when=time, tri=tri, side=sides, tp="flip", tri_tp=tri.type)

    elif time_edge_collapse is None and time_vertex_crash is not None:
        logging.debug(" case B, time vertex crash " + str(time_vertex_crash))

        if time_area_collapse is not None and time_area_collapse < time_vertex_crash:
            time = time_area_collapse
        else:
            time = time_vertex_crash

        dists = [
            math.sqrt(d.distance2_at(a, time)),
            math.sqrt(a.distance2_at(o, time)),
            math.sqrt(o.distance2_at(d, time)),
        ]
        unique_dists = [near_zero(_ - max(dists)) for _ in dists]
        unique_max_dists = unique_dists.count(True)
        logging.debug("uniq max dists (needs to be 1 for split) -- " + str(unique_max_dists))
        logging.debug("wavefront side -- " + str(wavefront_side))
        logging.debug("longest side -- " + str(dists.index(max(dists))))
        longest = []
        for i, _ in enumerate(dists):
            if near_zero(_ - max(dists)):
                longest.append(i)
        if wavefront_side in longest and len(longest) == 1:
            return Event(
                when=time_vertex_crash,
                tri=tri,
                side=(wavefront_side,),
                tp="split",
                tri_tp=tri.type,
            )
        else:
            zeros = [near_zero(_) for _ in dists]
            sides_collapse = zeros.count(True)
            if sides_collapse == 1:
                sides = (dists.index(min(dists)),)  # shortest side
                return Event(when=time, tri=tri, side=sides, tp="edge", tri_tp=tri.type)
            else:
                sides = (dists.index(max(dists)),)  # longest side
                return Event(when=time, tri=tri, side=sides, tp="flip", tri_tp=tri.type)

    elif time_edge_collapse is not None and time_vertex_crash is None:
        return Event(
            when=time_edge_collapse,
            tri=tri,
            side=(wavefront_side,),
            tp="edge",
            tri_tp=tri.type,
        )

    elif time_edge_collapse is not None and time_vertex_crash is not None:
        if time_edge_collapse < time_vertex_crash or time_edge_collapse == time_vertex_crash:
            time = time_edge_collapse

            dists_squared = [
                d.distance2_at(a, time),
                a.distance2_at(o, time),
                o.distance2_at(d, time),
            ]

            sides = [dists_squared.index(min(dists_squared))]
            assert len(sides) == 1
            tp = "edge"
            return Event(when=time, tri=tri, side=sides, tp=tp, tri_tp=tri.type)

        elif time_vertex_crash < time_edge_collapse:
            time = time_vertex_crash
            dists = [
                math.sqrt(d.distance2_at(a, time)),
                math.sqrt(a.distance2_at(o, time)),
                math.sqrt(o.distance2_at(d, time)),
            ]
            logging.debug("dists {}".format(dists))
            zeros = [near_zero(dist) for dist in dists]
            logging.debug("zeros {}".format(zeros))
            if True in zeros and zeros.count(True) == 1:
                tp = "edge"
                sides = (zeros.index(True),)
            elif True in zeros and zeros.count(True) == 3:
                tp = "edge"
                sides = list(range(3))
            else:
                max_dist = max(dists)
                max_dist_side = dists.index(max_dist)
                tp = "split" if tri.neighbours[max_dist_side] is None else "flip"
                sides = (max_dist_side,)
            return Event(when=time, tri=tri, side=sides, tp=tp, tri_tp=tri.type)
        else:
            raise NotImplementedError("Problem, unforeseen configuration")

    raise NotImplementedError("Problem, unforeseen configuration")


def compute_event_2triangle(tri, now, sieve):
    assert tri.neighbours.count(None) == 2
    o, d, a = tri.vertices
    times = []
    if tri.neighbours[2] is None:
        times.append(collapse_time_edge(o, d))
    if tri.neighbours[0] is None:
        times.append(collapse_time_edge(d, a))
    if tri.neighbours[1] is None:
        times.append(collapse_time_edge(a, o))

    times = get_unique_times(times)
    time = sieve(times, now)
    if time is None:
        time = sieve(area_collapse_times(o, d, a), now)

    if time is not None:
        dists = [
            math.sqrt(d.distance2_at(a, time)),
            math.sqrt(a.distance2_at(o, time)),
            math.sqrt(o.distance2_at(d, time)),
        ]
        dists = [d - min(dists) for d in dists]
        logging.debug("distances at time = {1}: {0}".format(dists, time))
        zeros = [near_zero(dist) for dist in dists]
        logging.debug("near_zero = {}".format(zeros))
        sides_collapse = zeros.count(True)
        if sides_collapse == 3:
            sides = tuple(range(3))
            return Event(when=time, tri=tri, side=sides, tp="edge", tri_tp=tri.type)
        elif sides_collapse == 2:
            raise ValueError("This is not possible with this type of triangle [{}]".format(tri.info))
        elif sides_collapse == 1:
            side = dists.index(min(dists))
            return Event(when=time, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)
        elif sides_collapse == 0:
            return None
    else:
        return None


def compute_event_3triangle(tri, now, sieve):
    a, o, d = tri.vertices
    t_e_c = [
        collapse_time_edge(o, d),
        collapse_time_edge(d, a),
        collapse_time_edge(a, o),
    ]

    dists = [
        o.distance2_at(d, t_e_c[0]),
        d.distance2_at(a, t_e_c[1]),
        a.distance2_at(o, t_e_c[2]),
    ]

    indices = []
    for i, _ in enumerate([math.sqrt(_) for _ in dists]):
        if near_zero(_):
            indices.append(i)

    assert tri.neighbours.count(None) == 3

    time_edge_collapse = sieve(t_e_c, now)
    time_area_collapse = sieve(area_collapse_times(o, d, a), now)

    if time_edge_collapse:
        sides = tuple(indices)
        if len(sides) == 0:
            logging.error(
                "3-triangle: override # of sides collapsing -- instead of 0!, we collapse all 3"
            )
            sides = list(range(3))
        if len(sides) == 2:
            logging.warning(
                "3-triangle: override # of sides collapsing -- instead of 2, we collapse all 3"
            )
            sides = list(range(3))
        return Event(when=time_edge_collapse, tri=tri, side=sides, tp="edge", tri_tp=tri.type)
    elif time_area_collapse:
        logging.error("3-triangle: using area collapse time as fall back not to miss out")
        sides = list(range(3))
        return Event(when=time_area_collapse, tri=tri, side=sides, tp="edge", tri_tp=tri.type)
    else:
        return None


def compute_event_inftriangle(tri, now, sieve):
    for inf_idx, v in enumerate(tri.vertices):
        if isinstance(v, InfiniteVertex):
            break
    side = inf_idx
    o, d, a = tri.vertices[cw(side)], tri.vertices[ccw(side)], tri.vertices[side]
    if tri.neighbours[side] is None:  # wavefront edge on the hull that collapses
        assert tri.type == 1, tri.type
        time = find_gt([collapse_time_edge(o, d)], now)
        logging.debug("time of closest approach {}".format(time))
        if time:
            dist = o.distance2_at(d, time)
            if near_zero(dist):
                return Event(when=time, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)
            else:
                return None
    else:
        time = sieve(area_collapse_times(o, d, a), now)
        logging.debug("time = {}".format(time))
        if time:
            dist = o.distance2_at(d, time)
            if near_zero(dist):
                return Event(when=time, tri=tri, side=(side,), tp="edge", tri_tp=tri.type)
            else:
                tp = "flip"
                dists = []
                for func in cw, ccw:
                    start, end = Edge(tri, func(side)).segment
                    dists.append(start.distance2_at(end, time))
                idx = dists.index(min(dists))
                min_dist_side = [cw, ccw][idx](side)
                return Event(when=time, tri=tri, side=(min_dist_side,), tp=tp, tri_tp=tri.type)
    return None


def compute_collapse_time(tri, now=0, sieve=find_gte):
    """Computes Event that represents how a triangle collapses at a given time."""
    event = None
    if tri.stops_at is not None:
        return event

    if tri.is_finite:
        logging.debug("=-=-= finite triangle #{} [{}]=-=-= ".format(id(tri), tri.info))
        tp = tri.type
        if tp == 0:
            event = compute_event_0triangle(tri, now, sieve)
        elif tp == 1:
            event = compute_event_1triangle(tri, now, sieve)
        elif tp == 2:
            event = compute_event_2triangle(tri, now, sieve)
        elif tp == 3:
            event = compute_event_3triangle(tri, now, sieve)

        if event is not None and all(not v.inf_fast for v in tri.vertices) is True:
            assert event is not None
            verts = [v.position_at(((event.time - now) * 0.5) + now) for v in tri.vertices]

            if orient2d(verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1]) < 0:
                logging.warning(
                    "TRIANGLE TAKES POSSIBLY WRONG ORIENTATION -- we may miss out on handling an event -- "
                    "triangle turns wrong way"
                )

        if event is None:
            verts = [v.position_at(now + 10) for v in tri.vertices]

            if orient2d(verts[0][0], verts[0][1], verts[1][0], verts[1][1], verts[2][0], verts[2][1]) < 0:
                logging.error(
                    "TRIANGLE HAS WRONG ORIENTATION and will not collapse -- we missed out on handling an event -- "
                    "triangle turns wrong way"
                )
    else:
        logging.debug("=-=-= infinite triangle #{} [{}] =-=-=".format(id(tri), tri.info))
        event = compute_event_inftriangle(tri, now, sieve)

    if event is not None:
        tri.event = event

    return event
def compute_new_edge_collapse_event(tri, time):
    """Compute new edge event for triangle that collapse at time.

    Somehow we know that one or more of the edges of this triangle do collapse at this moment.
    """
    o, d, a = tri.vertices
    dists = list(
        map(
            math.sqrt,
            [
                d.distance2_at(a, time),
                a.distance2_at(o, time),
                o.distance2_at(d, time),
            ],
        )
    )
    zeros = [near_zero(dist - min(dists)) for dist in dists]
    sides = []
    for i, zero in enumerate(zeros):
        if zero is True:
            sides.append(i)
    return Event(when=time, tri=tri, side=sides, tp="edge", tri_tp=tri.type)


def collapse_time_edge(v1, v2):
    """Returns the time when the given 2 kinetic vertices are closest to each other."""
    s1 = v1.velocity
    s2 = v2.velocity
    o1 = v1.origin
    o2 = v2.origin
    dv = sub(s1, s2)
    denominator = dot(dv, dv)
    if not near_zero(denominator):
        w0 = sub(o2, o1)
        nominator = dot(dv, w0)
        collapse_time = nominator / denominator
        logging.debug("edge collapse time: " + str(collapse_time))
        return collapse_time
    else:
        logging.debug(str(v1) + "|" + str(v2))
        return -1.0


def solve_quadratic(A, B, C):
    """Get the roots for quadratic equation, defined by A, B and C."""
    if near_zero(A) and not near_zero(B):
        return [-C / B]
    elif near_zero(A) and near_zero(B):
        return []

    T = -B / A  # a + d
    D = C / A  # a*d - b*c
    centre = T * 0.5
    under = 0.25 * pow(T, 2) - D
    if near_zero(under):
        return [centre]
    elif under < 0:
        return []
    else:
        plus_min = math.sqrt(under)
        return [centre - plus_min, centre + plus_min]


def area_collapse_time_coeff(kva, kvb, kvc):
    """Returns coefficients of quadratic in t."""
    pa, shifta = kva.origin, kva.velocity
    pb, shiftb = kvb.origin, kvb.velocity
    pc, shiftc = kvc.origin, kvc.velocity

    xaorig, yaorig = pa[0], pa[1]
    xborig, yborig = pb[0], pb[1]
    xcorig, ycorig = pc[0], pc[1]

    dxa, dya = shifta[0], shifta[1]
    dxb, dyb = shiftb[0], shiftb[1]
    dxc, dyc = shiftc[0], shiftc[1]

    A = dxa * dyb - dxb * dya + dxb * dyc - dxc * dyb + dxc * dya - dxa * dyc
    B = (
        xaorig * dyb
        - xborig * dya
        + xborig * dyc
        - xcorig * dyb
        + xcorig * dya
        - xaorig * dyc
        + dxa * yborig
        - dxb * yaorig
        + dxb * ycorig
        - dxc * yborig
        + dxc * yaorig
        - dxa * ycorig
    )
    C = (
        xaorig * yborig
        - xborig * yaorig
        + xborig * ycorig
        - xcorig * yborig
        + xcorig * yaorig
        - xaorig * ycorig
    )
    ret = (A, B, C)
    return ret


def visualize_collapse(tri, T=0):
    with open("/tmp/bisectors.wkt", "w") as bisector_fh:
        bisector_fh.write("wkt\n")
        for kvertex in tri.vertices:
            p1 = kvertex.position_at(T)
            bi = kvertex.velocity
            bisector_fh.write(
                "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})\n".format(p1, add(p1, bi))
            )

    with open("/tmp/ktris.wkt", "w") as fh:
        output_triangles_at_T([tri], T, fh)

    with open("/tmp/kvertices.wkt", "w") as fh:
        output_vertices_at_T(tri.vertices, T, fh)

    with open("/tmp/wavefront.wkt", "w") as fh:
        try:
            edges = []
            for i, n in enumerate(tri.neighbours):
                if n is None:
                    edges.append(Edge(tri, i))
            output_edges_at_T(edges, T, fh)
        except Exception:
            pass

    with open("/tmp/rays.wkt", "w") as bisector_fh:
        bisector_fh.write("wkt\n")
        for kvertex in tri.vertices:
            p1 = kvertex.origin
            bi = kvertex.velocity
            bineg = mul(bi, -10000.0)
            bipos = mul(bi, 10000.0)
            bisector_fh.write(
                "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})\n".format(
                    add(p1, bipos), add(p1, bineg)
                )
            )