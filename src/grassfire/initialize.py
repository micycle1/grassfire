import logging

from triangle.delaunay.iter import RegionatedTriangleIterator, StarEdgeIterator, Edge
from triangle.delaunay.tds import cw, ccw, orient2d

from grassfire.primitives import Skeleton, SkeletonNode
from grassfire.primitives import InfiniteVertex, KineticTriangle, KineticVertex
from grassfire.line2d import WaveFront, WaveFrontIntersector


def rotate_until_not_in_candidates(t, v, direction, candidates):
    """Rotate around a vertex, starting from triangle t until a triangle is
    found that is not in the candidates list. This triangle is then returned.
    """
    seen = set()
    while t is not None and t not in seen:
        seen.add(t)
        side = t.vertices.index(v)
        t = t.neighbours[direction(side)]
        if t not in candidates:
            return t
    return None


def make_support_line(edge):
    """Given an edge object, make a wavefront support line for this edge.

    Returns:
        WaveFront, if edge is constrained
        None, if edge is not constrained
    """
    if edge.constrained:
        start = (edge.segment[0].x, edge.segment[0].y)
        end = (edge.segment[1].x, edge.segment[1].y)
        return WaveFront(start, end)
    return None


def split_star(v):
    """Splits the edges of the star in groups by constrained edges."""
    around = [e for e in StarEdgeIterator(v)]
    groups = []
    group = []
    for edge in around:
        t, s = edge.triangle, edge.side
        group.append(edge)
        # if the edge ahead is constrained, we make a new group
        if Edge(t, ccw(s)).constrained:
            groups.append(group)
            group = []
    if group:
        groups.append(group)

    # if we have only one group, we do not have to change the first and the last
    if len(groups) <= 1:
        return groups

    # merge first and last group
    # this is necessary when the group does not start at a constrained edge
    edge = groups[0][0]
    if not edge.triangle.constrained[cw(edge.side)]:
        last = groups.pop()
        last.extend(groups[0])
        groups[0] = last

    # -- post condition checks
    # the first and last triangle in the group have a constrained
    # and the rest of the triangles in the middle of every group do not
    for group in groups:
        first, last = group[0], group[-1]
        assert first.triangle.constrained[cw(first.side)]
        assert last.triangle.constrained[ccw(last.side)]
        for middle in group[1:-1]:
            assert not middle.triangle.constrained[cw(middle.side)]
            assert not middle.triangle.constrained[ccw(middle.side)]

    return groups


def init_skeleton(dt):
    """Initialize a data structure that can be used for making the straight
    skeleton.
    """
    skel = Skeleton()

    # make skeleton nodes: every triangulation vertex becomes a skeleton node
    nodes = {}
    avg_x = 0.0
    avg_y = 0.0
    for v in dt.vertices:
        if v.is_finite:
            nodes[v] = SkeletonNode(pos=(v.x, v.y), step=-1, info=v.info)
            avg_x += v.x / len(dt.vertices)
            avg_y += v.y / len(dt.vertices)

    centroid = InfiniteVertex()
    centroid.origin = (avg_x, avg_y)

    ktriangles = []  # all kinetic triangles

    internal_triangles = set()
    for _, depth, triangle in RegionatedTriangleIterator(dt):
        if depth == 1:
            # FIXME: why not put the depth here as identifier into the triangle
            # holes can then be separated from others (and possibly calculated in parallel)
            internal_triangles.add(triangle)

    triangle2ktriangle = {}
    for idx, t in enumerate(dt.triangles, start=1):
        k = KineticTriangle()
        k.info = idx
        triangle2ktriangle[t] = k
        ktriangles.append(k)
        k.internal = t in internal_triangles  # whether triangle is internal to a polygon
    del idx

    link_around = []
    # set up properly the neighbours of all kinetic triangles
    # blank out the neighbour, if a side is constrained
    unwanted = []
    for t in dt.triangles:
        k = triangle2ktriangle[t]

        for i in range(3):
            edge = Edge(t, i)
            k.wavefront_support_lines[i] = make_support_line(edge)

        for j, n in enumerate(t.neighbours):
            # set neighbour pointer to None if constrained side
            if t.constrained[j]:
                continue
            # skip linking to non-existing triangle
            if n is None or n.vertices[2] is None:
                unwanted.append(k)
                continue
            k.neighbours[j] = triangle2ktriangle[n]

    # make kinetic vertices and link them to kinetic triangles
    # also make sure that every kinetic vertex is related to a skeleton node
    kvertices = []
    ct = 0
    for v in dt.vertices:
        assert v.is_finite, "infinite vertex found"

        groups = split_star(v)
        if len(groups) == 1:
            raise NotImplementedError(
                "not yet dealing with PSLG in initial conversion"
            )

        for group in groups:
            first, last = group[0], group[-1]

            # compute turn type at vertex
            tail, mid1 = Edge(last.triangle, ccw(last.side)).segment
            mid2, head = Edge(first.triangle, cw(first.side)).segment
            assert mid1 is mid2
            turn = orient2d((tail.x, tail.y), (mid1.x, mid1.y), (head.x, head.y))
            # left : + [ = ccw ]
            # straight : 0.
            # right : - [ = cw ]

            if turn < 0:
                turn_type = "RIGHT - REFLEX"
            elif turn > 0:
                turn_type = "LEFT - CONVEX"
            else:
                turn_type = "STRAIGHT"

            right = triangle2ktriangle[first.triangle].wavefront_support_lines[cw(first.side)]
            left = triangle2ktriangle[last.triangle].wavefront_support_lines[ccw(last.side)]

            assert left is not None
            assert right is not None

            intersector = WaveFrontIntersector(left, right)
            bi = intersector.get_bisector()

            ur = right.line
            ul = left.line

            ct += 1
            kv = KineticVertex()
            kv.turn = turn_type
            kv.info = ct
            kv.origin = (v.x, v.y)
            kv.velocity = bi
            kv.start_node = nodes[v]
            kv.starts_at = 0

            kv.ul = ul
            kv.ur = ur

            kv.wfl = left
            kv.wfr = right

            for edge in group:
                ktriangle = triangle2ktriangle[edge.triangle]
                ktriangle.vertices[edge.side] = kv
                kv.internal = ktriangle.internal

            kvertices.append(kv)

            # link vertices to each other in circular list
            link_around.append(
                ((last.triangle, cw(last.side)), kv, (first.triangle, ccw(first.side)))
            )

    # link vertices in circular list
    for left, kv, right in link_around:  # left is cw, right is ccw
        assert left is not None
        assert right is not None

        cwv = triangle2ktriangle[left[0]].vertices[left[1]]
        kv.left = cwv, 0

        ccwv = triangle2ktriangle[right[0]].vertices[right[1]]
        kv.right = ccwv, 0

    for left, kv, right in link_around:  # left is cw, right is ccw
        assert kv.left.wfr is kv.wfl, "{} vs\n {}".format(kv.left.wfr, kv.wfl)
        assert kv.wfr is kv.right.wfl
        assert kv.is_stopped is False

    # -- copy infinite vertices into the kinetic triangles
    # make dico of infinite vertices (lookup by coordinate value)
    infinites = {}
    for t in triangle2ktriangle:
        for v in t.vertices:
            if v is not None and not v.is_finite:
                infv = InfiniteVertex()
                infv.origin = (v[0], v[1])
                infinites[(v[0], v[1])] = infv
    assert len(infinites) == 3

    # link infinite triangles to the infinite vertex
    for (t, kt) in triangle2ktriangle.items():
        for i, v in enumerate(t.vertices):
            if v is not None and not v.is_finite:
                kt.vertices[i] = infinites[(v[0], v[1])]

    # there are 3 infinite triangles that are supposed to be removed
    # these triangles were already stored in the unwanted list
    remove = []
    for kt in ktriangles:
        if [isinstance(v, InfiniteVertex) for v in kt.vertices].count(True) == 2:
            remove.append(kt)
    assert len(remove) == 3
    assert len(unwanted) == 3
    assert remove == unwanted

    # remove the 3 unwanted triangles and link their neighbours together
    link = []
    for kt in unwanted:
        v = kt.vertices[kt.neighbours.index(None)]
        assert isinstance(v, KineticVertex)

        neighbour_cw = rotate_until_not_in_candidates(kt, v, cw, unwanted)
        neighbour_ccw = rotate_until_not_in_candidates(kt, v, ccw, unwanted)

        side_cw = ccw(neighbour_cw.vertices.index(v))
        side_ccw = cw(neighbour_ccw.vertices.index(v))

        link.append((neighbour_cw, side_cw, neighbour_ccw))
        link.append((neighbour_ccw, side_ccw, neighbour_cw))

    for ngb, side, new_ngb in link:
        ngb.neighbours[side] = new_ngb

    for kt in unwanted:
        kt.vertices = [None, None, None]
        kt.neighbours = [None, None, None]
        ktriangles.remove(kt)

    # replace the infinite vertices by one point in the center of the PSLG
    # (this could be the origin (0,0) if we would scale input to [-1,1] range)
    for kt in ktriangles:
        for i, v in enumerate(kt.vertices):
            if isinstance(v, InfiniteVertex):
                kt.vertices[i] = centroid

    assert check_ktriangles(ktriangles)

    ktriangles.sort(
        key=lambda t: (t.vertices[0].origin[1], t.vertices[0].origin[0])
    )
    skel.sk_nodes = list(nodes.values())
    skel.triangles = ktriangles
    skel.vertices = kvertices

    return skel


def internal_only_skeleton(skel):
    """Filter a skeleton and only maintain internal elements.

    Internal means enclosed by a set of boundaries, this way the skeleton
    to the interior of a polygon will be skeletonized.
    """
    new = Skeleton()
    new.sk_nodes = skel.sk_nodes[:]
    new.triangles = [t for t in skel.triangles if t.internal]
    new.vertices = [v for v in skel.vertices if v.internal]
    return new


def check_ktriangles(L, now=0):
    """Check whether kinetic triangles are all linked up properly."""
    valid = True

    # check if neighbours are properly linked
    for ktri in L:
        if ktri.stops_at is not None:
            continue

        for ngb in ktri.neighbours:
            if ngb is not None and ktri not in ngb.neighbours:
                logging.warning(
                    "Non neighbouring triangles: %s and %s",
                    id(ktri),
                    id(ngb),
                )
                valid = False

        for v in ktri.vertices:
            if ktri.is_finite:
                ok = (
                    (v.starts_at <= now and v.stops_at is not None and v.stops_at > now)
                    or (v.starts_at <= now and v.stops_at is None)
                )
                if not ok:
                    logging.warning(
                        "Triangle %s has invalid kinetic vertex %s for time=%s (starts_at=%s stops_at=%s)",
                        id(ktri),
                        id(v),
                        now,
                        v.starts_at,
                        v.stops_at,
                    )
                    valid = False

    # check if the sides of a triangle share the correct vertex at begin / end
    if False:  # FIXME: enable!!!
        for ktri in L:
            for i in range(3):
                ngb = ktri.neighbours[i]
                if ngb is not None:
                    j = ngb.neighbours.index(ktri)
                    if not ngb.vertices[cw(j)] is ktri.vertices[ccw(i)]:
                        print("something wrong with vertices 1")
                        valid = False
                    if not ngb.vertices[ccw(j)] is ktri.vertices[cw(i)]:
                        print("something wrong with vertices 2")
                        valid = False

    # FIXME: check orientation of triangles ????
    # -- could be little difficult with initial needle triangles at terminal
    # vertices of PSLG
    return valid