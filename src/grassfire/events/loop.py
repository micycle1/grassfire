import logging
from collections import deque

from tri.delaunay.tds import Edge
from grassfire.ordered_sequence import OrderedSequence

from grassfire.calc import near_zero
from grassfire.collapse import compute_collapse_time, find_gt

from grassfire.events.edge import (
    handle_edge_event,
    handle_edge_event_1side,
    handle_edge_event_3sides,
)
from grassfire.events.flip import handle_flip_event
from grassfire.events.split import handle_split_event
from grassfire.events.check import check_active_triangles_orientation, check_bisectors

from grassfire.inout import interactive_visualize, visualize


def choose_next_event(queue):
    """Choose a next event from the queue.

    NOTE: Event selection order is currently driven by the queue ordering.
    """
    it = iter(queue)
    evt = next(it)
    queue.remove(evt)
    return evt


def log_queue_content(step, immediate, queue):
    for i, e in enumerate(immediate):
        logging.debug("{0:5d} {1}".format(i, e))
    for i, e in enumerate(queue):
        logging.debug("{0:5d} {1}".format(i, e))
        if i >= 20:
            break
    if len(queue) >= 20:
        logging.debug("... skipping display of {} events".format(len(queue) - 20))
def event_loop(queue, skel, pause=False, stop_after=0, make_video=False, video_digits=3):
    """The main event loop.

    Args:
        queue: Event queue
        skel: Skeleton structure
        pause: Whether to pause for interactive visualization
        stop_after: Stop after this many steps (0 = no limit) - for testing/debugging
        make_video: Whether to generate video frames - for testing/debugging
        video_digits: Number of decimal digits for video timing - for testing/debugging
    """
    if stop_after != 0:
        logging.debug("Stopping for the first time after step#{}".format(stop_after))

    if make_video:
        from grassfire.tests2.test_event_loop_debugging import make_frames

    NOW = 0.0
    step = 0

    if pause:
        logging.getLogger().setLevel(logging.DEBUG)
        interactive_visualize(queue, skel, step, NOW)

    immediate = deque([])

    for i, e in enumerate(immediate):
        logging.debug("{0:5d} {1}".format(i, e))
    for i, e in enumerate(queue):
        logging.debug("{0:5d} {1}".format(i, e))

    if make_video:
        make_frames(NOW, video_digits, skel, queue, immediate)

    check_bisectors(skel, 0.0)
    check_active_triangles_orientation(skel.triangles, 0)

    guard = 0
    while queue or immediate:
        guard += 1
        if guard > 50000:
            raise ValueError("loop with more than 50_000 events stopped")

        step += 1
        if pause:
            log_queue_content(step, immediate, queue)

        if immediate:
            evt = immediate.popleft()
        else:
            evt = choose_next_event(queue)
            NOW = evt.time

        logging.debug(
            "About to handle event %s %s %s [%s] at time %.28g",
            evt.tp,
            evt.triangle.type,
            id(evt.triangle),
            evt.triangle.info,
            evt.time,
        )

        if pause and step >= stop_after:

            def check_direct(event):
                """Check the direct neighbours whether they will collapse."""
                current = event.triangle
                seen = {current}
                visit = [current]
                also = []
                while visit:
                    current = visit.pop()
                    for n in current.neighbours:
                        if n is not None and n.event is not None and n not in seen:
                            seen.add(n)
                            if near_zero(n.event.time - event.time):
                                also.append(n)
                                visit.append(n)
                if also:
                    logging.debug([n.info for n in also])

            check_direct(evt)
            with open("/tmpfast/current_event.wkt", "w") as fh:
                fh.write(
                    "pos;wkt;evttype;evttime;tritype;id;n0;n1;n2;finite;info;wavefront_directions\n"
                )
                fh.write(
                    "{0};{1};{2};{3};{4};{5};{6};{7};{8};{9};{10}\n".format(
                        0,
                        evt.triangle.visualize_at(NOW),
                        evt.tp,
                        evt.time,
                        evt.triangle.type,
                        id(evt.triangle),
                        id(evt.triangle.neighbours[0]),
                        id(evt.triangle.neighbours[1]),
                        id(evt.triangle.neighbours[2]),
                        evt.triangle.is_finite,
                        evt.triangle.info,
                        evt.triangle.wavefront_directions,
                    )
                )

            interactive_visualize(queue, skel, step, NOW)
        if evt.triangle.stops_at is not None:
            logging.warning("Already stopped %s, but still queued", id(evt.triangle))
            continue

        if evt.tp == "edge":
            if len(evt.side) == 3:
                handle_edge_event_3sides(evt, step, skel, queue, immediate)
            elif len(evt.side) == 1 and evt.triangle.type == 3:
                handle_edge_event_1side(
                    evt,
                    step,
                    skel,
                    queue,
                    immediate,
                    pause and step >= stop_after,
                )
            elif len(evt.side) == 2:
                raise ValueError(
                    "Impossible configuration, triangle [{}] with 2 sides collapsing cannot happen".format(
                        evt.triangle.info
                    )
                )
            else:
                handle_edge_event(
                    evt,
                    step,
                    skel,
                    queue,
                    immediate,
                    pause and step >= stop_after,
                )
        elif evt.tp == "flip":
            handle_flip_event(evt, step, skel, queue, immediate)
        elif evt.tp == "split":
            handle_split_event(evt, step, skel, queue, immediate, pause and step >= stop_after)

        logging.debug("=" * 80)

        if make_video:
            make_frames(NOW, video_digits, skel, queue, immediate)

    if pause:
        visualize(queue, skel, NOW)

    if make_video:
        make_frames(NOW, video_digits, skel, queue, immediate)

    if pause:
        interactive_visualize(queue, skel, step, NOW)
    not_stopped_tris = []
    for tri in skel.triangles:
        if all(v.internal for v in tri.vertices) and tri.stops_at is None:
            not_stopped_tris.append(tri.info)
    if not_stopped_tris:
        raise ValueError("triangles not stopped at end: {}".format(not_stopped_tris))

    return NOW


def compare_event_by_time(one, other):
    """Compare two events, first by time, in case they are equal by triangle type
    (first 2-triangle, then 1-triangle, then 0-triangle), as last resort by
    identifier of triangle.
    """
    if one.time < other.time:
        return -1
    elif one.time > other.time:
        return 1
    else:
        if -one.triangle_tp < -other.triangle_tp:
            return -1
        elif -one.triangle_tp > -other.triangle_tp:
            return 1
        else:
            if id(one.triangle) < id(other.triangle):
                return -1
            elif id(one.triangle) > id(other.triangle):
                return 1
            else:
                return 0


def init_event_list(skel):
    """Compute for all kinetic triangles when they will collapse and put them in
    an OrderedSequence, so that events are ordered properly for further processing.
    """
    q = OrderedSequence(cmp=compare_event_by_time)
    for tri in skel.triangles:
        res = compute_collapse_time(tri, 0, find_gt)
        if res is not None:
            q.add(res)
    return q
