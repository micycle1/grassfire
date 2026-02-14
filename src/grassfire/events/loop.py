import logging

from collections import deque

from oseq import OrderedSequence

from tri.delaunay.tds import Edge

from grassfire.regression import are_residuals_near_zero
from grassfire.calc import near_zero
from grassfire.collapse import compute_collapse_time, find_gt

from grassfire.events.flip import handle_flip_event
from grassfire.events.edge import handle_edge_event, handle_edge_event_3sides, handle_edge_event_1side
from grassfire.events.split import handle_split_event
from grassfire.events.check import check_wavefront_links, check_kinetic_vertices, check_active_triangles, check_bisectors, check_active_triangles_orientation

from grassfire.inout import visualize, interactive_visualize


def choose_next_event(queue):
    """Choose a next event from the queue

    In case we have all flip events, we pick the longest edge to flip

    *Incorrect:*
    In case there are multiple events happening more or less at the same time
    we pick a non-flipping event first. -- This is not okay
    -> we do not handle items in their correct order, and hence this leads to problems
    """
    # it = iter(queue)
    # first = next(it)
    # events = [first]
    # item = events[0]
    # queue.remove(item)
    # return item    

        # # return the first item, sorted on type (non-flipping events first)
        # for item in sorted(events, key=sort_key):
        #     break


    def sort_key(evt):
        """Key for sorting; first by event type, when these equal by time"""
        types = {'flip': 2, 'split': 1, 'edge': 0}
#         types = {'flip': 0, 'split': 1, 'edge': 2}
        # x = evt.triangle.vertices[0].position_at(evt.time)[0]
        # Do sort by event type (first split, than edge, than flip)
        # in case multiple same types, then sort on time, if time is exactly the same, we sort on info (stable id of kinetic triangle)
        # (otherwise we might miss events that should happen)
        return (types[evt.tp], evt.time, -evt.triangle.type, evt.triangle.info)
    it = iter(queue)
    first = next(it)


    if False:
        events = [first]
        for e in it:
            # if e.tp == "flip" and 
            if near_zero(e.time - first.time):
                events.append(e)
            else:
                break

        logging.debug("Multiple Events to pick next event from ...")
        # logging.debug("\n - " +
        #               "\n - ".join(map(str, events)))
        if len(events) == 1:
            # just one event happening now
            item = events[0]
        elif all([_.tp == 'flip' for _ in events]):

            logging.debug("Only flip events, picking flip event with longest side, should guarantee progress")

            # # only flip events happening "now", pick event with longest side
            dist_events = []
            for evt in events:
                # compute distance of side to flip
                orig, dest = Edge(evt.triangle, evt.side[0]).segment
                dist = orig.distance2_at(dest, evt.time)
                dist_events.append((dist, evt))
            # take event with longest side
            dist_events.sort(reverse=True)
            item = dist_events[0][1]
            # item = events[0]
        else:
    ##        pts = []
    ##        for evt in events:
    ##            pts.extend([v.position_at(first.time)
    ##                        for v in evt.triangle.vertices])
    ##        on_straight_line = are_residuals_near_zero(pts)
    ##         print pts
    ##         for pt in pts:
    ##             print "POINT({0[0]} {0[1]})".format(pt)
    ##         metric = r_squared(pts)
    ##         logging.debug('metric {0}'.format(metric))
    ## #        import numpy as np
    ##         x = np.array([pt[0] for pt in pts])
    ##         y = np.array([pt[1] for pt in pts])
    ##         A = np.vstack([x, np.ones(len(x))]).T
    ##         res = np.linalg.lstsq(A, y)
    ##         print res

            # if not on_straight_line:
            #     logging.debug('pick first event')
            #     # points should be on straight line to exhibit infinite flipping
            #     # if this is not the case, we return the first event
            #     item = events[0]
            # else:
    #        logging.debug('pick non-flip event (vertices on straight line)')
            

            ##item = events[0]

            # # return the first item, sorted on type (non-flipping events first)
            for item in sorted(events, key=sort_key):
                break


    ## FIXME: 
    ## BUG: this code leads to problems in removing events from the queue
    ## my assumption is that:
    ## 1. we can navigate onto a triangle that already has been dealt with
    ## 2. the event on the triangle has not been updated properly (old event)
    def check_direct(event):
        """ check the direct neighbours whether they will collapse """
        
        current = event.triangle
        seen = set([current])
        visit = [current]
        also = []
        # FIXME: maybe we should use 'geometric inliers' to the support line here
        # this way we could also resolve 'flip event loops'
        while visit:
            current = visit.pop()
            for n in current.neighbours:
                if n is not None and n.event is not None and n not in seen:
                    seen.add(n)
                    if near_zero(n.event.time - event.time):
                        also.append(n)
                        visit.append(n)
                # elif n is not None and n.event is None:
                #     # not sure -- should we visit 'uncollapsable' triangles to check if their neighbours can collapse?
                #     seen.add(n)
                #     visit.append(n)
        events = [event] + [n.event for n in also]
        if events and all([_.tp == 'flip' for _ in events]):
            logging.debug("Only flip events, picking flip event with longest side, should guarantee progress")
            logging.debug(events)
            # # only flip events happening "now", pick event with longest side
            dist_events = []
            for evt in events:
                # compute distance of side to flip
                orig, dest = Edge(evt.triangle, evt.side[0]).segment
                dist = orig.distance2_at(dest, evt.time)
                dist_events.append((dist, evt))
            # take event with longest side
            dist_events.sort(reverse=True)
            return dist_events[0][1]
        else:
            return event

    evt = first
    queue.remove(evt)

    # evt = check_direct(evt)
    # try:
    #     queue.remove(evt)
    # except:
    #     print(evt)
    #     raise
    return evt





def log_queue_content(step, immediate, queue):
    logging.debug("")
    logging.debug("STEP := " + str(step))
    logging.debug("")
    logging.debug("=" * 80)
    for i, e in enumerate(immediate):
        logging.debug("{0:5d} {1}".format(i, e))
    logging.debug("-" * 80)
    for i, e in enumerate(queue):
        logging.debug("{0:5d} {1}".format(i, e))
        if i >= 20:
            break
    if len(queue) >= 20:
        logging.debug("... skipping display of {} events".format(len(queue) - 20))
    logging.debug("=" * 80)

# Main event loop
# -----------------------------------------------------------------------------
def event_loop(queue, skel, pause=False, stop_after=0, make_video=False, video_digits=3):
    """ The main event loop 
    
    Args:
        queue: Event queue
        skel: Skeleton structure
        pause: Whether to pause for interactive visualization
        stop_after: Stop after this many steps (0 = no limit) - for testing/debugging
        make_video: Whether to generate video frames - for testing/debugging  
        video_digits: Number of decimal digits for video timing - for testing/debugging
    """
    STOP_AFTER = stop_after
    if STOP_AFTER != 0:
        logging.debug("Stopping for the first time after step#{}".format(STOP_AFTER))
    VIDEO_DIGITS = video_digits
    # -- clean out files for visualization
    if pause:
        for file_nm in [
            "/tmpfast/sknodes_progress.wkt",
            "/tmpfast/bisectors_progress.wkt",
            "/tmpfast/segments_progress.wkt",
            '/tmpfast/queue.wkt',
            # also contains next triangle to be visualised!
            "/tmpfast/vertices0_progress.wkt",
            "/tmpfast/vertices1_progress.wkt",
            '/tmpfast/ktri_progress.wkt',
        ]:
            with open(file_nm, 'w') as fh:
                pass
    # -- visualize
    NOW = prev_time = 0. #= 5e-6
    step = 0
    if pause:
        logging.getLogger().setLevel(logging.DEBUG)
        # visualize(queue, skel, prev_time)
        # notify_qgis()
        interactive_visualize(queue, skel, step, NOW)
    immediate = deque([])
    logging.debug("=" * 80)
    logging.debug("Immediate / Queue at start of process")
    logging.debug("=" * 80)
    for i, e in enumerate(immediate):
        logging.debug("{0:5d} {1}".format(i, e))
    logging.debug("-" * 80)
    for i, e in enumerate(queue):
        logging.debug("{0:5d} {1}".format(i, e))
    logging.debug("=" * 80)
    
    if make_video:
        from grassfire.tests2.test_event_loop_debugging import make_frames
        make_frames(NOW, VIDEO_DIGITS, skel, queue, immediate)
#     step = prev = 0.025
#     FILTER_CT = 220
    check_bisectors(skel, 0.)
    check_active_triangles_orientation(skel.triangles, 0)
    guard = 0
    while queue or immediate:
        guard += 1
        if guard > 50000:
            raise ValueError('loop with more than 50_000 events stopped')
        #         if len(queue) < FILTER_CT:
        #             pause = True
        step += 1
        if pause:
            log_queue_content(step, immediate, queue)

        if immediate:
            evt = immediate.popleft()
        else:
            evt = choose_next_event(queue)
            NOW = evt.time

        # -- decide what to do based on event type
        logging.debug("About to handle event " +
                      str(evt.tp) +
                      " " +
                      str(evt.triangle.type) +
                      " " +
                      str(id(evt.triangle)) +
                      " [" +
                      str(evt.triangle.info) +
                      "] " +
                      " at time " +
                      "{0:.28g}".format(evt.time))

        if pause and step >= STOP_AFTER:

            def check_direct(event):
                """ check the direct neighbours whether they will collapse """
                current = event.triangle
                seen = set([current])
                visit = [current]
                also = []
                # FIXME: maybe we should use 'geometric inliers' to the support line here
                # this way we could also resolve 'flip event loops'
                while visit:
                    current = visit.pop()
                    for n in current.neighbours:
                        if n is not None and n.event is not None and n not in seen:
                            seen.add(n)
                            if near_zero(n.event.time - event.time):
                                also.append(n)
                                visit.append(n)
                        # elif n is not None and n.event is None:
                        #     # not sure -- should we visit 'uncollapsable' triangles to check if their neighbours can collapse?
                        #     seen.add(n)
                        #     visit.append(n)
                if also:
                    logging.debug([n.info for n in also])
            check_direct(evt)
            
            with open('/tmpfast/current_event.wkt', 'w') as fh:
                fh.write("pos;wkt;evttype;evttime;tritype;id;n0;n1;n2;finite;info;wavefront_directions\n")
                fh.write("{0};{1};{2};{3};{4};{5};{6};{7};{8};{9};{10}\n".format(
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
                    evt.triangle.wavefront_directions
                )
                )

            interactive_visualize(queue, skel, step, NOW)
        # precondition: this triangle has not yet been dealt with before
        if evt.triangle.stops_at is not None:
            logging.warn("Already stopped {}, but still queued".format(
                                                             id(evt.triangle)))
            # assert evt.triangle.stops_at is None, \
            #     "already stopped {}".format(id(evt.triangle))
            continue
        if evt.tp == "edge":
            if len(evt.side) == 3:
                # collapse to *single* point, as all 3 sides collapse
                handle_edge_event_3sides(evt, step, skel, queue, immediate)
            elif len(evt.side) == 1 and evt.triangle.type == 3:
                # collapse 3-triangle to line, as just 1 side collapses
                handle_edge_event_1side(evt, step, skel, queue, immediate, pause and step >= STOP_AFTER)
            elif len(evt.side) == 2:
                raise ValueError("Impossible configuration, triangle [{}] with 2 sides collapsing cannot happen".format(evt.triangle.info))
            else:
                handle_edge_event(evt, step, skel, queue, immediate, pause and step >= STOP_AFTER)
        elif evt.tp == "flip":
            handle_flip_event(evt, step, skel, queue, immediate)
        elif evt.tp == "split":
            handle_split_event(evt, step, skel, queue, immediate, pause and step >= STOP_AFTER)

        if False:
            log_queue_content(step, immediate, queue)

        if False and ((step % 100 == 0) or step >= STOP_AFTER):
            check_active_triangles_orientation(skel.triangles, NOW)

        if False:  # len(queue) < FILTER_CT:
            logging.debug("=" * 80)
            logging.debug("Immediate / Queue at end of handling event")
            logging.debug("=" * 80)
            for i, e in enumerate(immediate):
                logging.debug("{0:5d} {1}".format(i, e))
            logging.debug("-" * 80)
            for i, e in enumerate(queue):
                if i > 5 and i < len(queue) - 5:
                    continue
                logging.debug("{0:5d} {1}".format(i, e))
                logging.debug(repr(e.triangle))
                if i == 5 and len(queue) > 5:
                    logging.debug("...")
        logging.debug("=" * 80)
        if make_video:
            from grassfire.tests2.test_event_loop_debugging import make_frames
            make_frames(NOW, VIDEO_DIGITS, skel, queue, immediate)
        if False and not immediate:
            # if we have immediate events, the linked list will not be
            # ok for a while
            try:
                check_active_triangles(skel.triangles)
            except AssertionError:
                print ("{}".format(err))
                if True:
                    visualize(queue, skel, NOW - 5e-4)
                    # Interactive debugging removed - use pause mode instead
                    break
                else:
                    raise
    if pause:
        visualize(queue, skel, NOW)
    if make_video:
        from grassfire.tests2.test_event_loop_debugging import make_frames
        make_frames(NOW, VIDEO_DIGITS, skel, queue, immediate)

    if pause:
        interactive_visualize(queue, skel, step, NOW)

    # post condition - all finite triangles have stopped at the end of the wavefront propagation process
    not_stopped_tris = []
    for tri in skel.triangles:
        if all(v.internal for v in tri.vertices) and tri.stops_at is None:
            not_stopped_tris.append(tri.info)
    if not_stopped_tris:
        raise ValueError("triangles not stopped at end: {}".format(not_stopped_tris))

    return NOW


def compare_event_by_time(one, other):
    """ Compare two events, first by time, in case they are equal
    by triangle type (first 2-triangle, then 1-triangle, then 0-triangle),
    as last resort by identifier of triangle.
    """
    # compare by time
    if one.time < other.time:  # lt
        return -1
    elif one.time > other.time:  # gt
        return 1
    # in case times are equal, compare by id of triangle
    # to be able to find the correct triangle back
    else:  # eq
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
    """Compute for all kinetic triangles when they will collapse and
    put them in an OrderedSequence, so that events are ordered properly
    for further processing
    """
    q = OrderedSequence(cmp=compare_event_by_time)
    logging.debug("Calculate initial events")
    logging.debug("=" * 80)
    for tri in skel.triangles:
        res = compute_collapse_time(tri, 0, find_gt)
        if res is not None:
            q.add(res)
    logging.debug("=" * 80)
    return q
