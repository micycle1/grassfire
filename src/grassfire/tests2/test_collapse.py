"""Tests for collapse module"""

import logging
import sys

from grassfire.collapse import (
    compute_collapse_time,
    visualize_collapse,
    solve_quadratic
)
from grassfire.primitives import KineticTriangle, KineticVertex, InfiniteVertex
from grassfire.calc import near_zero


def test_compute_collapse_times():
    from grassfire.primitives import KineticTriangle, KineticVertex
    cases = [
        (0,
         KineticTriangle(InfiniteVertex((2., 4.)),
                         KineticVertex((2., 0.), (-0.5, -0.5)),
                         KineticVertex((1., 1.), (0.5, 0.)), True, True, True),
         (1.211102550928, "flip")
         ),
        (0,
         KineticTriangle(InfiniteVertex((1., 4.)),
                         KineticVertex((2., 0.), (-0.5, -0.5)),
                         KineticVertex((0., 0.), (0.5, -0.5)), True, True, True),
         (2.0, "edge")
         ),
        (0,
         KineticTriangle(InfiniteVertex((2., 4.)),
                         KineticVertex((4., 0.), (0.5, -0.5)),
                         KineticVertex((0., 0.), (-0.5, -0.5)), True, True, True),
            None),
        (0,
         KineticTriangle(InfiniteVertex((1., 4.)),
                         KineticVertex((2., 0.), (-0.5, -0.5)),
                         KineticVertex((0., 0.), (0.5, -0.5)), None, True, True),
            (2.0, "edge")
         ),
        (0,
            KineticTriangle(KineticVertex((2., 4.), (0., -0.5)),
                            KineticVertex((0., 0.), (0., +1)),
                            KineticVertex((2., 0.), (0., +1)), None, True, True),
            (2.666666666667, "edge")),
        (0,
         KineticTriangle(KineticVertex((1., 4.), (0., -0.5)),
                         KineticVertex((0., 0.), (0., 1.)),
                         KineticVertex((2., 0.), (0., 1.)), None, True, True),
         (2.666666666667, "split")),
        (0,
         KineticTriangle(KineticVertex((1., 4.), (0., -0.5)),
                         KineticVertex((0., 0.), (0., +1.0)),
                         KineticVertex((2., 0.), (0., +1.0)),
                         None, True, True),
            (2.666666666666666666, "split")),
        (0,
            KineticTriangle(KineticVertex((1., 4.), (0., -0.5)),
                            KineticVertex((0., 0.), (0.5, .5)),
                            KineticVertex((2., 0.), (-.5, .5)),
                            None, True, True),
            (2., "edge")),
        (0,
         KineticTriangle(KineticVertex((3., 4.), (0., -0.5)),
                         KineticVertex((0., 0.), (0., +1.0)),
                         KineticVertex((2., 0.), (0., +1.0)),
                         None, True, True),
         (2.666666666667, "flip")),
        (0,
         KineticTriangle(
             KineticVertex(
                 (11.1,
                  0.0),
                 (-0.41421356237309614,
                  1.0)),
             KineticVertex(
                 (14.0,
                  10.0),
                 (-0.41434397951188867,
                  -1.0828510849683515)),
             KineticVertex(
                 (-33.307692307692705,
                  2.115384615384554),
                 (9.300198345114286,
                  0.536239302469343)),
             None,
             True,
             True),
            (4.569105324086005, "split")),
        (0,
         KineticTriangle(KineticVertex((0., 0.), (0.5, 0.5)),
                         KineticVertex((2., 0.), (-0.5, 0.5)),
                         KineticVertex((1., 1.), (0., -0.5)),
                         None, None, None),
            (1.2, "edge")),
        (0,
         KineticTriangle(
             KineticVertex(
                 (11.1, 0.0), (-0.41421356237309614, 1.0)), KineticVertex(
                 (9.0, 0.0), (0.024984394500786274, 1.0)), KineticVertex(
                 (11.0, -0.1), (-0.39329937395053033, 1.0209141884225659)), None, None, True),
            (4.781443007949, "edge")),

        (0,
         KineticTriangle(
             KineticVertex(
                 (-0.866025403784, 0.5), (0.3660254037831927, -1.3660254037847726)), KineticVertex(
                 (-0.866025403784, 0.5), (-1.3660254037847726, -0.3660254037831926)), InfiniteVertex(
                 (2.4999999999993334, 1.4433756729733331)), True, True, None),
            (0, "edge")
         ),

        (0,
         KineticTriangle(
             KineticVertex(
                 (3.36602540378, 4.83012701892), (-1.366025403784519, 0.366025403784139)), KineticVertex(
                 (1.63397459622, 4.83012701892), (1.3660254037847726, 0.366025403783193)), KineticVertex(
                 (5.86602540378, 0.5), (-0.366025403784139, -1.3660254037845188)), True, None, True),
            (0.6339745962123428, "edge")),

        (0,
         KineticTriangle(KineticVertex((-
                                        0.9510565162951535, 0.3090169943749475), (1.6827457682352098, -
                                                                                  0.5467572438521933)), KineticVertex((6.123233995736766e-17, 1.0), (-
                                                                                                                                                     2.3811458388420067e-16, -
                                                                                                                                                     1.7693436082961256)), KineticVertex((-
                                                                                                                                                                                          1.6180339887498947, 1.1755705045849465), (0.18250881904109725, 1.4023874396799996)), True, None, True),
         (0.56518134482820892, "edge")),

        (0,
         KineticTriangle(
             KineticVertex(
                 (-0.2, -0.06666666666666667), (1.8695058979924137, 0.04448684017361237)), KineticVertex(
                 (0.3333333333333333, -0.6), (-0.4142135623730951, 1.0)), KineticVertex(
                 (0.3333333333333333, 0.6), (-0.2506607793572, -1.0800962240008098)), True, True, True),
            (0.241024119, "flip")),

        (0,
         KineticTriangle(KineticVertex((-
                                        0.9872805569585875, 0.12537065799230004), (-
                                                                                   1.024866765183742, -
                                                                                   1.1114119702481062)), KineticVertex((-
                                                                                                                        0.8437938978984622, -
                                                                                                                        0.19040991430610052), (21.666935207544306, 44.860914394127434)), KineticVertex((-
                                                                                                                                                                                                        0.7653348342872799, -
                                                                                                                                                                                                        0.10537661431275783), (-
                                                                                                                                                                                                                               1.0428319428196307, -
                                                                                                                                                                                                                               0.08857245780155545)), True, True, True),
            (0.00242630813253, "flip")),

        (0,
         KineticTriangle(KineticVertex((-
                                        0.25, 0.75), (-
                                                      2.4142135623730945, -
                                                      0.9999999999999996)), KineticVertex((-
                                                                                           0.25, -
                                                                                           0.75), (2.4142135623730945, 0.9999999999999996)), KineticVertex((0.25, -
                                                                                                                                                            0.75), (-
                                                                                                                                                                    2.4142135623730945, 0.9999999999999996)), True, True, True),
            (0.103553390593274, "edge")
         ),
        (0.01,
         KineticTriangle(KineticVertex((-
                                        0.2514204545452013, -
                                        0.43678977272891734), (-
                                                               1.1213425200822953, 0.8653503174771875)), KineticVertex((-
                                                                                                                        0.39342794594029157, 0.34274872190648226), (1.1112498367417047, -
                                                                                                                                                                    0.7879730785326)), KineticVertex((-
                                                                                                                                                                                                      0.39346590909065293, 0.3430397727259427), (1.119825823231893, -
                                                                                                                                                                                                                                                 0.8537223082927444)), None, True, True),
            (0.019994907, "flip")
         ),
        (0,
         KineticTriangle(KineticVertex((0.4348145985280354, -
                                        0.45254225871408715), (-
                                                               0.9629719544552265, 0.27936526970827974)), KineticVertex((0.45503984798207725, -
                                                                                                                         0.35731302646132906), (-
                                                                                                                                                0.8099200307905356, 0.9999999999999997)), KineticVertex((0.31446758274483105, -
                                                                                                                                                                                                         0.6518477551830741), (-
                                                                                                                                                                                                                               1.413908790141726, -
                                                                                                                                                                                                                               0.02935869819935455)), True, True, None),
            None)
    ]
    do_test = True
    if do_test:
        for i, (now, tri, expected) in enumerate(cases, start=0):
            print()
            print(( "Case", i))
            print(( "=" * 10))
            visualize_collapse(tri, now)
            evt = compute_collapse_time(tri, now)
            print(( evt, expected))
            if evt is not None:
                if expected is None:
                    assert False, "incorrect None"

                else:
                    time, tp = expected
                    assert near_zero(evt.time - time)
                    assert evt.tp == tp
            else:
                assert evt is expected

    now = 0.
    tri = cases[-1][1]

    try:
        evt = compute_collapse_time(tri, now)
        print (evt)
    except:
        pass
    visualize_collapse(tri, now)


def test_solve():
    A, B, C = 3.0, 4.5, 9.0
    print((solve_quadratic(A, B, C) == []))

    A, B, C = 2.0, 0.0, 0.0
    print((solve_quadratic(A, B, C) == [0.0]))

    A, B, C = 1.0, 3.0, 2.0
    print((solve_quadratic(A, B, C) == [-2.0, -1.0]))


def main():
    test_one_collapse()


def test_one_collapse():


    from grassfire.primitives import KineticTriangle, KineticVertex

    tri = KineticTriangle(KineticVertex((-0.25890832526488483,
    1.0),
    (0.8629392738894351,
    -1.2541589253803134),
    (-0.24881895085820574,
    -0.9685500140384191),
    (0.993509292081015,
    -0.11375098482510204)),

    KineticVertex((-0.4180575328215234,
    -0.3900206383088658),
    (-0.4172703730965519,
    -12.435602163414288),
    (0.993509292081015,
    -0.11375098482510204),
    (-0.9988997443815248,
    -0.04689670217109577)),

    KineticVertex((-0.4159764930624239,
    -0.4312893163239368),
    (-1.0459590484062025,
    0.9554664616505931),
    (-0.9988997443815248,
    -0.04689670217109577),
    (-0.04343900602459881,
    0.9990560808861507)),
    None,
    True,
    None)

    for info, v in enumerate(tri.vertices, start = 1):
        v.info = info

    now = 0.003107819648575184823952044511
    print(now)
    evt = compute_collapse_time(tri, now)
    print(evt)
    times = [now - 0.1, now+0, now+0.05, now+0.1]
    for time in sorted(times):
        visualize_collapse(tri, time)
        input("paused at " + str(time))


if __name__ == "__main__":

    import logging
    import sys
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    ch.setFormatter(formatter)
    root.addHandler(ch)
    main()
