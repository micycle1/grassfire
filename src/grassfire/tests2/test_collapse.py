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
        # infinite 0-triangle
        (0,
         KineticTriangle(InfiniteVertex((2., 4.)),
                         KineticVertex((2., 0.), (-0.5, -0.5)),
                         KineticVertex((1., 1.), (0.5, 0.)), True, True, True),
         (1.211102550928, "flip")
         ),
        # infinite 0-triangle
        (0,
         KineticTriangle(InfiniteVertex((1., 4.)),
                         KineticVertex((2., 0.), (-0.5, -0.5)),
                         KineticVertex((0., 0.), (0.5, -0.5)), True, True, True),
         (2.0, "edge")
            # None
         ),
        # infinite 0-triangle
        (0,
         KineticTriangle(InfiniteVertex((2., 4.)),
                         KineticVertex((4., 0.), (0.5, -0.5)),
                         KineticVertex((0., 0.), (-0.5, -0.5)), True, True, True),
            None),
        # infinite 1-triangle
        (0,
         KineticTriangle(InfiniteVertex((1., 4.)),
                         KineticVertex((2., 0.), (-0.5, -0.5)),
                         KineticVertex((0., 0.), (0.5, -0.5)), None, True, True),
            (2.0, "edge")
         ),
        # finite 1-triangle
        # miss edge
        (0,
            KineticTriangle(KineticVertex((2., 4.), (0., -0.5)),
                            KineticVertex((0., 0.), (0., +1)),
                            KineticVertex((2., 0.), (0., +1)), None, True, True),
            (2.666666666667, "edge")),
        # finite 1-triangle
        (0,
         KineticTriangle(KineticVertex((1., 4.), (0., -0.5)),
                         KineticVertex((0., 0.), (0., 1.)),
                         KineticVertex((2., 0.), (0., 1.)), None, True, True),
         (2.666666666667, "split")),
        # finite 1-triangle
        # rot-0
        (0,
         KineticTriangle(KineticVertex((1., 4.), (0., -0.5)),
                         KineticVertex((0., 0.), (0., +1.0)),
                         KineticVertex((2., 0.), (0., +1.0)),
                         None, True, True),
            (2.666666666666666666, "split")),
        # finite 1-triangle -- wavefront edge collapses
        (0,
            KineticTriangle(KineticVertex((1., 4.), (0., -0.5)),
                            KineticVertex((0., 0.), (0.5, .5)),
                            KineticVertex((2., 0.), (-.5, .5)),
                            None, True, True),
            (2., "edge")),
        # finite 1-triangle -- apex misses the wavefront
        (0,
         KineticTriangle(KineticVertex((3., 4.), (0., -0.5)),
                         KineticVertex((0., 0.), (0., +1.0)),
                         KineticVertex((2., 0.), (0., +1.0)),
                         None, True, True),
         (2.666666666667, "flip")),
        # finite 1 triangle that should split
        # we miss the event
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
        # 3-triangle
        (0,
         KineticTriangle(KineticVertex((0., 0.), (0.5, 0.5)),
                         KineticVertex((2., 0.), (-0.5, 0.5)),
                         KineticVertex((1., 1.), (0., -0.5)),
                         None, None, None),
            (1.2, "edge")),
        # 2-triangle collapse to point
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
            # None
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
#     now = 4.781443007949
    tri = cases[-1][1]

    try:
        evt = compute_collapse_time(tri, now)
        print (evt)
    except:
        pass
    visualize_collapse(tri, now)

#     print solve_quadratic(*area_collapse_time_coeff(*tri.vertices))
#     print solve_quadratic_old(*area_collapse_time_coeff(*tri.vertices))
#     import matplotlib.pyplot as plt
#     areas = []
#     times = range(-30, 20)
#     for t in times:
#         area = orient2d(tri.vertices[0].position_at(t),
#                         tri.vertices[1].position_at(t),
#                         tri.vertices[2].position_at(t)
#                     )
#         areas.append(area)
#
#     def distance(p0, p1):
#         return sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)
#     distances = []
#     for t in [t/ 10. for t in times]:
#         distances.append(distance(tri.vertices[0].origin, tri.vertices[1].position_at(t)))
#
#     plt.plot([t/ 10. for t in times], distances)
#     plt.grid(True)
#     plt.show()
#
#
#     if evt != None:
#         t = evt.time
#     else:
#         t = now
#     visualize_collapse(tri, 1.2111)


def test_solve():
    A, B, C = 3.0, 4.5, 9.0
    print((solve_quadratic(A, B, C) == []))

    A, B, C = 2.0, 0.0, 0.0
    print((solve_quadratic(A, B, C) == [0.0]))

    A, B, C = 1.0, 3.0, 2.0
    print((solve_quadratic(A, B, C) == [-2.0, -1.0]))


def main():
    # test_compute_collapse_times()
    test_one_collapse()


def test_one_collapse():


    from grassfire.primitives import KineticTriangle, KineticVertex
    # the vertex below should split the triangle in pieces!!!
    # -> we miss this event with our current calculation!!!

#     tri = KineticTriangle(KineticVertex((3.36602540378, 4.83012701892), (-1.366025403784519, 0.366025403784139)), KineticVertex((1.63397459622, 4.83012701892), (1.3660254037847726, 0.366025403783193)), KineticVertex((5.86602540378, 0.5), (-0.366025403784139, -1.3660254037845188)), True, None, True)
#     tri = KineticTriangle(KineticVertex((11.1, 0.0), (-0.41421356237309614, 1.0)), KineticVertex((14.0, 10.0), (-0.41434397951188867, -1.0828510849683515)), KineticVertex((-33.307692307692705, 2.115384615384554), (9.300198345114286, 0.536239302469343)), None, True, True)
    #tri = KineticTriangle(KineticVertex((-0.9510565162951535, 0.3090169943749475), (1.6827457682352098, -0.5467572438521933)), KineticVertex((6.123233995736766e-17, 1.0), (-2.3811458388420067e-16, -1.7693436082961256)), KineticVertex((-1.6180339887498947, 1.1755705045849465), (0.18250881904109725, 1.4023874396799996)), True, None, True)
    # tri = KineticTriangle(KineticVertex((-0.587785252292473, 0.8090169943749475), (1.0399940791944127, -1.4314290480002558)), KineticVertex((-0.9510565162951535, 0.3090169943749475), (1.6827457682352098, -0.5467572438521933)), KineticVertex((-0.5877852522924732, -0.8090169943749473), (1.0399940791944131, 1.4314290480002556)), True, True, True)

#     tri = KineticTriangle(KineticVertex([3.36092084343274, 4.81953957828363], (0, 0)), KineticVertex((0.0, 0.0), (1.0, 1.0)), KineticVertex((14.0, 10.0), (-0.41434397951188867, -1.0828510849683515)), True, None, None)
    # tri = KineticTriangle(KineticVertex((2.0, 0.5), (2.4142135623730945, -0.9999999999999996)), KineticVertex((1.0, 0.0), (2.4142135623730945, 0.9999999999999996)), KineticVertex((2.0, 0.0), (-2.4142135623730945, 0.9999999999999996)), True, True, True)
    #tri = KineticTriangle(KineticVertex((2.0, 0.5), (2.4142135623730945, -0.9999999999999996)), KineticVertex((1.0, 0.0), (2.4142135623730945, 0.9999999999999996)), KineticVertex((2.0, 0.0), (-2.4142135623730945, 0.9999999999999996)), True, True, True)

    #tri = KineticTriangle(KineticVertex((-0.2, -0.06666666666666667), (1.8695058979924137, 0.04448684017361237)), KineticVertex((0.3333333333333333, -0.6), (-0.4142135623730951, 1.0)), KineticVertex((0.3333333333333333, 0.6), (-0.2506607793572, -1.0800962240008098)), True, True, True)
    #tri = KineticTriangle(KineticVertex((-0.6388129767693403, -0.31936665619917365), (0.8030798568408047, -0.5996411436857993)), KineticVertex((-0.6231441738161808, -0.46672335701255496), (0.7799096015853287, 0.6333147760081715)), KineticVertex((-0.4881354888333025, -0.4609454744153349), (-0.3283637586937928, 0.94926819247684)), True, True, True)

    #tri = KineticTriangle(KineticVertex((-0.9872805569585875, 0.12537065799230004), (-1.024866765183742, -1.1114119702481062)), KineticVertex((-0.8437938978984622, -0.19040991430610052), (21.666935207544306, 44.860914394127434)), KineticVertex((-0.7653348342872799, -0.10537661431275783), (-1.0428319428196307, -0.08857245780155545)), True, True, True)

    #tri = KineticTriangle(KineticVertex((141.02031145318313, 15.11765139459854), (-1.09679614706904, 0.8919821935206798)), InfiniteVertex((0.17181892901794138, -0.016845276247151694)), KineticVertex((-0.11765994741487949, 0.576486269356072), (0.892675917072541, 1.09692483546051)), True, None, True)
    #tri = KineticTriangle(KineticVertex((141.02031145318313, 15.11765139459854), (-1.09679614706904, 0.8919821935206798)), InfiniteVertex((0.17181892901794127, -0.016845276247151646)), KineticVertex((-0.11765994741487949, 0.576486269356072), (0.892675917072541, 1.09692483546051)), True, None, True)
    #tri = KineticTriangle(KineticVertex((0.37664329535438995, -0.20429813029318875), (-1.0914440709287114, 0.8995592632828218)), KineticVertex((0.37941863862056496, 0.2502008472118177), (-0.8920557213641949, -1.0976424754293432)), KineticVertex((0.2307551855093853, -0.21963555360760362), (-1.0915314233653777, 0.8995177965775215)), True, True, True)

    #tri = KineticTriangle(KineticVertex((-0.25, 0.75), (-2.4142135623730945, -0.9999999999999996)), KineticVertex((-0.25, -0.75), (2.4142135623730945, 0.9999999999999996)), KineticVertex((0.25, -0.75), (-2.4142135623730945, 0.9999999999999996)), True, True, True)

    # infinite triangle that needs to flip (so around 0.14?
    #tri = KineticTriangle(InfiniteVertex((0.03563328767477883, -0.06225382499568152)), KineticVertex((0.8308871493803005, -0.35731302646132906), (1.248735516282837, 1.0000000000000002)), KineticVertex((0.7971434930656011, -0.5079794287841405), (0.9543989274069261, -0.3142213903649677)), None, True, True)
    # tri = KineticTriangle(KineticVertex((0.8534771264074428, -0.06782665418758393), (1.0, -1.0)), InfiniteVertex((0.03563328767477881, -0.06225382499568156)), KineticVertex((0.8534771264074428, 0.10897926770469894), (1.0, 0.03649041282240072)), True, None, True)

    #tri = KineticTriangle(KineticVertex((0.8308871493803005, -0.35731302646132906), (1.248735516282837, 1.0000000000000002)), InfiniteVertex((0.0356332876747788, -0.06225382499568155)), KineticVertex((0.8534771264074428, -0.06782665418758393), (1.0, -1.0)), True, True, True)

    # infinite 0-triangle should flip at: 0.14372166332027514
    # tri = KineticTriangle(KineticVertex((0.8308871493803005, -
    # 0.35731302646132906), (1.248735516282837, 1.0000000000000002)), InfiniteVertex((0.035633287674778816, -
#                                                                                                                          0.06225382499568152)), KineticVertex((0.8534771264074428, -
#                                                                                                                                                                0.06782665418758393), (1.0, -
#                                                                                                                                                                                       1.0)), True, True, True)
#
#     tri = KineticTriangle(KineticVertex((0.2, 0.1), (-0.5, -1.0), (-0.0, -1.0), (-0.8, -0.6000000000000001)),
# KineticVertex((0.2, -0.1), (-0.5, 1.0), (-0.8, 0.6000000000000001), (-0.0, 1.0)),
# KineticVertex((1.0, 0.1), (-0.9999999999999999, -0.9999999999999999), (-1.0, 0.0), (-0.0, -1.0)), True, None, True)




    # tri = KineticTriangle(KineticVertex((0.2, 0.1), (-0.5, -1.0), (-0.0, -1.0), (-0.8, -0.6000000000000001)),
    #                       KineticVertex((0.2, -0.1), (-0.5, 1.0),
    #                                     (-0.8, 0.6000000000000001), (-0.0, 1.0)),
    #                       KineticVertex((1.0, 0.1), (-0.9999999999999999, -0.9999999999999999), (-1.0, 0.0), (-0.0, -1.0)), True, None, True)
    # now = 0.  # 0.6339745962123428

#    def KineticVertex(origin=None, velocity=None, ul=None, ur=None)

##    tri = KineticTriangle(KineticVertex((-0.2941460749857284, -1.9944028440429675),
##        #(-31.916444710007777, 625.4387471202283),
##        (0, 0), # static
##        (-0.9987805884988515, -0.049369383608547396),
##        (0.998617829332686, 0.05255883311941421)),

##        KineticVertex((-0.3764671693980948, -0.4303110631701819),
##        (0.9986179074408387,
##        0.052557349064511094),
##        (0.998617829332686,
##        0.05255883311941421),
##        (0.9986179855445806,
##        0.05255586500937582)),

##        KineticVertex((-0.37147762485121205, -0.4299261332026524),
##        (-0.9986662884239297,
##        -0.05168176200486036),
##        (-0.9985412823874946,
##        -0.05399358635928628),
##        (-0.9987805884988515,
##        -0.049369383608547396)),
##        True,
##        None,
##        None)

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


#    tri = KineticTriangle(KineticVertex((0.2173277124280331, 0.698416760614311), (0.41124654232084884, 0.919545641592878), (0.29546368892668795, 0.9553539702779454), (0.5151235666701907, 0.8571159262672592)),
#KineticVertex((0.21178657399529294, 0.9831390013046308), (-0.18724127375184607, -0.9872847139880617), (-0.28216706807814973, -0.9593652827219572), (-0.08868469195702208, -0.9960597499209011)),
#KineticVertex((-1.07262137319989, 1.0973611594703172), (9.311325418444167, -1.8329944831033405), (-0.08868469195702208, -0.9960597499209011), (0.29546368892668795, 0.9553539702779454)), None, None, True)
#    now = 0.1304109494714286443528550308

    for info, v in enumerate(tri.vertices, start = 1):
        v.info = info

    now = 0.003107819648575184823952044511
    print(now)
    evt = compute_collapse_time(tri, now)
    print(evt)
    # 0.00242630813252781, 0.6340506109731798, 0.004284474881621788,
    # 0.0022096098886525, 0.22933526207436553]
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
