import unittest
from grassfire.collapse import compute_event_3triangle, find_gt,\
    visualize_collapse
from grassfire.primitives import KineticTriangle, KineticVertex


class TestCollapseTime(unittest.TestCase):
    def setUp(self):
        self.now = 0
        self.triangle = KineticTriangle(KineticVertex((-5.586845556126167, 12.264500723992885), (42.96862677742887, -96.56675461024079), (-0.9097484670494593, -0.4151598808906743), (0.9174410515320596, 0.39787173431113276)),
            KineticVertex((0.07562289915966396, -0.792423319327731), (0.892367111503632, 0.45568901541445345), (0.9174410515320596, 0.39787173431113276), (0.8602330667847619, 0.5099010402127879)),
            KineticVertex((1.3619372870019553, -2.9625114016712693), (-8.665884533180343, 16.58102212313996), (0.8602330667847619, 0.5099010402127879), (-0.9097484670494593, -0.4151598808906743)), None, None, None)

        for t in [0., 0.13]:#, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.0947089087]:
             visualize_collapse(self.triangle,T=t)

    def test_bug(self):
        evt = compute_event_3triangle(self.triangle, now=self.now, sieve=find_gt)
        print(evt)





if __name__ == "__main__":
    if True:
        import logging
        import sys
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

    unittest.main(verbosity=2)