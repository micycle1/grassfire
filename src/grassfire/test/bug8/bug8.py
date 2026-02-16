"""
"""


if True:
    import logging
    import sys
    root = logging.getLogger()
    root.setLevel(logging.WARNING)

from tri.delaunay.helpers import ToPointsAndSegments
from grassfire import calc_skel

import timeit
with open('bug_exterior.wkt') as fh:
    fh.readline()  # skip header in file (line with "wkt")
    wkt = fh.readline()
    from simplegeom.wkt import loads
    poly = loads(wkt)
    poly = poly[0]


conv = ToPointsAndSegments()
for start, end in zip(poly[:-1], poly[1:]):
    conv.add_point(start)
    conv.add_point(end)
    conv.add_segment(start, end)
start = timeit.default_timer()
print(start)
skel = calc_skel(conv, pause=False, output=True, shrink=True, internal_only=True)
now = timeit.default_timer()
print now - start, "sec(s)"