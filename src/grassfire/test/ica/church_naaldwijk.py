from tri.delaunay.helpers import ToPointsAndSegments
from grassfire import calc_skel
from simplegeom.wkt import loads

"""Church in Naaldwijk
"""

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

with open("/home/martijn/Documents/work/archive/2016-01_grassfire_for_building_generalization/data/naaldwijk_church/in_out/in.geojson") as fh:
    s = fh.read()

import json
x = json.loads(s)
segments = []
for y in x['features']:
    segments.append(tuple(map(tuple, y['geometry']['coordinates'])))
conv = ToPointsAndSegments()
for line in segments:
    conv.add_point(line[0])
    conv.add_point(line[1])
    conv.add_segment(*line)
skel = calc_skel(conv, pause=False, output=True, shrink=False, internal_only=False)
