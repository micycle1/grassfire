import logging

import pytest
import requests
from tri.delaunay.helpers import ToPointsAndSegments

from grassfire import calc_skel
from grassfire.test.intersection import segments_intersecting


LOG = logging.getLogger(__name__)

JSON_EXPECTED_SEGMENTS = {
    "eberly-10.json": 36,
    "eberly-14.json": 39,
    "elgindy-1.json": 31,
    "gray-embroidery.json": 61,
    "held-1.json": 122,
    "held-12.json": 63,
    "held-3.json": 118,
    "held-7a.json": 125,
    "held-7b.json": 123,
    "held-7c.json": 125,
    "held-7d.json": 119,
    "mapbox-building.json": 26,
    "mapbox-dude.json": 211,
    "matisse-alga.json": 494,
    "matisse-blue.json": 394,
    "matisse-icarus.json": 178,
    "matisse-nuit.json": 245,
    "mei-2.json": 49,
    "mei-3.json": 81,
    "mei-4.json": 231,
    "mei-5.json": 555,
    "mei-6.json": 732,
    "meisters-3.json": 57,
    "misc-discobolus.json": 291,
    "misc-fu.json": 352,
    "seidel-3.json": 37,
    "skimage-horse.json": 235,
    "toussaint-1a.json": 247,
}

BASE_URL = "https://raw.githubusercontent.com/LingDong-/interesting-polygon-archive/master/json/"


def _load_coords(name: str):
    url = BASE_URL + name
    return requests.get(url).json()


def _calc_segments(coords):
    conv = ToPointsAndSegments()
    for ring in coords:
        for p in ring:
            conv.add_point(tuple(p))
        for i in range(len(ring)):
            start = tuple(ring[i])
            end = tuple(ring[(i + 1) % len(ring)])
            conv.add_segment(start, end)

    skel = calc_skel(conv, internal_only=True)
    return skel.segments()


def _segment_pairs(segments):
    return [segment for segment, _infos in segments]


class _OnlyThisLoggerFilter(logging.Filter):
    def filter(self, record):
        return record.name == __name__


@pytest.fixture(autouse=True, scope="module")
def _filter_logging():
    root = logging.getLogger()
    filt = _OnlyThisLoggerFilter()
    root.addFilter(filt)
    yield
    root.removeFilter(filt)


@pytest.mark.parametrize("name,expected", JSON_EXPECTED_SEGMENTS.items())
def test_segment_counts(name, expected):
    coords = _load_coords(name)
    segments = _calc_segments(coords)
    LOG.info("%s expected=%s actual=%s", name, expected, len(segments))
    assert len(segments) == expected


@pytest.mark.parametrize("name", JSON_EXPECTED_SEGMENTS.keys())
def test_segments_do_not_cross(name):
    coords = _load_coords(name)
    segments = _calc_segments(coords)
    assert not segments_intersecting(_segment_pairs(segments)), (
        "intersection between straight skeleton segments found"
    )
