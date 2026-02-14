import argparse
import time

import requests
from tri.delaunay.helpers import ToPointsAndSegments

from grassfire import calc_skel


INPUT_NAMES = (
    "eberly-10.json",
    "eberly-14.json",
    "elgindy-1.json",
    "gray-embroidery.json",
    "held-1.json",
    "held-12.json",
    "held-3.json",
    "held-7a.json",
    "held-7b.json",
    "held-7c.json",
    "held-7d.json",
    "mapbox-building.json",
    "mapbox-dude.json",
    "matisse-alga.json",
    "matisse-blue.json",
    "matisse-icarus.json",
    "matisse-nuit.json",
    "mei-2.json",
    "mei-3.json",
    "mei-4.json",
    "mei-5.json",
    "mei-6.json",
    "meisters-3.json",
    "misc-discobolus.json",
    "misc-fu.json",
    "seidel-3.json",
    "skimage-horse.json",
    "toussaint-1a.json",
)

BASE_URL = "https://raw.githubusercontent.com/LingDong-/interesting-polygon-archive/master/json/"


def load_coords(name):
    return requests.get(BASE_URL + name).json()


def calc_segments(coords):
    conv = ToPointsAndSegments()
    for ring in coords:
        for p in ring:
            conv.add_point(tuple(p))
        for i in range(len(ring)):
            start = tuple(ring[i])
            end = tuple(ring[(i + 1) % len(ring)])
            conv.add_segment(start, end)
    return calc_skel(conv, internal_only=True).segments()


def benchmark_total_skeleton_time(
    names=INPUT_NAMES,
    repeats=3,
    load_coords_fn=load_coords,
    calc_segments_fn=calc_segments,
    timer=time.perf_counter,
):
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    totals = []
    for _ in range(repeats):
        start = timer()
        for name in names:
            coords = load_coords_fn(name)
            calc_segments_fn(coords)
        totals.append(timer() - start)
    return sum(totals) / len(totals), totals


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark skeleton generation on polygon archive segment inputs."
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Number of full benchmark runs to average (default: 3).",
    )
    args = parser.parse_args()
    average_total, totals = benchmark_total_skeleton_time(repeats=args.repeats)
    print(f"inputs={len(INPUT_NAMES)} repeats={args.repeats}")
    for i, total in enumerate(totals, start=1):
        print(f"run {i}: total_time={total:.6f}s")
    print(f"average_total_time={average_total:.6f}s")


if __name__ == "__main__":
    main()
