import argparse
import cProfile
import pstats
import time
from statistics import mean

import requests
from triangle.delaunay.helpers import ToPointsAndSegments

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
    try:
        response = requests.get(BASE_URL + name, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as exc:
        raise RuntimeError(f"failed to load polygon archive input '{name}'") from exc


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
    profiler=None,
):
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    
    # Pre-load coordinates to exclude from timing and profiling
    loaded_coords = [load_coords_fn(name) for name in names]

    totals = []
    
    if profiler:
        profiler.enable()
        
    for _ in range(repeats):
        start = timer()
        for coords in loaded_coords:
            calc_segments_fn(coords)
        totals.append(timer() - start)
        
    if profiler:
        profiler.disable()
        
    return mean(totals), totals


def run_benchmark(
    repeats=3,
    profile=False,
    benchmark_fn=benchmark_total_skeleton_time,
    profiler_factory=cProfile.Profile,
):
    profiler = profiler_factory() if profile else None
    
    # Pass profiler into benchmark function so it can control scope
    average_total, totals = benchmark_fn(repeats=repeats, profiler=profiler)
    
    return average_total, totals, profiler


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
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Enable cProfile while running the benchmark.",
    )
    args = parser.parse_args()
    average_total, totals, profiler = run_benchmark(
        repeats=args.repeats, profile=args.profile
    )
    print(f"inputs={len(INPUT_NAMES)} repeats={args.repeats}")
    for i, total in enumerate(totals, start=1):
        print(f"run {i}: total_time={total:.6f}s")
    print(f"average_total_time={average_total:.6f}s")
    if profiler is not None:
        print("\nprofile_stats_cumulative")
        pstats.Stats(profiler).sort_stats("cumulative").print_stats(30)


if __name__ == "__main__":
    main()
