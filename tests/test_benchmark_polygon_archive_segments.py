import pytest

from grassfire.benchmark_polygon_archive_segments import benchmark_total_skeleton_time


def test_benchmark_total_skeleton_time_average():
    times = iter((1.0, 3.0, 10.0, 16.0))
    avg, totals = benchmark_total_skeleton_time(
        names=("a", "b"),
        repeats=2,
        load_coords_fn=lambda name: [],
        calc_segments_fn=lambda coords: [],
        timer=lambda: next(times),
    )
    assert totals == [2.0, 6.0]
    assert avg == 4.0


def test_benchmark_total_skeleton_time_repeats_validation():
    with pytest.raises(ValueError, match="repeats must be >= 1"):
        benchmark_total_skeleton_time(
            repeats=0,
            load_coords_fn=lambda name: [],
            calc_segments_fn=lambda coords: [],
        )
