import pytest

from grassfire.benchmark_polygon_archive_segments import (
    benchmark_total_skeleton_time,
    run_benchmark,
)


def test_benchmark_total_skeleton_time_average():
    timer_calls = [0]
    times = iter((1.0, 3.0, 10.0, 16.0))
    def timer():
        timer_calls[0] += 1
        return next(times)
    avg, totals = benchmark_total_skeleton_time(
        names=("a", "b"),
        repeats=2,
        load_coords_fn=lambda name: [],
        calc_segments_fn=lambda coords: [],
        timer=timer,
    )
    assert totals == [2.0, 6.0]
    assert avg == 4.0
    assert timer_calls[0] == 4


def test_benchmark_total_skeleton_time_repeats_validation():
    with pytest.raises(ValueError, match="repeats must be >= 1"):
        benchmark_total_skeleton_time(
            names=(),
            repeats=0,
            load_coords_fn=lambda name: [],
            calc_segments_fn=lambda coords: [],
        )


def test_run_benchmark_without_profiling():
    avg, totals, profiler = run_benchmark(
        repeats=1,
        profile=False,
        benchmark_fn=lambda repeats: (1.25, [1.25]),
    )
    assert avg == 1.25
    assert totals == [1.25]
    assert profiler is None


def test_run_benchmark_with_profiling():
    class FakeProfiler:
        def __init__(self):
            self.enabled = False
            self.disabled = False

        def enable(self):
            self.enabled = True

        def disable(self):
            self.disabled = True

    avg, totals, profiler = run_benchmark(
        repeats=1,
        profile=True,
        benchmark_fn=lambda repeats: (2.5, [2.5]),
        profiler_factory=FakeProfiler,
    )
    assert avg == 2.5
    assert totals == [2.5]
    assert profiler.enabled
    assert profiler.disabled
