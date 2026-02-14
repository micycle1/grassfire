"""Testing and debugging helpers for event loop visualization.

This module contains code extracted from the events module that was used
for testing, debugging, and visualization purposes. These utilities are
not part of the core library functionality.
"""

import random
import time


def make_frames(now, digits, skel, queue, immediate):
    """Generate frames for video visualization during event loop execution.
    
    This is a testing/debugging utility that was embedded in the event loop.
    It creates visualization frames by writing to /tmpfast/ directory.
    
    Args:
        now: Current time in the event loop
        digits: Number of digits for rounding time values
        skel: Skeleton structure being computed
        queue: Event queue
        immediate: Immediate events queue
        
    Note:
        This function assumes T >= N (next event time >= current time).
        If delta is negative, no frames will be generated.
    """
    from grassfire.inout import visualize
    
    if immediate:
        return
    
    scale = pow(10, digits)
    N = round(now, digits)
    try:
        peek = next(iter(queue))
        T = round(peek.time, digits)
    except StopIteration:
        T = N + 0.2
    
    delta = T - N
    if delta <= 0:
        return  # No frames to generate if next event is not in the future
    
    times = int(delta * scale)
    for t in range(1, times):
        print( ".")
        cur = N + float(t) / float(scale)
        time.sleep(0.25)
        visualize(queue, skel, cur)
        time.sleep(0.5)
        with open("/tmpfast/signal", "w") as fh:
            fh.write("{0}".format(random.randint(0, int(1e6))))
        time.sleep(0.25)


# Testing constants that were hardcoded in loop.py
TEST_STOP_AFTER_VALUES = [14830, 14569, 14851, 104700]
"""List of STOP_AFTER values used for testing specific scenarios."""

TEST_VIDEO_DIGITS = 3
"""Number of decimal digits for video frame timing."""

TEST_MAKE_VIDEO = False
"""Flag to enable video generation during testing."""


if __name__ == "__main__":
    print("Testing and debugging utilities for grassfire event loop")
    print(f"Available STOP_AFTER test values: {TEST_STOP_AFTER_VALUES}")
    print(f"Video generation disabled by default: make_video={TEST_MAKE_VIDEO}")
