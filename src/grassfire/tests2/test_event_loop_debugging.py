"""Testing and debugging helpers for event loop visualization.

This module contains code extracted from the events module that was used
for testing, debugging, and visualization purposes. These utilities are
not part of the core library functionality.
"""

import random
import time
import logging

# Constants for video frame generation
DEFAULT_TIME_DELTA = 0.2  # Default time delta when queue is empty
SIGNAL_RANGE = int(1e6)   # Random signal range for frame synchronization


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
        # Use a small default time delta when queue is empty
        T = N + DEFAULT_TIME_DELTA
    
    delta = T - N
    if delta <= 0:
        return  # No frames to generate if next event is not in the future
    
    times = int(delta * scale)
    for t in range(1, times):
        print(".")
        cur = N + float(t) / float(scale)
        time.sleep(0.25)
        visualize(queue, skel, cur)
        time.sleep(0.5)
        # Write random signal for frame synchronization with external visualization
        with open("/tmpfast/signal", "w") as fh:
            fh.write("{0}".format(random.randint(0, SIGNAL_RANGE)))
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
