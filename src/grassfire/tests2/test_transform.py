"""Tests for transform module"""

from grassfire.transform import get_transform, get_box


def _test_transform():
    box = ((70000., 100000.), (75000., 125000.))
    t = get_transform(box)

    assert t.forward(box[0]) == (-0.2, -1)
    assert t.forward((72500., 112500.)) == (0, 0)
    assert t.forward(box[1]) == (0.2, 1)

    assert t.backward((-0.2, -1)) == box[0]
    assert t.backward((0.2, 1)) == box[1]
    assert t.backward((0, 0)) == (72500, 112500)


def _test_box():
    pts = [(78000., 100000.), (75000., 125000.)]
    assert get_box(pts) == ((75000.0, 100000.0), (78000.0, 125000.0))


if __name__ == "__main__":
    _test_transform()
    _test_box()
