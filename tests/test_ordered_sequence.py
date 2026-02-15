import importlib.util
from pathlib import Path


_MODULE_PATH = Path(__file__).resolve().parents[1] / "src" / "grassfire" / "ordered_sequence.py"
_SPEC = importlib.util.spec_from_file_location("grassfire.ordered_sequence", _MODULE_PATH)
_MODULE = importlib.util.module_from_spec(_SPEC)
assert _SPEC is not None and _SPEC.loader is not None
_SPEC.loader.exec_module(_MODULE)
OrderedSequence = _MODULE.OrderedSequence


def _cmp_int(one, other):
    if one < other:
        return -1
    if one > other:
        return 1
    return 0


def test_ordered_sequence_sorts_and_allows_duplicates():
    queue = OrderedSequence(cmp=_cmp_int)
    queue.add(3)
    queue.add(1)
    queue.add(2)
    queue.add(2)

    assert list(queue) == [1, 2, 2, 3]


def test_ordered_sequence_remove_removes_single_item():
    queue = OrderedSequence(cmp=_cmp_int)
    queue.add(2)
    queue.add(2)
    queue.remove(2)

    assert list(queue) == [2]


def test_ordered_sequence_discard_ignores_missing_item():
    queue = OrderedSequence(cmp=_cmp_int)
    queue.add(1)
    queue.discard(2)

    assert list(queue) == [1]
