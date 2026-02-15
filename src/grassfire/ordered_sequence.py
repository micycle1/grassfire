from bisect import bisect_right
from functools import cmp_to_key


class OrderedSequence:
    """Ordered sequence where duplicates are allowed."""

    def __init__(self, cmp):
        self._cmp = cmp
        self._key = cmp_to_key(cmp)
        self._items = []
        self._keys = []

    def add(self, item):
        key_item = self._key(item)
        index = bisect_right(self._keys, key_item)
        self._items.insert(index, item)
        self._keys.insert(index, key_item)

    def remove(self, item):
        for index, existing in enumerate(self._items):
            if existing == item:
                del self._items[index]
                del self._keys[index]
                return
        raise ValueError("item not found in ordered sequence")

    def discard(self, item):
        try:
            self.remove(item)
        except ValueError:
            pass

    def __iter__(self):
        return iter(self._items)

    def __len__(self):
        return len(self._items)

    def __bool__(self):
        return bool(self._items)
