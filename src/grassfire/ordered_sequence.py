from functools import cmp_to_key


class OrderedSequence:
    """Ordered sequence where duplicates are allowed."""

    def __init__(self, cmp):
        self._cmp = cmp
        self._items = []

    def add(self, item):
        self._items.append(item)
        self._items.sort(key=cmp_to_key(self._cmp))

    def remove(self, item):
        self._items.remove(item)

    def __iter__(self):
        return iter(self._items)

    def __len__(self):
        return len(self._items)

    def __bool__(self):
        return bool(self._items)
