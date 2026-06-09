from collections.abc import Mapping
from typing import Any, Iterator


class FrozenDict(Mapping):
    """
    Immutable (at the shallow level) version of a dict.
    """

    def __init__(self, data: Mapping) -> None:
        self._data = dict(data)

    def __getitem__(self, key: Any) -> Any:
        return self._data[key]

    def __iter__(self) -> Iterator:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __repr__(self) -> str:
        return repr(self._data)
