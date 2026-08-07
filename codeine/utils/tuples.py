from typing import Any, Tuple, Type, Union


def tuplify(value: Any, item_type: Union[Type, Tuple[Type, ...]]) -> Tuple:
    """
    Convert a value to a tuple, treating values of the specified type as single items.

    If ``value`` is an instance of ``item_type``, it is wrapped in a one-item tuple.
    Otherwise, ``value`` is converted to a tuple. If ``value`` is None, return an empty tuple.

    Parameters
    ----------
    value
        The value to convert.
    item_type
        The type or types that should be treated as a single item.
    """
    if value is None:
        return ()

    if isinstance(value, item_type):
        return (value,)

    return tuple(value)
