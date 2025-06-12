def bounded(value: int | float, min_value: int | float, max_value: int | float) -> int | float:
    """
    Bound a value to be within a specified range.

    Parameters
    ----------
    value : int | float
        The value to be bounded.
    min_value : int | float
        The minimum value of the range.
    max_value : int | float
        The maximum value of the range.

    Returns
    -------
    int | float
        The bounded value, which is guaranteed to be within the range [min_value, max_value].
    """
    return max(min(value, max_value), min_value)
