from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from numpy import ndarray


def bounded(
    value: int | float, min_value: int | float, max_value: int | float
) -> int | float:
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


def get_x_of_min_y(x: "ndarray", y: "ndarray") -> float:
    """
    Get the x value at which the y is at its minimum.

    Parameters
    ----------
    x : np.ndarray
        The x data.
    y : np.ndarray
        The y data.

    Returns
    -------
    float
        The x value at which the y is at its minimum.
    """
    from numpy import argmin

    return x[argmin(y)]
