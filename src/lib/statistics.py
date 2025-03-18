from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from numpy import ndarray
    from numpy.typing import ArrayLike
    from pandas import Series


def rolling_mean(s: 'ArrayLike', window: int) -> 'Series':
    """Computes the rolling mean of a pandas Series.

    Parameters:
    ----------
    s : array-like
        The input series.
    window : int
        The window size for the rolling mean.

    Returns
    -------
    pandas.Series
        The rolling mean of the input series.
    """
    import pandas as pd

    s = pd.Series(s)
    return s.rolling(window=window, min_periods=1, center=True).mean()


def mean_square_difference(y1: 'ndarray', y2: 'ndarray') -> float:
    """
    Computes the mean square difference between y1 and y2.

    Parameters
    ----------
    y1 : ndarray
        The first array.
    y2 : ndarray
        The second array.

    Returns
    -------
    float
        The mean square difference between y1 and y2.

    Raises
    ------
    ValueError
        If the two arrays do not have the same length.

    Notes
    -----
    The mean square difference is computed as:

    .. math:: \\frac{1}{N} \\sum_{i=1}^{N} (y1_i - y2_i)^2

    where :math:`N` is the length of the arrays.
    """
    if len(y1) != len(y2):
        raise ValueError('The two arrays must have the same length')
    import numpy as np
    return np.sum((y1 - y2)**2) / len(y1)
