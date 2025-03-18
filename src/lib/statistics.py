from typing import TYPE_CHECKING

if TYPE_CHECKING:
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
