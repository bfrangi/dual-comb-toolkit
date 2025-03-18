from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

    from numpy import ndarray
    from numpy.typing import ArrayLike


def interpolate_onto_common_grid(x1: 'ArrayLike', y1: 'ArrayLike', x2: 'ArrayLike',
                                 y2: 'ArrayLike', step: 'Optional[complex]' = 1000j,
                                 remove_nan: 'Optional[bool]' = True) -> 'tuple[ndarray, ndarray, ndarray, ndarray]':
    """
    Interpolates the data sets y1 and y2 onto a common grid. This is useful when the data sets have
    different x-values, and we want to compare them.

    Parameters
    ----------
    x1 : array-like
        The x-values of the first data set.
    y1 : array-like
        The y-values of the first data set.
    x2 : array-like
        The x-values of the second data set.
    y2 : array-like
        The y-values of the second data set.
    step : complex, optional
        The step size of the grid. Defaults to 1000j.
    remove_nan : bool, optional
        Whether to remove NaN values from the interpolated data. Defaults to True.

    Returns
    -------
    tuple[ndarray, ndarray, ndarray, ndarray]
        A tuple containing the x-values and y-values of the interpolated data sets. Tuple
        is of the form (new_x1, new_y1, new_x2, new_y2).
    """
    import numpy as np
    from scipy.interpolate import griddata
    x1_mod, y1_mod, x2_mod, y2_mod = np.array(x1), np.array(y1), np.array(x2), np.array(y2)

    uneven_grid_x = np.unique(np.concatenate((x1_mod, x2_mod)))
    grid_x = np.mgrid[uneven_grid_x[0]:uneven_grid_x[-1]:step]

    new_y1 = griddata(x1_mod, y1_mod, grid_x, method='cubic')
    new_y2 = griddata(x2_mod, y2_mod, grid_x, method='cubic')

    if not remove_nan:
        return grid_x, new_y1, grid_x, new_y2

    is_nan = np.isnan(new_y1)
    new_y1 = new_y1[~is_nan]
    new_x1 = grid_x[~is_nan]

    is_nan = np.isnan(new_y2)
    new_y2 = new_y2[~is_nan]
    new_x2 = grid_x[~is_nan]

    return new_x1, new_y1, new_x2, new_y2


def intersect_onto_common_grid(x1: 'ArrayLike', y1: 'ArrayLike', x2: 'ArrayLike', y2: 'ArrayLike',
                               step: complex = 1000j) -> 'tuple[ndarray, ndarray, ndarray, ndarray]':
    """
    Interpolates the data sets y1 and y2 onto a common grid only keeping the values where both data
    sets are not NaN.

    Parameters
    ----------
    x1 : array-like
        The x-values of the first data set.
    y1 : array-like
        The y-values of the first data set.
    x2 : array-like
        The x-values of the second data set.
    y2 : array-like
        The y-values of the second data set.
    step : complex, optional
        The step size of the grid. Defaults to 1000j.

    Returns
    -------
    tuple[ndarray, ndarray, ndarray, ndarray]
        A tuple containing the x-values and y-values of the interpolated data sets. Tuple
        is of the form (new_x1, new_y1, new_x2, new_y2).
    """
    import numpy as np
    new_x1, new_y1, new_x2, new_y2 = interpolate_onto_common_grid(x1, y1, x2, y2, step=step, remove_nan=False)

    is_nan_1 = np.isnan(new_y1)
    is_nan_2 = np.isnan(new_y2)

    is_nan = np.logical_or(is_nan_1, is_nan_2)

    new_y1 = new_y1[~is_nan]
    new_x1 = new_x1[~is_nan]
    new_y2 = new_y2[~is_nan]
    new_x2 = new_x2[~is_nan]

    return new_x1, new_y1, new_x2, new_y2