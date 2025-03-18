from typing import TYPE_CHECKING

from lib.fitting.interpolation import interpolate_onto_common_grid
from lib.statistics import mean_square_difference

if TYPE_CHECKING:
    from numpy import ndarray
    from numpy.typing import ArrayLike


def get_offset(x1: 'ArrayLike', y1: 'ArrayLike', x2: 'ArrayLike', y2: 'ArrayLike',
               step: complex = 1000j) -> float:
    """
    Returns the horzontal shift between the data y1 and y2. This shift is taken as the distance that
    y1 should be shifted to the right in order for the two data sets to overlap as much as possible.
    Note that for this function to work, one of the data sets must be sufficiently shorter than the
    other, so that it can be compared to as many overlapping sections of the longer data set as
    possible.

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
        The step size of the common grid built for the two data sets. Defaults to 1000j.

    Returns
    -------
    float
        The horizontal shift between the two data sets.
    """
    x1_mod, y1_mod, x2_mod, y2_mod = interpolate_onto_common_grid(x1, y1, x2, y2, step=step)
    step = x1_mod[1] - x1_mod[0]

    longest = y1_mod if len(y1_mod) >= len(y2_mod) else y2_mod
    shortest = y2_mod if len(y1_mod) >= len(y2_mod) else y1_mod

    def target(i):
        return mean_square_difference(longest[i:i+len(shortest)], shortest)

    min_distance = None
    min_i = 0
    for i in range(len(longest) - len(shortest)):
        res = target(i)
        if min_distance is None or res < min_distance:
            min_distance = res
            min_i = i

    result = step * min_i

    if len(y1_mod) >= len(y2_mod):
        return -x1_mod[0] + x2_mod[0] - result
    return -x1_mod[0] + x2_mod[0] + result


def overlap_data(x_to_fit: 'ndarray', y_to_fit: 'ArrayLike', x_reference: 'ArrayLike',
                 y_reference: 'ArrayLike', step: complex = 1000j) -> 'ndarray':
    """
    Overlaps the data y_to_fit with the reference data y_reference. This is done by shifting 
    y_to_fit horizontally so that it overlaps with y_reference as much as possible. The function
    returns the x-values of the shifted data y_to_fit (the y-values remain the same).

    Parameters
    ----------
    x_to_fit : ndarray
        The x-values of the data to be shifted.
    y_to_fit : array-like
        The y-values of the data to be shifted.
    x_reference : array-like
        The x-values of the reference data.
    y_reference : array-like
        The y-values of the reference data.
    step : complex, optional
        The step size of the common grid built for the two data sets. Defaults to 1000j.

    Returns
    -------
    ndarray
        The x-values of the shifted data.

    Example
    -------
    ```
    import matplotlib.pyplot as plt

    x1 = [1, 1.5, 2, 3, 4, 4.5, 5, 6, 7, 8, 9]
    y1 = [-10, -4, 1, 2, 3, 3.8, 4, 5, 4, 3, 2]
    x2 = [7, 8.5, 9, 10]
    y2 = [1, 2, 3, 4 ]

    x1 = overlap_data(x1, y1, x2, y2)

    plt.plot(x2, y2, label=f'shifted y2')
    plt.plot(x1, y1, label='y1')
    plt.legend()
    plt.show()
    ```
    """
    offset = get_offset(x_to_fit, y_to_fit, x_reference, y_reference, step=step)
    return x_to_fit + offset


def overlap_transmission(x_to_fit: 'ndarray', y_to_fit: 'ndarray', x_reference: 'ndarray',
                         y_reference: 'ndarray', step: complex = 1000j) -> 'ndarray':
    """
    Overlaps y_to_fit with y_reference.
    
    Parameters
    ----------
    x_to_fit : ndarray
        The x-values of the data to be shifted.
    y_to_fit : ndarray
        The y-values of the data to be shifted.
    x_reference : ndarray
        The x-values of the reference data.
    y_reference : ndarray
        The y-values of the reference data.
    step : complex, optional
        The step size of the common grid built for the two data sets. Defaults to 1000j.

    Returns
    -------
    ndarray
        The x-values of the shifted data.
    """
    y_to_fit_mod: 'ndarray' = y_to_fit - y_to_fit.min()
    if y_to_fit_mod.max() == 0:
        return x_to_fit
    y_to_fit_mod /= y_to_fit_mod.max()

    y_reference_mod: 'ndarray' = y_reference - y_reference.min()
    if y_reference_mod.max() == 0: 
        return x_to_fit
    y_reference_mod /= y_reference_mod.max()
    return overlap_data(x_to_fit, y_to_fit_mod, x_reference, y_reference_mod, step=step)
