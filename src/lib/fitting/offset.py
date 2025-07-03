from typing import TYPE_CHECKING

import numpy as np

from lib.simulations import closest_value_indices

if TYPE_CHECKING:
    from numpy import ndarray


def overlap_transmission(
    x_to_fit: "ndarray",
    y_to_fit: "ndarray",
    x_reference: "ndarray",
    y_reference: "ndarray",
) -> "ndarray":
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
    if len(x_to_fit) < 2:
        return x_to_fit

    y_to_fit_mod: "ndarray" = y_to_fit - y_to_fit.min()
    y_to_fit_mod_max = y_to_fit_mod.max()
    if y_to_fit_mod_max == 0:
        return x_to_fit
    y_to_fit_mod /= y_to_fit_mod_max

    y_reference_mod: "ndarray" = y_reference - y_reference.min()
    y_reference_mod_max = y_reference_mod.max()
    if y_reference_mod_max == 0:
        return x_to_fit
    y_reference_mod /= y_reference_mod_max

    if len(x_to_fit) >= len(x_reference):
        raise ValueError("`x_to_fit` must not be longer than `x_reference`.")

    ref_len = len(x_reference)
    ref_step = x_reference[1] - x_reference[0]
    to_fit_len = len(x_to_fit)
    to_fit_step = x_to_fit[1] - x_to_fit[0]
    max_nr_shifts = int(
        (ref_step * (ref_len - 1) - to_fit_step * (to_fit_len - 1)) // ref_step + 1
    )

    offset_x_to_fit = np.tile(x_to_fit - x_to_fit[0], (max_nr_shifts, 1)) + np.tile(
        x_reference[:max_nr_shifts, np.newaxis], (1, to_fit_len)
    )

    min_difference = np.inf
    index = None
    for i, row in enumerate(offset_x_to_fit):
        indices = closest_value_indices(x_reference, row)
        diff = ((y_reference_mod[indices] - y_to_fit_mod) ** 2).mean()
        if diff < min_difference:
            min_difference = diff
            index = i

    return x_to_fit + offset_x_to_fit[index][0] - x_to_fit[0]
