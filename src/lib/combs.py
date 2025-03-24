from typing import TYPE_CHECKING

from lib.constants import c

if TYPE_CHECKING:
    from numpy import ndarray


def approximate_high_frequency_spectrum(f: 'ndarray', a: 'ndarray', f0: float, fs: float, fS: float,
                                        laser_wl: float = 1.645560e-6) -> 'tuple[ndarray, ndarray]':
    """
    Performs an approximate transformation of the low frequency comb into the high frequency comb,
    which contains the frequencies that actually interact with and are absorbed by the gas. This is 
    an approximation because it assumes a specific wavelength for the laser, which in real life can
    vary slightly accross measurements due to temperature/current fluctuations and other 
    perturbations. The real high frequency comb will be slightly shifted in frequency 
    from the one obtained by this function; nevertheless, this function is useful for a first 
    approximation of the high frequency comb.

    Parameters
    ----------
    f : ndarray
        Frequency array of the low-frequency comb.
    a : ndarray
        Amplitude array of the low-frequency comb.
    f0 : float
        Center frequency of the low-frequency comb.
    fs : float
        Frequency spacing of the low-frequency comb.
    fS : float
        Frequency spacing of the high-frequency comb.
    laser_wl : float
        Wavelength of the laser, should be as close as possible to the actual laser wavelength used
        in the experiment.

    Returns
    -------
    tuple[ndarray, ndarray]
        A tuple containing the frequency and amplitude arrays of the high-frequency comb.
    """
    laser_f = c / laser_wl
    return laser_f + (f - f0)*fS/fs, a


def to_frequency(wl: 'ndarray', amplitude: 'ndarray') -> 'tuple[ndarray, ndarray]':
    """
    Converts a wavelength spectrum to a frequency spectrum. The input wavelength array is assumed to
    be in nm.

    Parameters
    ----------
    wl : ndarray
        Wavelength array in nm.
    amplitude : ndarray
        Amplitude array.

    Returns
    -------
    tuple[ndarray, ndarray]
        A tuple containing the frequency and amplitude arrays.
    """
    import numpy as np
    return c * 1e9 / np.array(wl[::-1]), amplitude[::-1]


def to_wavelength(f: 'ndarray', amplitude: 'ndarray') -> 'tuple[ndarray, ndarray]':
    """
    Converts a frequency spectrum to a wavelength spectrum. The output wavelength array is in nm.

    Parameters
    ----------
    f : ndarray
        Frequency array.
    amplitude : ndarray
        Amplitude array.

    Returns
    -------
    tuple[ndarray, ndarray]
        A tuple containing the wavelength and amplitude arrays.
    """
    return to_frequency(f, amplitude)


def nan_shift(a: 'ndarray', n: int) -> 'ndarray':
    """
    Shifts the array elements by n positions, filling the empty positions with nan.

    Parameters
    ----------
    a : ndarray
        Array to shift.
    n : int
        Number of positions to shift. If n is negative, the array is shifted to the left.

    Returns
    -------
    ndarray
        Shifted array.
    """
    import numpy as np
    if n < 0:
        return np.concatenate([a[-n:], np.full(-n, np.nan)])
    return np.concatenate([np.full(n, np.nan), a[:-n]])


def normalise_transmission(f: 'ndarray', a: 'ndarray',
                           replace_outliers: bool = False) -> 'tuple[ndarray, ndarray]':
    """
    Normalises the transmission spectrum by dividing it by the maximum value of its rolling mean
    Optionally replaces any outliers with the rolling mean.

    Parameters
    ----------
    f : ndarray
        Frequency array.
    a : ndarray
        Amplitude array.
    replace_outliers : bool
        Whether to replace outliers with the rolling mean or not.

    Returns
    -------
    tuple[ndarray, ndarray]
        A tuple containing the normalised frequency and amplitude arrays.

    Note
    ----
    The rolling mean is calculated with a window of 3 points after removing outliers.
    Outliers are defined as points that deviate more than 20% up from the average of their
    nearest neighbours. This is necessary because the transmission spectrum can have spikes
    on its side wings due to noisy comb teeth that can affect the normalisation.
    """
    import numpy as np

    from lib.statistics import rolling_mean

    prv = nan_shift(a, 1)
    nxt = nan_shift(a, -1)
    neighbor_avg = (prv + nxt) / 2
    relative_diff = (a - neighbor_avg) / neighbor_avg
    relative_diff[np.isnan(relative_diff)] = 0  # set to 0 where nan:
    outliers = relative_diff > 0.2
    a_without_outliers = a[~outliers]

    if replace_outliers:
        a[outliers] = neighbor_avg[outliers]

    rolling_avg = rolling_mean(a_without_outliers, window=3)

    return f, a / max(rolling_avg)


def get_comb_frequencies(center_freq: float, freq_spacing: float, number_of_teeth: int) -> 'ndarray':
    """
    Returns an array of the frequencies of the comb teeth.

    Parameters
    ----------
    center_freq : float
        Center frequency of the comb.
    freq_spacing : float
        Frequency spacing between teeth.
    number_of_teeth : int
        Number of teeth in the comb.

    Returns
    -------
    ndarray
        Array of the frequencies of the comb teeth.
    """
    import numpy as np

    n_teeth_left = number_of_teeth // 2
    n_teeth_right = number_of_teeth - n_teeth_left

    return np.arange(center_freq - n_teeth_left * freq_spacing,
                     center_freq + n_teeth_right * freq_spacing, freq_spacing)
