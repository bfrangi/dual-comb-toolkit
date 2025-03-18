from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from numpy import ndarray


def fft(x: 'ndarray', y: 'ndarray') -> 'tuple[ndarray, ndarray]':
    """
    Compute the Fast Fourier Transform (FFT) of a time series.

    Parameters
    ----------
    x : ndarray
        Array of x values of the series
    y : ndarray
        Array of y values of the series

    Returns
    -------
    tuple[ndarray, ndarray]
        Tuple containing the x and y values of the FFT.
    """
    from scipy.fft import fft as scipy_fft
    from scipy.fft import fftfreq

    N = len(x)  # Number of sample points
    T = x[1] - x[0]  # Sample spacing
    yf = scipy_fft(y)
    xf = fftfreq(N, T)[:N//2]
    return xf, yf


def fft_plot_data(x: 'ndarray', y: 'ndarray') -> 'tuple[ndarray, ndarray]':
    """
    Obtain the positive Fast Fourier Transform frequency and amplitude of a time series.

    Parameters
    ----------
    x : ndarray
        Array of x values of the series
    y : ndarray
        Array of y values of the series

    Returns
    -------
    tuple[ndarray, ndarray]
        Tuple containing the frequency and amplitude values of the FFT.
    """
    import numpy as np

    xf, yf = fft(x, y)
    N = len(x)  # Number of sample points
    x_data = xf
    y_data = 2.0/N * np.abs(yf[0:N//2])
    return x_data, y_data
