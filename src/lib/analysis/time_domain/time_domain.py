from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

    import matplotlib.pyplot as plt
    from numpy import ndarray


class FFTCalculator:
    """Compute the Fast Fourier Transform of a signal.

    Parameters
    ----------
    t : ndarray
        Time values.
    amplitude : ndarray
        Amplitude values.
    """
    def __init__(self, t: 'ndarray', amplitude: 'ndarray'):
        self.t: 'ndarray' = t
        self.amplitude: 'ndarray' = amplitude

        self._fft_x: 'Optional[ndarray]' = None
        self._fft_y: 'Optional[ndarray]' = None
    
    @property
    def fft_x(self) -> 'ndarray':
        if self._fft_x is None:
            self.calculate_fft()
        return self._fft_x
    
    @property
    def fft_y(self) -> 'ndarray':
        if self._fft_y is None:
            self.calculate_fft()
        return self._fft_y

    def calculate_fft(self) -> 'tuple[ndarray, ndarray]':
        from lib.fft import fft_plot_data
        self._fft_x, self._fft_y = fft_plot_data(self.t, self.amplitude)
        return self._fft_x, self._fft_y

    # Plot

    def generate_fft_plot(self) -> 'plt':
        from lib.plots import spectrum_plot
        return spectrum_plot(self.fft_x, self.fft_y, 'Spectrum of the Signal', 
                             'Frequency (Hz)', 'Amplitude')

    def show_fft_plot(self) -> None:
        self.generate_fft_plot().show()
