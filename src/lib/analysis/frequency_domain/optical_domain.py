from typing import TYPE_CHECKING

from .common import CombAnalyser

if TYPE_CHECKING:
    from typing import Optional, Tuple

    import matplotlib.pyplot as plt
    from numpy import ndarray

__all__ = ['OpticalFrequencyCombExtractor', 'OpticalFrequencyCombAnalyser']


class OpticalFrequencyCombExtractor:
    """Obtains the optical comb teeth from the rf comb teeth.

    Parameters
    ----------
    freq : ndarray
        The frequency array of the rf comb.
    amp : ndarray
        The amplitude array of the rf comb.

    Keyword Arguments
    -----------------
    center_freq : float, optional
        The central frequency of the comb in Hz.
    freq_spacing : float, optional
        The frequency spacing of the comb in Hz.
    laser_wavelength : float, optional
        The wavelength of the laser in meters.
    optical_comb_spacing : float, optional
        Spacing of the optical optical comb in Hz.
    """

    def __init__(self, rf_domain_frequency: 'ndarray', rf_domain_amplitude: 'ndarray',
                 **kwargs) -> None:
        from lib.defaults import (
            CENTER_FREQ,
            FREQ_SPACING,
            OPTICAL_COMB_SPACING,
            LASER_WAVELENGTH,
        )

        self.rf_domain_frequency: 'ndarray' = rf_domain_frequency
        self.rf_domain_amplitude: 'ndarray' = rf_domain_amplitude

        self._optical_domain_frequency: 'Optional[ndarray]' = None
        self._optical_domain_amplitude: 'Optional[ndarray]' = None

        self.center_freq: float = kwargs.get('center_freq', CENTER_FREQ)
        self.freq_spacing: float = kwargs.get('freq_spacing', FREQ_SPACING)
        self.laser_wavelength: float = kwargs.get('laser_wavelength', LASER_WAVELENGTH)
        self.optical_comb_spacing: float = kwargs.get('optical_comb_spacing', OPTICAL_COMB_SPACING)

    @property
    def comb_freq(self) -> 'ndarray':
        if self._optical_domain_frequency is None:
            self.extract_comb()
        return self._optical_domain_frequency

    @property
    def comb_amp(self) -> 'ndarray':
        if self._optical_domain_amplitude is None:
            self.extract_comb()
        return self._optical_domain_amplitude

    def extract_comb(self) -> 'Tuple[ndarray, ndarray]':
        from lib.combs import approximate_high_frequency_spectrum

        # Transform the rf comb spectrum to the optical domain
        self._optical_domain_frequency, self._optical_domain_amplitude = approximate_high_frequency_spectrum(
            self.rf_domain_frequency, self.rf_domain_amplitude, f0=self.center_freq,
            fs=self.freq_spacing, fS=self.optical_comb_spacing, laser_wl=self.laser_wavelength)

        return self._optical_domain_frequency, self._optical_domain_amplitude

    # Plots

    def generate_spectrum_plot(self) -> 'plt':
        from lib.plots import stem_plot

        rep = self.optical_comb_spacing
        laser_wl = self.laser_wavelength
        title = 'Comb Spectrum showing extracted optical teeth and their amplitudes.\n' + \
                f'Central wavelength ${laser_wl*1e9}$ nm and repetition rate ${rep}$ Hz.'

        return stem_plot(self.extracted_frequency, self.extracted_amplitude, title=title,
                         xlabel='Frequency (Hz)', ylabel='Amplitude')

    def show_spectrum_plot(self) -> None:
        self.generate_spectrum_plot().show()


class OpticalFrequencyCombAnalyser(CombAnalyser):
    """Generates the optical frequency absorption and transmission curves of a pair of sample and
    reference optical comb spectra. Frequency arrays must be the same for both sample and reference
    data.

    Parameters
    ----------
    f_sample : ndarray 
        The frequency array of the sample optical comb teeth in Hz.
    a_sample : ndarray 
        The amplitude of the sample optical comb teeth.
    f_reference : ndarray 
        The frequency array of the reference optical comb teeth in Hz.
    a_reference : ndarray 
        The amplitude of the reference optical comb teeth.
    mean_sample_threshold : float, optional
        The threshold for the mean of the sample data. Defaults to 1e-4.
        Spectra with a mean below this threshold are considered to have no
        absorption features.
    """
    pass
