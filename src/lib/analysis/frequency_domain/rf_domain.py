from typing import TYPE_CHECKING

from .common import CombAnalyser

if TYPE_CHECKING:
    from typing import Optional

    import matplotlib.pyplot as plt
    from numpy import ndarray

__all__ = ['RadioFrequencyCombExtractor', 'RadioFrequencyCombAnalyser']


class RadioFrequencyCombExtractor:
    """Extracts the comb teeth from a given frequency and amplitude data.

    Parameters
    ----------
    freq : ndarray
        The frequency array of the measured rf spectrum in Hz.
    amp : ndarray
        The amplitude array of the measured rf spectrum in Hz.

    Keyword Arguments
    -----------------
    center_freq : float
        The central frequency of the rf comb in Hz.
    freq_spacing : float
        The frequency spacing of the rf comb in Hz.
    number_of_teeth : int
        The number of teeth to extract.
    """

    def __init__(self, frequency: 'ndarray', amplitude: 'ndarray', **kwargs) -> None:
        from lib.defaults import CENTER_FREQ, FREQ_SPACING, NUMBER_OF_TEETH

        self.frequency: 'ndarray' = frequency
        self.amplitude: 'ndarray' = amplitude

        self.center_freq: float = kwargs.get('center_freq', CENTER_FREQ)
        self.freq_spacing: float = kwargs.get('freq_spacing', FREQ_SPACING)
        self.number_of_teeth: int = kwargs.get('number_of_teeth', NUMBER_OF_TEETH)

        self._extracted_frequency: 'Optional[ndarray]' = None
        self._extracted_amplitude: 'Optional[ndarray]' = None

    @property
    def comb_freq(self) -> 'ndarray':
        if self._extracted_frequency is None:
            self.extract_comb()
        return self._extracted_frequency

    @property
    def comb_amp(self) -> 'ndarray':
        if self._extracted_amplitude is None:
            self.extract_comb()
        return self._extracted_amplitude

    @property
    def comb_frequencies(self) -> list:
        n_teeth = min(
            int(self.center_freq // self.freq_spacing * 2), self.number_of_teeth)
        return [float(self.center_freq + (i - n_teeth//2) * self.freq_spacing)
                for i in range(n_teeth)]

    def _closest_frequency(self, f: float) -> float:
        from numpy import argmin
        return self.frequency[argmin(abs(self.frequency - f))]

    def extract_comb(self) -> 'tuple[ndarray, ndarray]':
        from numpy import array, nonzero

        approx_comb_frequencies = array(
            [self._closest_frequency(f) for f in self.comb_frequencies])
        indices = nonzero(approx_comb_frequencies[:, None] == self.frequency)[1]

        y = self.amplitude[indices]

        self._extracted_frequency = approx_comb_frequencies
        self._extracted_amplitude = y
        return self._extracted_frequency, self._extracted_amplitude

    # Plots

    def generate_raw_spectrum_plot(self) -> 'plt':
        from lib.plots import scatter_plot

        return scatter_plot(self.frequency, self.amplitude, title='Raw Spectrum',
                            xlabel='Frequency (Hz)', ylabel='Amplitude')

    def show_raw_spectrum_plot(self) -> None:
        self.generate_raw_spectrum_plot().show()

    def generate_spectrum_plot(self) -> 'plt':
        from lib.plots import stem_plot

        rep = self.freq_spacing
        f0 = self.center_freq
        title = 'Comb Spectrum showing extracted rf teeth and their amplitudes.\n' + \
                f'Central frequency ${f0}$ Hz and repetition rate ${rep}$ Hz.'

        return stem_plot(self._extracted_frequency, self._extracted_amplitude, title=title,
                         xlabel='Frequency (Hz)', ylabel='Amplitude')

    def show_spectrum_plot(self) -> None:
        self.generate_spectrum_plot().show()


class RadioFrequencyCombAnalyser(CombAnalyser):
    """Generates the radio frequency (rf) absorption and transmission curves of a pair of sample and
    reference rf comb spectra. Frequency arrays must be the same for both sample and reference data.

    Parameters
    ----------
    f_sample : ndarray 
        The frequency array of the sample rf comb teeth in Hz.
    a_sample : ndarray 
        The amplitude of the sample rf comb teeth.
    f_reference : ndarray 
        The frequency array of the reference rf comb teeth in Hz.
    a_reference : ndarray 
        The amplitude of the reference rf comb teeth.
    mean_sample_threshold : float, optional
        The threshold for the mean of the sample data. Defaults to 1e-4.
        Spectra with a mean below this threshold are considered to have no
        absorption features.
    """
    pass
