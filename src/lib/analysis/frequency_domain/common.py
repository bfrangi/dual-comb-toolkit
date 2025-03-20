from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

    import matplotlib.pyplot as plt
    from numpy import ndarray

__all__ = ['CombAnalyser']


class CombAnalyser:
    """Generates the frequency absorption and transmission curves of a pair of sample and
    reference comb spectra. Frequency arrays must be the same for both sample and reference data.

    Parameters
    ----------
    f_sample : ndarray
        The frequency array of the sample comb teeth in Hz.
    a_sample : ndarray
        The amplitude of the sample comb teeth.
    f_reference : ndarray
        The frequency array of the reference comb teeth in Hz.
    a_reference : ndarray
        The amplitude of the reference comb teeth.
    mean_sample_threshold : float, optional
        The threshold for the mean of the sample data. Defaults to 1e-4.
        Spectra with a mean below this threshold are considered to have no
        absorption features.
    """

    def __init__(self, f_sample: 'ndarray', a_sample: 'ndarray',
                 f_reference: 'ndarray', a_reference: 'ndarray',
                 mean_sample_threshold: float = 1e-4) -> None:

        if not (f_sample == f_reference).all():
            raise ValueError('The frequencies of the sample and reference data must be the same.')

        self.f_sample: 'ndarray' = f_sample
        self.a_sample: 'ndarray' = a_sample
        self.f_reference: 'ndarray' = f_reference
        self.a_reference: 'ndarray' = a_reference

        self._absorption_amp: 'Optional[ndarray]' = None
        self._absorption_freq: 'Optional[ndarray]' = None

        self.scaling_factor: float = 1
        self.mean_sample_threshold: float = mean_sample_threshold
    
    @property
    def absorption_amp(self) -> 'ndarray':
        if self._absorption_amp is None:
            self.get_absorption_curve()
        return self._absorption_amp

    @property
    def absorption_freq(self) -> 'ndarray':
        if self._absorption_freq is None:
            self.get_absorption_curve()
        return self._absorption_freq

    @property
    def transmission_amp(self) -> 'ndarray':
        return self.scaling_factor - self.absorption_amp

    @property
    def transmission_freq(self) -> 'ndarray':
        return self.absorption_freq

    def get_absorption_curve(self) -> 'tuple[ndarray, ndarray]':
        if self.a_sample.mean() < self.mean_sample_threshold:
            from numpy import zeros_like
            self._absorption_amp = zeros_like(self.a_sample) * self.scaling_factor
        else:
            self._absorption_amp = (self.a_reference - self.a_sample) / \
                self.a_reference * self.scaling_factor
        self._absorption_freq = self.f_sample
        return self._absorption_freq, self._absorption_amp

    def get_transmission_curve(self) -> 'tuple[ndarray, ndarray]':
        return self.transmission_freq, self.transmission_amp

    # Plots

    def generate_absorption_plot(self, interp: bool = False) -> 'plt':
        from lib.plots import scatter_plot

        return scatter_plot(self.absorption_freq, self.absorption_amp, interp=interp,
                            title='Absorption Spectrum', xlabel='Frequency (Hz)',
                            ylabel='Absorption')

    def show_absorption_plot(self, interp: bool = False) -> None:
        return self.generate_absorption_plot(interp=interp).show()

    def generate_transmission_plot(self, interp: bool = False) -> 'plt':
        from lib.plots import scatter_plot

        return scatter_plot(self.transmission_freq, self.transmission_amp,
                            interp=interp, title='Transmission Spectrum', xlabel='Frequency (Hz)',
                            ylabel='Transmission')

    def show_transmission_plot(self, interp: bool = False) -> None:
        return self.generate_transmission_plot(interp=interp).show()

    def generate_reference_spectrum_plot(self) -> 'plt':
        from lib.plots import stem_plot

        return stem_plot(self.f_reference, self.a_reference, title='Reference Spectrum',
                         xlabel='Frequency (Hz)', ylabel='Amplitude')

    def show_reference_spectrum_plot(self) -> None:
        return self.generate_reference_spectrum_plot().show()

    def generate_sample_spectrum_plot(self) -> 'plt':
        from lib.plots import stem_plot

        return stem_plot(self.f_sample, self.a_sample, title='Sample Spectrum',
                         xlabel='Frequency (Hz)', ylabel='Amplitude')

    def show_sample_spectrum_plot(self) -> None:
        return self.generate_sample_spectrum_plot().show()
