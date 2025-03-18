
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

    from numpy import ndarray

__all__ = ['TransmissionAnalyser']


class TransmissionAnalyser:
    """Analyse the transmission of a sample and reference time-series data.

    Parameters
    ----------
    t : ndarray
        The time array of the time series data.
    sample_amplitude : ndarray
        The amplitude of the sample time series data.
    reference_amplitude : ndarray
        The amplitude of the reference time series data.

    Keyword Arguments
    -----------------
    center_freq : float, optional
        The central frequency of the rf comb in Hz.
    freq_spacing : float, optional
        The frequency spacing of the rf comb in Hz.
    number_of_teeth : int, optional
        The number of teeth to extract.
    laser_wavelength : float, optional
        The wavelength of the laser in meters.
    high_freq_modulation : float, optional
        The high frequency modulation in Hz.
    """

    def __init__(self, t: 'ndarray', sample_amplitude: 'ndarray', reference_amplitude: 'ndarray',
                 **kwargs) -> None:
        from lib.defaults import (CENTER_FREQ, FREQ_SPACING,
                                  HIGH_FREQ_MODULATION, LASER_WAVELENGTH,
                                  NUMBER_OF_TEETH)

        # Time series data
        self.t: 'ndarray' = t
        self.sample_amplitude: 'ndarray' = sample_amplitude
        self.reference_amplitude: 'ndarray' = reference_amplitude

        # Transmission analysis parameters
        self.center_freq: float = kwargs.get('center_freq', CENTER_FREQ)
        self.freq_spacing: float = kwargs.get('freq_spacing', FREQ_SPACING)
        self.number_of_teeth: int = kwargs.get('number_of_teeth', NUMBER_OF_TEETH)
        self.laser_wavelength: float = kwargs.get('laser_wavelength', LASER_WAVELENGTH)
        self.high_freq_modulation: float = kwargs.get('high_freq_modulation', HIGH_FREQ_MODULATION)

        # Transmission data
        self._transmission_freq: 'Optional[ndarray]' = None
        self._transmission_amp: 'Optional[ndarray]' = None

    @property
    def kwargs(self) -> dict:
        return {
            'center_freq': self.center_freq,
            'freq_spacing': self.freq_spacing,
            'number_of_teeth': self.number_of_teeth,
            'laser_wavelength': self.laser_wavelength,
            'high_freq_modulation': self.high_freq_modulation
        }

    @property
    def transmission_freq(self) -> 'ndarray':
        if self._transmission_freq is None:
            self.compute_transmission()
        return self._transmission_freq

    @property
    def transmission_amp(self) -> 'ndarray':
        if self._transmission_amp is None:
            self.compute_transmission()
        return self._transmission_amp

    @property
    def absorption_freq(self) -> 'ndarray':
        return self.transmission_freq

    @property
    def absorption_amp(self) -> 'ndarray':
        return 1 - self.transmission_amp

    def compute_transmission(self) -> 'tuple[ndarray, ndarray]':
        from .frequency_domain.optical_domain import (
            OpticalFrequencyCombAnalyser, OpticalFrequencyCombExtractor)
        from .frequency_domain.rf_domain import RadioFrequencyCombExtractor
        from .time_domain import FFTCalculator

        kwargs = self.kwargs

        fftc_sample = FFTCalculator(self.t, self.sample_amplitude)
        fftc_reference = FFTCalculator(self.t, self.reference_amplitude)

        rfce_sample = RadioFrequencyCombExtractor(fftc_sample.fft_x, fftc_sample.fft_y, **kwargs)
        rfce_reference = RadioFrequencyCombExtractor(
            fftc_reference.fft_x, fftc_reference.fft_y, **kwargs)

        ofce_sample = OpticalFrequencyCombExtractor(
            rfce_sample.comb_freq, rfce_sample.comb_amp, **kwargs)
        ofce_reference = OpticalFrequencyCombExtractor(
            rfce_reference.comb_freq, rfce_reference.comb_amp, **kwargs)

        ofca = OpticalFrequencyCombAnalyser(
            ofce_sample.comb_freq, ofce_sample.comb_amp, ofce_reference.comb_freq,
            ofce_reference.comb_amp)
        self._transmission_freq = ofca.transmission_freq
        self._transmission_amp = ofca.transmission_amp

        return self._transmission_freq, self._transmission_amp
