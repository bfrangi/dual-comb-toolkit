from typing import TYPE_CHECKING

from lib.measurements import Measurement

if TYPE_CHECKING:
    from typing import Any, Optional

    from matplotlib import pyplot as plt
    from numpy import ndarray


class Baseline:
    """Baseline used to correct the transmission spectrum.

    Parameters
    ----------
    measurement_names : str
        The name of measurements that are used to obtain the baseline.
        The baseline is obtained by averaging the transmission spectra of 
        the given measurements.

    Keyword Arguments
    -----------------
    center_freq : float, optional
        The center frequency of the comb.
    freq_spacing : float, optional
        The frequency spacing of the comb.
    number_of_teeth : int, optional
        The number of teeth in the comb.
    laser_wavelength : float, optional
        The wavelength of the laser.
    optical_comb_spacing : float, optional
        Spacing of the optical optical comb in Hz.
    acq_freq : float, optional
        The acquisition frequency of the measurement.
    sample_tag : str
        The tag of the sample file.
    reference_tag : str
        The tag of the reference file.
    """

    def __init__(self, measurement_names: list[str], **kwargs) -> None:
        self.measurement_names = measurement_names
        self.kwargs: 'dict[str, Any]' = kwargs
        self._baseline_amps: 'list[ndarray]' = []
        self._baseline_freq: 'Optional[ndarray]' = None
        self._baseline_amp: 'Optional[ndarray]' = None
        self._std_dev: 'Optional[ndarray]' = None

        self.kwargs['normalize'] = True

    @property
    def baseline_amps(self) -> 'list[ndarray]':
        """
        The transmission spectra of the measurements used to obtain the baseline.
        """
        if self._baseline_amps is None:
            self.compute_baseline()
        return self._baseline_amps

    @property
    def baseline_amp(self) -> 'ndarray':
        """
        The averaged transmission spectrum of the measurements used to obtain 
        the baseline.
        """
        if self._baseline_amp is None:
            self.compute_baseline()
        return self._baseline_amp

    @property
    def baseline_freq(self) -> 'ndarray':
        """
        The frequency array of the averaged transmission spectrum of the 
        measurements used to obtain the baseline.
        """
        if self._baseline_freq is None:
            self.compute_baseline()
        return self._baseline_freq

    @property
    def std_dev(self) -> 'ndarray':
        """
        The standard deviation of the transmission spectra of the measurements
        used to obtain the baseline.
        """
        if self._std_dev is None:
            self._std_dev = self.compute_standard_deviation()
        return self._std_dev

    def compute_baseline(self) -> None:
        """
        Compute the baseline transmission spectrum by averaging the transmission
        spectra of the given measurements.
        """
        baseline_amp = 0

        for measurment in self.measurement_names:
            m = Measurement(measurment, **self.kwargs)
            self._baseline_amps.append(m.transmission_spectrum.y_hz)
            baseline_amp += m.transmission_spectrum.y_hz

        self._baseline_freq = m.transmission_spectrum.x_hz
        self._baseline_amp = baseline_amp / len(self.measurement_names)

    def compute_standard_deviation(self) -> 'ndarray':
        """
        Compute the standard deviation of the transmission spectra of the given
        measurements.
        """
        from numpy import std

        return std(self.baseline_amps, axis=0)
    
    def correct_transmission(self, fr: 'ndarray', tr: 'ndarray') -> 'tuple[ndarray, ndarray]':
        """
        Correct the transmission spectrum by subtracting the baseline.
        
        Parameters
        ----------
        fr : ndarray
            The frequency array of the transmission spectrum.
        tr : ndarray
            The transmission spectrum.

        Returns
        -------
        tuple[ndarray, ndarray]
            The frequency array and corrected transmission spectrum.
        """
        if len(fr) != len(self.baseline_freq):
            raise ValueError("The frequency array of the transmission spectrum " +
                             "must have the same length as the baseline frequency array.")
        
        # If all the transmission values are 1, do not correct
        if all(tr == 1):
            return fr, tr
        
        scale = tr.max()
        tr = tr / self.baseline_amp
        tr = tr / tr.max() * scale

        return fr, tr

    def generate_baseline_plot(self) -> 'plt':
        from lib.combs import to_wavelength
        from lib.plots import scatter_plot

        wl, baseline = to_wavelength(self.baseline_freq, self.baseline_amp)

        return scatter_plot(wl, baseline, title='Baseline spectrum',
                            xlabel='Wavelength [nm]', ylabel='Transmittance [-]')

    def show_baseline_plot(self) -> None:
        self.generate_baseline_plot().show()

    def generate_baselines_plot(self) -> 'plt':
        from matplotlib import pyplot as plt

        from lib.combs import to_wavelength
        from lib.plots import article_tight

        for name, amp in zip(self.measurement_names, self.baseline_amps):
            wl, baseline = to_wavelength(self.baseline_freq, amp)
            plt.scatter(wl, baseline, label=f'Baseline {name}')

        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Transmittance [-]')
        plt.title('Baselines spectrum')
        plt.legend()
        plt.tight_layout(**article_tight)
        return plt

    def show_baselines_plot(self) -> None:
        self.generate_baselines_plot().show()
    
    def generate_standard_deviation_plot(self) -> 'plt':
        from lib.combs import to_wavelength
        from lib.plots import scatter_plot

        wl, std_dev = to_wavelength(self.baseline_freq, self.std_dev)

        return scatter_plot(wl, std_dev, title='Standard deviation',
                            xlabel='Wavelength [nm]', ylabel='Standard deviation')
    
    def show_standard_deviation_plot(self) -> None:
        self.generate_standard_deviation_plot().show()
