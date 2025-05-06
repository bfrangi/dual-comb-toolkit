import warnings
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

    import matplotlib.pyplot as plt
    from numpy import ndarray

    from lib.entities import Baseline, MeasuredSpectrum


class LVMReader:
    """Reads a LabView Measurement file (.lvm) and extracts the data.

    Parameters
    ----------
    filename : str
        The name of the file to read.

    Other Parameters
    ----------------
    acq_freq : float, optional
        The acquisition frequency of the measurement. If given, it will be used instead of the one
        found in the file.
    default_acq_freq : float, optional
        The default acquisition frequency to use if not found in the file and not given as an
        argument.
    """

    def __init__(self, filename: str, **kwargs) -> None:
        import os

        from lib.defaults import ACQ_FREQ

        script_dir = os.path.dirname(__file__)

        self.filename: str = filename
        self.acq_freq: float | None = kwargs.get('acq_freq', None)
        self.data: dict | None = None

        self._full_path: str = os.path.join(script_dir, f'../../measurements/{self.filename}')
        self._default_acq_freq: float = kwargs.get('default_acq_freq', ACQ_FREQ)

    def _read(self) -> list[str]:
        with open(self._full_path, 'r') as f:
            file_lines = f.readlines()

        if not self.acq_freq:
            try:
                self.acq_freq = float(file_lines[0].split('=')[0].strip())
            except ValueError:
                pass

        if not self.acq_freq:
            warnings.warn(
                f'Acquisition frequency not found in file {self.filename}. Must be specified in ' +
                'the first line as `f=<value in Hz>`. Taking default value of ' +
                f'{self._default_acq_freq} Hz.', RuntimeWarning)
            self.acq_freq = self._default_acq_freq
            return file_lines

        return file_lines[1::]

    def _clean_data(self) -> 'ndarray':
        data = self._read()
        data = [d.replace(',', '.') for d in data]
        return [float(d.strip()) for d in data]

    def extract_data(self) -> 'dict[str, ndarray]':
        """
        Extract the data from the file.

        Returns
        -------
        dict[str, ndarray]
            The data extracted from the file. The keys are 'time' and 'amplitude'.
        """
        import numpy as np

        amplitude = np.array(self._clean_data())
        time = np.linspace(0, len(amplitude) / self.acq_freq, len(amplitude))

        self.data = {'time': time, 'amplitude': amplitude}
        return self.data


def read_lvm(filename: str, **kwargs) -> 'tuple[ndarray, ndarray]':
    """
    Read a LabView Measurement file (.lvm) and return the time series data.

    Parameters
    ----------
    filename : str
        The name of the file to read.

    Keyword Args
    ------------
    acq_freq : float, optional
        The acquisition frequency of the measurement.

    Returns
    -------
    t : ndarray
        The time array.
    amplitude : ndarray
        The amplitude array.
    """
    acq_freq = kwargs.get('acq_freq', None)
    reader = LVMReader(filename, acq_freq=acq_freq)
    data = reader.extract_data()
    return data['time'], data['amplitude']


def read_sample_and_reference_lvm(sample_filename: str, reference_filename: str,
                                  **kwargs: dict) -> 'tuple[ndarray, ndarray, ndarray]':
    """
    Read the sample and reference files and return the time series data.

    Parameters
    ----------
    sample_filename : str
        The name of the sample file.
    reference_filename : str
        The name of the reference file.

    Keyword Args
    ------------
    acq_freq : float, optional
        The acquisition frequency of the measurement.

    Returns
    -------
    t : ndarray
        The time array.
    sample_amplitude : ndarray
        The sample amplitude array.
    reference_amplitude : ndarray
        The reference amplitude array.
    """
    acq_freq = kwargs.get('acq_freq', None)
    sample_t, sample_amplitude = read_lvm(sample_filename, acq_freq=acq_freq)
    reference_t, reference_amplitude = read_lvm(reference_filename, acq_freq=acq_freq)

    # Check that time arrays are the same in both time series
    if not (sample_t == reference_t).all():
        raise ValueError('Time arrays are not the same in both time series')

    return sample_t, sample_amplitude, reference_amplitude


class Measurement:
    """Represents a measurement of a gas sample and a reference sample.

    Parameters
    ----------
    measurement_name : str
        The name of the measurement.

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
        Spacing between optical comb teeth.

    Other Parameters
    ----------------
    acq_freq : float, optional
        The acquisition frequency of the measurement.
    normalize : bool, optional
        Whether to normalize or not the transmission spectrum.
    sample_tag : str, optional
        The tag of the sample file.
    reference_tag : str, optional
        The tag of the reference file.
    baseline: list[str], optional
        The baseline to correct the measurement.
    molecule : str, optional
        The molecule measured.
    pressure : float, optional
        The pressure in Pa.
    temperature : float, optional
        The temperature in K.
    length : float, optional
        The length of the absorption path in m.    
    concentration : float, optional
        The concentration of the molecule in VMR.
    """
    time_series_properties = ['t', 'sample_amplitude', 'reference_amplitude']
    transmission_properties = ['transmission_freq', 'transmission_amp']

    def __init__(self, measurement_name: str, **kwargs) -> None:
        from lib.defaults import (
            CENTER_FREQ,
            FREQ_SPACING,
            OPTICAL_COMB_SPACING,
            LASER_WAVELENGTH,
            NUMBER_OF_TEETH,
            REFERENCE_TAG,
            SAMPLE_TAG,
        )

        # Measurement files
        self.measurement_name = measurement_name
        self._sample_tag = kwargs.get('sample_tag', SAMPLE_TAG)
        self._reference_tag = kwargs.get('reference_tag', REFERENCE_TAG)
        self.sample_filename = f'{self.measurement_name}-{self._sample_tag}.lvm'
        self.reference_filename = f'{self.measurement_name}-{self._reference_tag}.lvm'
        self.acq_freq = kwargs.get('acq_freq', None)

        # Time series data
        self._t = None
        self._sample_amplitude = None
        self._reference_amplitude = None

        # Transmission data
        self._transmission_freq: 'Optional[ndarray]' = None
        self._transmission_amp: 'Optional[ndarray]' = None
        self._transmission_spectrum: 'Optional[MeasuredSpectrum]' = None

        # Transmission analysis parameters
        self.center_freq: float = kwargs.get('center_freq', CENTER_FREQ)
        self.freq_spacing: float = kwargs.get('freq_spacing', FREQ_SPACING)
        self.number_of_teeth: int = kwargs.get('number_of_teeth', NUMBER_OF_TEETH)
        self.laser_wavelength: float = kwargs.get('laser_wavelength', LASER_WAVELENGTH)
        self.optical_comb_spacing: float = kwargs.get('optical_comb_spacing', OPTICAL_COMB_SPACING)
        self.normalize: bool = kwargs.get('normalize', True)
        self.baseline: 'Optional[Baseline]' = kwargs.get('baseline', None)

        # Other parameters
        self.molecule: str = kwargs.get('molecule', None)
        self.pressure: float = kwargs.get('pressure', None)
        self.temperature: float = kwargs.get('temperature', None)
        self.length: float = kwargs.get('length', None)
        self.concentration: float = kwargs.get('concentration', None)


    @property
    def kwargs(self) -> dict:
        return {
            'center_freq': self.center_freq,
            'freq_spacing': self.freq_spacing,
            'number_of_teeth': self.number_of_teeth,
            'laser_wavelength': self.laser_wavelength,
            'optical_comb_spacing': self.optical_comb_spacing
        }

    def get(self, property_name: str) -> 'ndarray':
        property_name = property_name.lower().replace(' ', '_').strip('_ ')
        if property_name not in self.time_series_properties + self.transmission_properties:
            raise ValueError(f'Property {property_name} not found.')

        if getattr(self, f'_{property_name}') is None:
            if property_name in self.time_series_properties:
                self._t, self._sample_amplitude, self._reference_amplitude = \
                    read_sample_and_reference_lvm(self.sample_filename, self.reference_filename,
                                                  acq_freq=self.acq_freq)
            else:
                self._compute_transmission()
        return getattr(self, f'_{property_name}')

    # Time series properties

    @property
    def t(self) -> 'ndarray':
        """Time array of the sample and reference signals."""
        return self.get('t')

    @property
    def sample_signal(self) -> 'ndarray':
        """Amplitude of the sample signal."""
        return self.get('sample_amplitude')

    @property
    def reference_signal(self) -> 'ndarray':
        """Amplitude of the reference signal."""
        return self.get('reference_amplitude')

    # Transmission properties

    @property
    def transmission_spectrum(self) -> 'MeasuredSpectrum':
        if self._transmission_spectrum is None:
            self._compute_transmission()
        return self._transmission_spectrum

    def _compute_transmission(self) -> None:
        from lib.analysis import TransmissionAnalyser
        from lib.combs import normalize_transmission
        from lib.entities import MeasuredSpectrum

        ta = TransmissionAnalyser(self.t, self.sample_signal,
                                  self.reference_signal, **self.kwargs)
        tr_freq, tr_amp = ta.transmission_freq, ta.transmission_amp

        if self.normalize:
            tr_freq, tr_amp = \
                normalize_transmission(tr_freq, tr_amp, replace_outliers=False)

        if self.baseline:
            tr_freq, tr_amp = self.baseline.correct_transmission(tr_freq, tr_amp)

        self._transmission_spectrum = MeasuredSpectrum(tr_freq, tr_amp, xu='Hz',
                                                       meas_name=self.measurement_name,
                                                       center_freq=self.center_freq,
                                                       freq_spacing=self.freq_spacing,
                                                       number_of_teeth=self.number_of_teeth,
                                                       laser_wavelength=self.laser_wavelength,
                                                       optical_comb_spacing=self.optical_comb_spacing,
                                                       acq_freq=self.acq_freq,
                                                       molecule=self.molecule,
                                                       pressure=self.pressure,
                                                       temperature=self.temperature,
                                                       length=self.length,
                                                       concentration=self.concentration)

    # Plots

    def generate_time_series_plot(self) -> 'plt':
        import matplotlib.pyplot as plt

        plt.plot(self.t, self.sample_signal, label='Sample')
        plt.plot(self.t, self.reference_signal, label='Reference')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.title('Time Series')
        return plt

    def show_time_series_plot(self) -> None:
        self.generate_time_series_plot().show()

    def generate_transmission_plot(self, wu: str = 'nm') -> 'plt':
        import matplotlib.pyplot as plt

        from lib.combs import to_wavelength

        if wu == 'nm':
            freq, trans = to_wavelength(self.transmission_freq, self.transmission_amp)

        plt.scatter(freq, trans)
        plt.xlabel('Frequency (Hz)' if wu == 'hz' else 'Wavelength (nm)')
        plt.ylabel('Transmission')
        plt.title('Transmission Spectrum')
        return plt

    def show_transmission_plot(self, wu: str = 'nm') -> None:
        self.generate_transmission_plot(wu).show()
