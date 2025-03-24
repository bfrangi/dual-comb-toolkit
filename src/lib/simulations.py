from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

    from numpy import ndarray

import matplotlib.pyplot as plt
import numpy as np

# Spectrum simulation ##############################################################################


def absorption_spectrum(wl_min, wl_max, molecule: str, vmr: float, pressure: float,
                        temperature: float, length: float, isotopes: str = '1',
                        wavelength_step: float = 0.001, database: str = 'hitran') -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the absorption spectrum of a molecule in air.

    Parameters
    ----------
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.
    molecule : str
        The name of the molecule.
    vmr : float
        The volume mixing ratio of the molecule in air.
    pressure : float
        The pressure in Pa.
    temperature : float
        The temperature in K.
    length : float
        The length of the absorption path in m.
    isotopes : str, optional
        The isotopes of the molecule. Format: '1,2,3'.
    wavelength_step : float, optional
        The wavelength step in nm.
    database : str, optional
        The database to use. Either 'hitran' or 'hitemp'.

    Returns
    -------
    wavelength : np.ndarray
        The wavelength array in nm.
    absorbance : np.ndarray
        The absorption spectrum.
    """
    wl, tr = transmission_spectrum(wl_min, wl_max, molecule, vmr, pressure, temperature,
                                   length, isotopes, wavelength_step, database)
    return wl, 1 - tr


def transmission_spectrum(wl_min, wl_max, molecule: str, vmr: float, pressure: float,
                          temperature: float, length: float, isotopes: str = '1',
                          wavelength_step: float = 0.001, database: str = 'hitran') -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the transmission spectrum of a molecule in air.

    Parameters
    ----------
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.
    molecule : str
        The name of the molecule.
    vmr : float
        The volume mixing ratio of the molecule in air.
    pressure : float
        The pressure in Pa.
    temperature : float
        The temperature in K.
    length : float
        The length of the absorption path in m.
    isotopes : str, optional
        The isotopes of the molecule. Format: '1,2,3'.
    wavelength_step : float, optional
        The wavelength step in nm.
    database : str, optional
        The database to use. Either 'hitran' or 'hitemp'.

    Returns
    -------
    wavelength : np.ndarray
        The wavelength array in nm.
    transmission : np.ndarray
        The transmission spectrum.
    """
    from radis import calc_spectrum

    from lib.conversions import (
        delta_wavelength_to_delta_wavenumber,
        pa_to_bar,
        wavelength_to_wavenumber,
        wavenumber_to_wavelength,
    )

    pressure = pa_to_bar(pressure)
    wn_min = wavelength_to_wavenumber(wl_max)
    wn_max = wavelength_to_wavenumber(wl_min)
    central_wl = (wl_min + wl_max) / 2
    wavenumber_step = delta_wavelength_to_delta_wavenumber(wavelength_step, central_wl)

    s = calc_spectrum(wn_min, wn_max, molecule=molecule, isotope=isotopes, pressure=pressure,
                      Tgas=temperature, mole_fraction=vmr, path_length=length*100,
                      databank=database, verbose=False, wstep=wavenumber_step)

    wn, transmittance = s.get('transmittance_noslit')
    return wavenumber_to_wavelength(wn)[::-1], transmittance[::-1]


def wavelength_range(wl_min: float, wl_max: float, wavelength: np.ndarray,
                     spectrum: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Filter a spectrum to a wavelength range.

    Parameters
    ----------
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.
    wl : np.ndarray
        The wavelength array in nm.
    spectrum : np.ndarray
        The spectrum array.

    Returns
    -------
    wavelength : np.ndarray
        The filtered wavelength array in nm.
    spectrum : np.ndarray
        The filtered spectrum array.
    """
    if wl_min > wl_max:
        raise ValueError('`wl_min` must be less than `wl_max`.')
    mask = (wavelength >= wl_min) & (wavelength <= wl_max)
    return wavelength[mask], spectrum[mask]


class Simulator:
    def __init__(self, **kwargs: 'dict[str, str | float]') -> None:
        """Simulator for computing the transmission spectrum of a molecule.

        Parameters
        ----------
        molecule : str
            The name of the molecule.
        vmr : float
            The volume mixing ratio of the molecule in air.
        pressure : float
            The pressure in Pa.
        temperature : float
            The temperature in K.
        length : float
            The length of the absorption path in m.

        Other Parameters
        ----------------
        isotopes : str, optional
            The isotopes of the molecule. Format: '1,2,3'.
        database : str, optional
            The database to use. Either 'hitran' or 'hitemp'.
        wavelength_step : float, optional
            The wavelength step in nm.
        """
        required_kwargs = ['molecule', 'vmr', 'pressure', 'temperature', 'length']
        for key in required_kwargs:
            if key not in kwargs:
                raise ValueError(f"Missing specification for: {key}.")

        self._molecule: str = kwargs.get('molecule')
        self._vmr: float = kwargs.get('vmr')
        self._pressure: float = kwargs.get('pressure')
        self._temp: float = kwargs.get('temperature')
        self._length: float = kwargs.get('length')
        self._wl_min: 'Optional[float]' = None
        self._wl_max: 'Optional[float]' = None
        self._wavelength: 'Optional[ndarray]' = None
        self._transmission: 'Optional[ndarray]' = None
        self._isotopes: str = kwargs.get('isotopes', '1')
        self._database: str = kwargs.get('database', 'hitran')
        self._wavelength_step: float = kwargs.get('wavelength_step', 0.001)

    @property
    def molecule(self) -> str:
        """The name of the molecule."""
        return self._molecule

    @property
    def vmr(self) -> float:
        """The volume mixing ratio of the molecule in air."""
        return self._vmr

    @property
    def pressure(self) -> float:
        """The pressure in Pa."""
        return self._pressure

    @property
    def temperature(self) -> float:
        """The temperature in K."""
        return self._temp

    @property
    def length(self) -> float:
        """The length of the absorption path in m."""
        return self._length

    @property
    def wavelength(self) -> 'Optional[list]':
        """The wavelength array in nm."""
        return self._wavelength

    @property
    def transmission(self) -> 'Optional[list]':
        """The transmission spectrum."""
        return self._transmission

    def compute_transmission_spectrum(self, wl_min: float, wl_max: float) -> None:
        """
        Calculate the transmission spectrum of the molecule in air.

        Parameters
        ----------
        wl_min : float
            The minimum wavelength in nm.
        wl_max : float
            The maximum wavelength
        """
        if wl_min == self._wl_min and wl_max == self._wl_max:
            if self._wavelength is not None and self._transmission is not None:
                return self._wavelength, self._transmission

        self._wl_min = wl_min
        self._wl_max = wl_max

        self._wavelength, self._transmission = transmission_spectrum(
            self._wl_min, self._wl_max, self._molecule, self._vmr, self._pressure,
            self._temp, self._length, self._isotopes, self._wavelength_step, self._database)

    def get_transmission_spectrum(self, wl_min: 'Optional[float]' = None,
                                  wl_max: 'Optional[float]' = None) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the transmission spectrum of the molecule in air.

        Parameters
        ----------
        wl_min : float, optional
            The minimum wavelength in nm.
        wl_max : float, optional
            The maximum wavelength in nm.

        Returns
        -------
        wavelength : np.ndarray
            The wavelength array in nm.
        transmission : np.ndarray
            The transmission spectrum.
        """
        if self._wavelength is None or self._transmission is None:
            raise ValueError('Transmission spectrum not calculated.')

        wl_min = wl_min if wl_min is not None else self._wl_min
        wl_max = wl_max if wl_max is not None else self._wl_max

        if self._wl_min > wl_min or self._wl_max < wl_max:
            raise ValueError(
                'Transmission spectrum not calculated for the specified wavelength range.')

        return wavelength_range(wl_min, wl_max, self._wavelength, self._transmission)

    def plot_transmission_spectrum(self, wl_min: 'Optional[float]' = None,
                                   wl_max: 'Optional[float]' = None) -> plt:
        """
        Plot the transmission spectrum of the molecule in air.

        Parameters
        ----------
        wl_min : float, optional
            The minimum wavelength in nm.
        wl_max : float, optional
            The maximum wavelength in nm.

        Returns
        -------
        plt
            The matplotlib plot of the transmission spectrum.
        """
        from lib.plots import spectrum_plot

        wl, transmission = self.get_transmission_spectrum(wl_min, wl_max)

        return spectrum_plot(
            wl, transmission,
            f'Transmittance spectrum for {self._molecule} at {self._temp} K, ' +
            f'{self._pressure} Pa, {self._length} cm and {self.vmr} VMR',
            'Frequency [Hz]',
            'Transmittance [-]')

    def show_transmission_spectrum(self, wl_min: 'Optional[float]' = None,
                                   wl_max: 'Optional[float]' = None) -> None:
        """
        Show the transmission spectrum of the molecule in air.

        Parameters
        ----------
        wl_min : float, optional
            The minimum wavelength in nm.
        wl_max : float, optional
            The maximum wavelength in nm.
        """
        return self.plot_transmission_spectrum(wl_min, wl_max).show()


# Measurement simulation ###########################################################################


def closest_value_indices(array: 'ndarray', values: 'ndarray') -> 'ndarray':
    return np.abs(array-values[:, np.newaxis]).argmin(axis=1)


def simulate_measurement(molecule: str, wl_min: float, wl_max: float,
                         **kwargs: dict[str, float]) -> 'tuple[np.ndarray, np.ndarray]':
    """
    Simulate a measurement.

    Parameters
    ----------
    molecule : str
        The name of the molecule.
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.
    laser_wavelength : float
        The wavelength of the laser in nm.
    high_freq_modulation : float
        The frequency spacing of the high-frequency comb in Hz.
    number_of_teeth : int
        The number of teeth in the high-frequency comb.

    Other Parameters
    ----------------
    isotopes : str, optional
        The isotopes of the molecule. Format: '1,2,3'.
    database : str, optional
        The database to use. Either 'hitran' or 'hitemp'.
    wavelength_step : float, optional
        The wavelength step in nm used for the simulation. This 
        limits the resolution of the simulation.
    std_dev : float, optional
        The standard deviation of the noise added to the simulated spectrum.
    """
    from lib.combs import get_comb_frequencies, to_frequency, to_wavelength
    from lib.constants import c
    from lib.simulations import Simulator

    required_kwargs = ['vmr', 'pressure', 'temperature', 'length', 'laser_wavelength',
                       'high_freq_modulation', 'number_of_teeth']

    for key in required_kwargs:
        if key not in kwargs:
            raise ValueError(f"Missing specification for: {key}.")

    # Validate parameters.

    laser_wl = kwargs.get('laser_wavelength')
    if wl_min > laser_wl or wl_max < laser_wl:
        raise ValueError("The laser wavelength must be within the specified wavelength range.")

    # Simulate the transmission spectrum.

    database = kwargs.pop('database', 'hitran')
    s = Simulator(
        molecule=molecule,
        **kwargs,
        database=database,
    )
    s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max)
    wl_sim, tr_sim = s.get_transmission_spectrum(wl_min, wl_max)
    fr_sim, tr_sim = to_frequency(wl_sim, tr_sim)

    # Sample the simulated spectrum with a comb.

    laser_freq = c / laser_wl * 1e9
    fr_sam = get_comb_frequencies(
        center_freq=laser_freq,
        freq_spacing=kwargs.get('high_freq_modulation'),
        number_of_teeth=kwargs.get('number_of_teeth')
    )
    indices = closest_value_indices(fr_sim, fr_sam)
    tr_sam = tr_sim[indices]

    # Add noise to the sampled spectrum.
    noise_mean = 0
    std_dev = kwargs.get('std_dev', 0.01)
    noise = np.random.normal(noise_mean, std_dev, len(tr_sam))
    tr_sam += noise

    return to_wavelength(fr_sam, tr_sam)

    # import numpy as np

    # from lib.combs import to_frequency
    # from lib.constants import c

    # laser_wl = kwargs['laser_wavelength']
    # if wl_min > laser_wl or wl_max < laser_wl:
    #     raise ValueError("The laser wavelength must be within the specified wavelength range.")
    # laser_freq = c / laser_wl * 1e9

    # optical_comb_frequencies = get_comb_frequencies(
    #     center_freq=laser_freq,
    #     freq_spacing=kwargs['high_freq_modulation'],
    #     number_of_teeth=kwargs['number_of_teeth']
    # )
    # wl, tr = simulate_line(molecule, wl_min, wl_max, **kwargs)
    # fr, tr = to_frequency(wl, tr)

    # indices = closest_value_indices(fr, optical_comb_frequencies)
    # tr_sampled = tr[indices]

    # rf_comb_frequencies = get_comb_frequencies(
    #     center_freq=kwargs['center_freq'],
    #     freq_spacing=kwargs['freq_spacing'],
    #     number_of_teeth=kwargs['number_of_teeth']
    # )

    # power_per_tooth = kwargs['comb_power'] / kwargs['number_of_teeth']

    # # Calculate signal power and convert to dB
    # sig_avg_watts = power_per_tooth
    # sig_avg_db = 10 * np.log10(sig_avg_watts)

    # # Calculate noise
    # noise_avg_db = sig_avg_db - snr_db
    # noise_avg_watts = 10 ** (noise_avg_db / 10)

    # # Generate an sample of white noise
    # mean_noise = 0
    # noise_volts = np.random.normal(mean_noise, np.sqrt(noise_avg_watts), len(tr_sampled))
    # # Noise up the original signal
    # tr_sampled += noise_volts

    # return rf_comb_frequencies, tr_sampled
