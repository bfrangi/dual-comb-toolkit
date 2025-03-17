from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

import matplotlib.pyplot as plt
import numpy as np


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
    def __init__(self, molecule: str, vmr: float, pressure: float, temperature: float,
                 length: float, isotopes: str = '1', database: str = 'hitran',
                 wavelength_step: float = 0.001) -> None:
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
        """
        self._molecule = molecule
        self._vmr = vmr
        self._pressure = pressure
        self._temp = temperature
        self._length = length
        self._wl_min = None
        self._wl_max = None
        self._wavelength = None
        self._transmission = None
        self._isotopes = isotopes
        self._database = database
        self._wavelength_step = wavelength_step

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
