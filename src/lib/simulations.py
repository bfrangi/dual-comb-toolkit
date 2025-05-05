import warnings
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Callable, Optional

    from numpy import ndarray
    from radis import Spectrum

import matplotlib.pyplot as plt
import numpy as np

# Spectrum simulation ##############################################################################


def absorption_spectrum(
    wl_min: float,
    wl_max: float,
    molecule: str,
    vmr: float,
    pressure: float,
    temperature: float,
    length: float,
    isotopes: str = "1",
    wavelength_step: float = 0.001,
    database: str = "hitran",
    return_spectrum: bool = False,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, "Simulator"]:
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
    return_spectrum : bool, optional
        Whether to return the spectrum object. Default is False.

    Returns
    -------
    wavelength : np.ndarray
        The wavelength array in nm.
    absorbance : np.ndarray
        The absorption spectrum.
    spectrum : Spectrum, optional
        The spectrum object used for the calculation.
    """
    results = transmission_spectrum(
        wl_min,
        wl_max,
        molecule,
        vmr,
        pressure,
        temperature,
        length,
        isotopes,
        wavelength_step,
        database,
        return_spectrum,
    )
    if return_spectrum:
        wl, tr, simulator = results
        return wl, 1 - tr, simulator

    wl, tr = results
    return wl, 1 - tr


def absorption_spectrum_gpu(
    wl_min: float,
    wl_max: float,
    molecule: str,
    vmr: float,
    pressure: float,
    temperature: float,
    length: float,
    isotopes: str = "1",
    wavelength_step: float = 0.001,
    database: str = "hitran",
    spectrum: "Optional[Spectrum]" = None,
    exit_gpu: bool = True,
) -> tuple[np.ndarray, np.ndarray, "Optional[Spectrum]"]:
    """
    Calculate the absorption spectrum of a molecule in air using GPU.

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
    spectrum : Spectrum, optional
        The spectrum object to use. If None, a new spectrum object will be created.
    exit_gpu : bool, optional
        Whether to exit the GPU after calculation. Default is True.

    Returns
    -------
    wavelength : np.ndarray
        The wavelength array in nm.
    transmission : np.ndarray
        The transmission spectrum.
    spectrum : Spectrum
        The spectrum object used for the calculation.
    """
    wl, tr, spectrum = transmission_spectrum_gpu(
        wl_min,
        wl_max,
        molecule,
        vmr,
        pressure,
        temperature,
        length,
        isotopes,
        wavelength_step,
        database,
        spectrum,
        exit_gpu,
    )

    return wl, 1 - tr, spectrum


def transmission_spectrum(
    wl_min: float,
    wl_max: float,
    molecule: str,
    vmr: float,
    pressure: float,
    temperature: float,
    length: float,
    isotopes: str = "1",
    wavelength_step: float = 0.001,
    database: str = "hitran",
    return_spectrum: bool = False,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, "Spectrum"]:
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
    return_spectrum : bool, optional
        Whether to return the spectrum object. Default is False.

    Returns
    -------
    wavelength : np.ndarray
        The wavelength array in nm.
    transmission : np.ndarray
        The transmission spectrum.
    spectrum : Spectrum, optional
        The spectrum object used for the calculation.
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

    s = calc_spectrum(
        wn_min,
        wn_max,
        molecule=molecule,
        isotope=isotopes,
        pressure=pressure,
        Tgas=temperature,
        mole_fraction=vmr,
        path_length=length * 100,
        databank=database,
        verbose=False,
        wstep=wavenumber_step,
    )

    wn, transmittance = s.get("transmittance_noslit")

    if return_spectrum:
        return wavenumber_to_wavelength(wn)[::-1], transmittance[::-1], s

    return wavenumber_to_wavelength(wn)[::-1], transmittance[::-1]


def transmission_spectrum_gpu(
    wl_min: float,
    wl_max: float,
    molecule: str,
    vmr: float,
    pressure: float,
    temperature: float,
    length: float,
    isotopes: str = "1",
    wavelength_step: float = 0.001,
    database: str = "hitran",
    spectrum: "Optional[Spectrum]" = None,
    exit_gpu: bool = True,
) -> tuple[np.ndarray, np.ndarray, "Optional[Spectrum]"]:
    """
    Calculate the transmission spectrum of a molecule in air using GPU.

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
    spectrum : Spectrum, optional
        The spectrum object to use. If None, a new spectrum object will be created.
    exit_gpu : bool, optional
        Whether to exit the GPU after calculation. Default is True.

    Returns
    -------
    wavelength : np.ndarray
        The wavelength array in nm.
    transmission : np.ndarray
        The transmission spectrum.
    spectrum : Spectrum
        The spectrum object used for the calculation.
    """

    from lib.conversions import wavelength_to_wavenumber, wavenumber_to_wavelength

    wn_min = wavelength_to_wavenumber(wl_max)
    wn_max = wavelength_to_wavenumber(wl_min)

    pressure_bar = pressure / 1e5
    length_cm = length * 1e2

    if spectrum is None:
        import contextlib

        from radis import SpectrumFactory

        sf = SpectrumFactory(
            wn_min,
            wn_max,
            molecule=molecule,
            isotope=isotopes,
            wstep=wavelength_step,
        )

        with contextlib.redirect_stdout(None):
            sf.fetch_databank(database)

            spectrum = sf.eq_spectrum_gpu(
                Tgas=temperature,
                pressure=pressure_bar,
                mole_fraction=vmr,
                path_length=length_cm,
                exit_gpu=False,
            )
    else:
        spectrum.recalc_gpu(
            Tgas=temperature,
            pressure=pressure_bar,
            mole_fraction=vmr,
            path_length=length_cm,
        )

    wn, transmittance = spectrum.get("transmittance_noslit")

    if exit_gpu:
        from radis.gpu.gpu import gpu_exit

        gpu_exit()

    return wavenumber_to_wavelength(wn)[::-1], transmittance[::-1], spectrum


def wavelength_range(
    wl_min: float, wl_max: float, wavelength: np.ndarray, spectrum: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
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
        raise ValueError("`wl_min` must be less than `wl_max`.")
    mask = (wavelength >= wl_min) & (wavelength <= wl_max)
    return wavelength[mask], spectrum[mask]


def get_mac_interpolation_points(
    wl_min: float,
    wl_max: float,
    molecule: str,
    pressure: float,
    temperature: float,
    length: float,
    n_points: int,
    points: "Optional[list[float]]" = None,
    **conditions,
) -> "list[tuple[ndarray, ndarray]]":
    """
    Get the molar attenuation coefficient interpolation points.

    Parameters
    ----------
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.
    molecule : str
        The name of the molecule.
    pressure : float
        The pressure in Pa.
    temperature : float
        The temperature in K.
    length : float
        The length of the absorption path in m.
    n_points : int
        The number of points to use for the interpolation.
    points : list[float], optional
        The points to use for the interpolation. Default is None. If specified, takes precedence
        over `n_points`.
    conditions : dict
        The conditions for the simulation. Must include 'length', 'pressure', and
        'temperature'.

    Other Parameters
    ----------------
    isotopes : str, optional
        The isotopes of the molecule. Format: '1,2,3'.
    database : str, optional
        The database to use. Either 'hitran' or 'hitemp'.
    wavelength_step : float, optional
        The wavelength step in nm.

    Returns
    -------
    vmrs : np.ndarray
        The volume mixing ratio array.
    mac_by_parts : list[tuple[ndarray, ndarray]]
        The molar attenuation coefficient data.
    """
    if points is not None:
        n_points = None

    if n_points and n_points < 2:
        raise ValueError("`n_points` must be greater than 1.")

    if n_points:
        vmrs = np.linspace(0, 1, n_points)
    else:
        vmrs = np.array(points)

    mac_by_parts: "list[tuple[ndarray, ndarray]]" = []

    for vmr in vmrs:
        if vmr == 0:
            placeholder_vmr = 1.0e-7
            wl, tr = transmission_spectrum(
                wl_min,
                wl_max,
                molecule,
                placeholder_vmr,
                pressure,
                temperature,
                length,
                **conditions,
            )
            mac = -np.log10(tr) / (length * placeholder_vmr)
        else:
            wl, tr = transmission_spectrum(
                wl_min,
                wl_max,
                molecule,
                vmr,
                pressure,
                temperature,
                length,
                **conditions,
            )
            mac = -np.log10(tr) / (length * vmr)
        mac_by_parts.append((wl, mac))

    return vmrs, mac_by_parts


def curry_interpolated_molar_attenuation_coefficient(
    wl_min: float,
    wl_max: float,
    molecule: str,
    n_points: int = 3,
    points: "Optional[list[float]]" = None,
    **conditions,
) -> "Callable":
    """
    Return a function that computes the linearly interpolated molar attenuation coefficient
    for a given molecule and conditions.

    Parameters
    ----------
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.
    molecule : str
        The name of the molecule.
    n_points : int, optional
        The number of points to use for the interpolation. Default is 3.
    points : list[float], optional
        The points to use for the interpolation. Default is None. If specified, takes precedence
        over `n_points`.
    conditions : dict
        The conditions for the simulation. Must include 'length', 'pressure', and
        'temperature'.

    Other Parameters
    ----------------
    isotopes : str, optional
        The isotopes of the molecule. Format: '1,2,3'.
    database : str, optional
        The database to use. Either 'hitran' or 'hitemp'.
    wavelength_step : float, optional
        The wavelength step in nm.

    Returns
    -------
    interpolated_molar_attenuation_coefficient : Callable
        A function that computes the interpolated molar attenuation coefficient
        for a given concentration.
    """
    # Check parameters

    required_conditions = ["pressure", "temperature", "length"]
    for condition in required_conditions:
        if condition not in conditions:
            raise ValueError(f"Missing specification for: {condition}.")

    from lib.fitting import interpolate_onto_common_grid

    if points is not None:
        n_points = None

    length = conditions.pop("length")
    pressure = conditions.pop("pressure")
    temperature = conditions.pop("temperature")

    # Get interpolation points

    vmrs, mac_by_parts = get_mac_interpolation_points(
        wl_min,
        wl_max,
        molecule,
        pressure,
        temperature,
        length,
        n_points,
        points,
        **conditions,
    )

    # Define interpolation function

    def interpolated_molar_attenuation_coefficient(
        vmr: float,
    ) -> "tuple[ndarray, ndarray]":
        """
        Compute the interpolated molar attenuation coefficient for a given concentration.

        Parameters
        ----------
        concentration : float
            The concentration of the molecule in VMR.

        Returns
        -------
        wl : np.ndarray
            The wavelength array in nm.
        mac : np.ndarray
            The molar attenuation coefficient in m⁻¹.
        """
        if vmr < 0 or vmr > 1:
            raise ValueError("`vmr` must be between 0 and 1.")
        if vmr == 0:
            return mac_by_parts[0]
        if vmr == 1:
            return mac_by_parts[-1]

        idx_lower = np.argmax(vmrs > vmr) - 1
        idx_higher = idx_lower + 1

        wl_lower, mac_lower = mac_by_parts[idx_lower]
        wl_higher, mac_higher = mac_by_parts[idx_higher]
        vmr_lower = vmrs[idx_lower]
        vmr_higher = vmrs[idx_higher]

        wl_lower, mac_lower, wl_higher, mac_higher = interpolate_onto_common_grid(
            wl_lower, mac_lower, wl_higher, mac_higher, remove_nan=False
        )

        mac = mac_lower + (mac_higher - mac_lower) * (vmr - vmr_lower) / (
            vmr_higher - vmr_lower
        )

        mask = np.isnan(mac)
        return wl_lower[~mask], mac[~mask]

    return interpolated_molar_attenuation_coefficient


def curry_interpolated_transmission_curve(
    wl_min: float,
    wl_max: float,
    molecule: str,
    n_points: int = 5,
    points: "Optional[list[float]]" = None,
    **conditions,
) -> "Callable":
    """
    Return a function that computes the interpolated transmission curve
    for a given molecule and conditions.

    Parameters
    ----------
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.
    molecule : str
        The name of the molecule.
    n_points : int, optional
        The number of points to use for the interpolation. Default is 3.
    points : list[float], optional
        The points to use for the interpolation. Default is None. If specified, takes precedence
        over `n_points`.
    conditions : dict, optional
        The conditions for the simulation. Must include 'length', 'pressure', and
        'temperature'.

    Other Parameters
    ----------------
    isotopes : str, optional
        The isotopes of the molecule. Format: '1,2,3'.
    database : str, optional
        The database to use. Either 'hitran' or 'hitemp'.
    wavelength_step : float, optional
        The wavelength step in nm.
    """
    required_conditions = ["pressure", "temperature", "length"]
    for condition in required_conditions:
        if condition not in conditions:
            raise ValueError(f"Missing specification for: {condition}.")

    if points is not None:
        n_points = None

    length = conditions.get("length")

    interp_mac = curry_interpolated_molar_attenuation_coefficient(
        wl_min, wl_max, molecule, n_points, points, **conditions
    )

    def interpolated_transmission_curve(vmr: float) -> "tuple[np.ndarray, np.ndarray]":
        """
        Compute the interpolated transmission curve for a given concentration.

        Parameters
        ----------
        concentration : float
            The concentration of the molecule in VMR.

        Returns
        -------
        wl : np.ndarray
            The wavelength array in nm.
        transmission : np.ndarray
            The transmission spectrum.
        """
        if vmr < 0 or vmr > 1:
            raise ValueError("`vmr` must be between 0 and 1.")

        wl, mac = interp_mac(vmr)
        return wl, 10 ** (-mac * length * vmr)

    return interpolated_transmission_curve


class Simulator:
    def __init__(self, **kwargs: "dict[str, str | float]") -> None:
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
        use_gpu : bool, optional
            Whether to use GPU for the calculation. Default is False.
        spectrum: Spectrum, optional
            The spectrum object to use when `use_gpu` is True. If None given and `use_gpu` is True,
            a new spectrum object will be created.
        """
        required_kwargs = ["molecule", "vmr", "pressure", "temperature", "length"]
        for key in required_kwargs:
            if key not in kwargs:
                raise ValueError(f"Missing specification for: {key}.")

        self._molecule: str = kwargs.get("molecule")
        self.vmr: float = kwargs.get("vmr")
        self.pressure: float = kwargs.get("pressure")
        self.temperature: float = kwargs.get("temperature")
        self.length: float = kwargs.get("length")
        self._wl_min: "Optional[float]" = None
        self._wl_max: "Optional[float]" = None
        self._wavelength: "Optional[ndarray]" = None
        self._transmission: "Optional[ndarray]" = None
        self._computed_vmr = None
        self._computed_pressure = None
        self._computed_temperature = None
        self._computed_length = None
        self._isotopes: str = kwargs.get("isotopes", "1")
        self._database: str = kwargs.get("database", "hitran")
        self._wavelength_step: float = kwargs.get("wavelength_step", 0.01)
        self._use_gpu: bool = kwargs.get("use_gpu", False)
        self._spectrum: "Optional[Spectrum]" = None

    @property
    def spectrum(self) -> "Optional[Spectrum]":
        """The spectrum object used for the calculation."""
        return self._spectrum

    @property
    def molecule(self) -> str:
        """The name of the molecule."""
        return self._molecule

    @property
    def wavelength(self) -> "Optional[list]":
        """The wavelength array in nm."""
        return self._wavelength

    @property
    def transmission(self) -> "Optional[list]":
        """The transmission spectrum."""
        return self._transmission

    def compute_transmission_spectrum(
        self, wl_min: float, wl_max: float, exit_gpu: bool = True
    ) -> None:
        """
        Calculate the transmission spectrum of the molecule in air.

        Parameters
        ----------
        wl_min : float
            The minimum wavelength in nm.
        wl_max : float
            The maximum wavelength in nm.
        exit_gpu : bool, optional
            Whether to exit the GPU after calculation. Only applies if using GPU. Default is True.
        """
        sim_params = [
            self.vmr,
            self.pressure,
            self.temperature,
            self.length,
            wl_min,
            wl_max,
        ]
        sim_params_computed = [
            self._computed_vmr,
            self._computed_pressure,
            self._computed_temperature,
            self._computed_length,
            self._wl_min,
            self._wl_max,
        ]

        if all(
            param == computed_param
            for param, computed_param in zip(sim_params, sim_params_computed)
        ):
            if self._wavelength is not None and self._transmission is not None:
                return self._wavelength, self._transmission

        if self._use_gpu:
            spectrum = (
                self._spectrum
                if self._wl_min == wl_min and self._wl_max == wl_max
                else None
            )

            self._wavelength, self._transmission, self._spectrum = (
                transmission_spectrum_gpu(
                    wl_min,
                    wl_max,
                    self._molecule,
                    self.vmr,
                    self.pressure,
                    self.temperature,
                    self.length,
                    self._isotopes,
                    self._wavelength_step,
                    self._database,
                    spectrum,
                    exit_gpu=exit_gpu,
                )
            )
        else:
            self._wavelength, self._transmission, self._spectrum = (
                transmission_spectrum(
                    wl_min,
                    wl_max,
                    self._molecule,
                    self.vmr,
                    self.pressure,
                    self.temperature,
                    self.length,
                    self._isotopes,
                    self._wavelength_step,
                    self._database,
                    return_spectrum=True,
                )
            )

        self._wl_min = wl_min
        self._wl_max = wl_max
        self._computed_vmr = self.vmr
        self._computed_pressure = self.pressure
        self._computed_temperature = self.temperature
        self._computed_length = self.length

    def get_transmission_spectrum(
        self, wl_min: "Optional[float]" = None, wl_max: "Optional[float]" = None
    ) -> tuple[np.ndarray, np.ndarray]:
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
            raise ValueError("Transmission spectrum not calculated.")

        wl_min = wl_min if wl_min is not None else self._wl_min
        wl_max = wl_max if wl_max is not None else self._wl_max

        if self._wl_min > wl_min or self._wl_max < wl_max:
            raise ValueError(
                "Transmission spectrum not calculated for the specified wavelength range."
            )

        return wavelength_range(wl_min, wl_max, self._wavelength, self._transmission)

    def plot_transmission_spectrum(
        self, wl_min: "Optional[float]" = None, wl_max: "Optional[float]" = None
    ) -> plt:
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
            wl,
            transmission,
            f"Transmittance spectrum for {self._molecule} at {self.temperature} K, "
            + f"{self.pressure} Pa, {self.length} cm and {self.vmr} VMR",
            "Frequency [Hz]",
            "Transmittance [-]",
        )

    def show_transmission_spectrum(
        self, wl_min: "Optional[float]" = None, wl_max: "Optional[float]" = None
    ) -> None:
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


def closest_value_indices(array: "ndarray", values: "ndarray") -> "ndarray":
    return np.abs(array - values[:, np.newaxis]).argmin(axis=1)


def simulate_measurement(
    molecule: str, wl_min: float, wl_max: float, **kwargs: dict[str, float]
) -> "tuple[np.ndarray, np.ndarray]":
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
    vmr : float
        The concentration of the molecule in VMR.
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
        The wavelength step in nm used for the simulation. This
        limits the resolution of the simulation.
    std_dev : float, optional
        The standard deviation of the noise added to the simulated spectrum.
    number_of_teeth_for_std_dev : int, optional
        The number of teeth in the high-frequency comb used to calculate the standard deviation
        of the noise added to the simulated spectrum. If not specified, std_dev is used
        directly. If specified, the standard deviation is calculated as:
            std_dev = std_dev · number_of_teeth_for_std_dev / number_of_teeth
    scaling_std_dev : float, optional
        The standard deviation of the scaling factor applied to the simulated spectrum.
    x_shift_std_dev : float, optional
        The standard deviation of the wavelength shift applied to the simulated spectrum.
    laser_wavelength_slack : tuple(float), optional
        Range of the laser wavelength's random variation. Defaults to (-0.05, 0.05).

    References
    ----------
    Coddington, I., Newbury, N., & Swann, W. (2016). Dual-comb spectroscopy. Optica, 3(4), 414–426.
    https://doi.org/10.1364/OPTICA.3.000414
    """
    from lib.combs import get_comb_frequencies, to_frequency, to_wavelength
    from lib.constants import c
    from lib.simulations import Simulator

    if "simulator" not in kwargs:
        required_kwargs = [
            "vmr",
            "pressure",
            "temperature",
            "length",
            "laser_wavelength",
            "high_freq_modulation",
            "number_of_teeth",
        ]

        for key in required_kwargs:
            if key not in kwargs:
                raise ValueError(f"Missing specification for: {key}.")

        # Obtain comb frequencies.

        laser_wl = kwargs.get("laser_wavelength")
        laser_wl_shift_range = kwargs.get("laser_wavelength_slack", (-0.05, 0.05))
        laser_wl = np.random.uniform(
            laser_wl + laser_wl_shift_range[0], laser_wl + laser_wl_shift_range[1]
        )
        laser_freq = c / laser_wl * 1e9

        min_fr = c / wl_max * 1e9
        max_fr = c / wl_min * 1e9

        fr_sam = get_comb_frequencies(
            center_freq=laser_freq,
            freq_spacing=kwargs.get("high_freq_modulation"),
            number_of_teeth=kwargs.get("number_of_teeth"),
        )
        if fr_sam[-1] - fr_sam[0] > max_fr - min_fr:
            raise ValueError("The comb is too wide for the simulated frequency range.")

        left_overflow = np.abs(min(fr_sam[0] - min_fr, 0))
        right_overflow = np.abs(min(max_fr - fr_sam[-1], 0))
        if left_overflow + right_overflow > 0:
            warnings.warn(
                "Warning: The laser wavelength shift is too large and will be reduced "
                "to fit the range. If you require a larger range, please adjust `wl_min` and "
                "`wl_max`."
            )
        fr_sam = fr_sam + left_overflow - right_overflow

        # Simulate the transmission spectrum.

        database = kwargs.pop("database", "hitran")
        s = Simulator(molecule=molecule, **kwargs, database=database)
    else:
        s: Simulator = kwargs.get("simulator")

    s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max)
    wl_sim, tr_sim = s.get_transmission_spectrum(wl_min, wl_max)
    fr_sim, tr_sim = to_frequency(wl_sim, tr_sim)

    # Sample the simulated spectrum with a comb.

    indices = closest_value_indices(fr_sim, fr_sam)
    tr_sam = tr_sim[indices]

    # Add noise to the sampled spectrum.

    number_of_teeth_for_std_dev = kwargs.get("number_of_teeth_for_std_dev", None)
    if number_of_teeth_for_std_dev is not None:
        std_dev = (
            kwargs.get("std_dev", 0.01)
            / number_of_teeth_for_std_dev
            * kwargs.get("number_of_teeth")
        )
    else:
        std_dev = kwargs.get("std_dev", 0.01)
    noise = np.random.normal(0, std_dev, len(tr_sam))
    tr_sam += noise

    scaling_std_dev = kwargs.get("scaling_std_dev", 0.01)
    scaling_factor = max(np.abs(np.random.normal(1, scaling_std_dev, 1)[0]), 0.001)
    tr_sam *= scaling_factor

    wl_sam, tr_sam = to_wavelength(fr_sam, tr_sam)

    left_slack = wl_sam[0] - wl_sim[0]
    right_slack = wl_sim[-1] - wl_sam[-1]
    x_shift_std_dev = kwargs.get("x_shift_std_dev", 0.01)
    x_shift = min(
        max(np.random.normal(0, x_shift_std_dev, 1)[0], -left_slack), right_slack
    )
    wl_sam += x_shift

    return wl_sam, tr_sam
