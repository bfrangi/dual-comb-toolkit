from typing import TYPE_CHECKING

from lib.combs import to_frequency, to_wavelength

from numpy import ndarray

if TYPE_CHECKING:
    from typing import Optional

    from matplotlib import pyplot as plt


class Spectrum:
    """
    A spectrum with x and y values. The x values are always stored in nm,
    but can be given in nm or Hz.

    Parameters
    ----------
    x : ndarray
        X values of the spectrum in `xu` units.
    y : ndarray
        Transmission values.
    xu : str, optional
        Unit of the given values. Defaults to 'nm'. Can be 'nm' or 'Hz'.

    Other Parameters
    ----------------
    molecule : Optional[str]
        Name of the molecule.
    pressure : Optional[float]
        Pressure in Pa.
    temperature : Optional[float]
        Temperature in K.
    concentration : Optional[float]
        Concentration of the molecule.
    length : Optional[float]
        Length of the absorption path in m.
    """

    def __init__(self, x: "ndarray", y: "ndarray", xu: str = "nm", **kwargs) -> None:
        # Properties

        self.molecule: "Optional[str]" = kwargs.get("molecule", None)
        self.pressure: "Optional[float]" = kwargs.get("pressure", None)
        self.temperature: "Optional[float]" = kwargs.get("temperature", None)
        self.concentration: "Optional[float]" = kwargs.get("concentration", None)
        self.length: "Optional[float]" = kwargs.get("length", None)

        # Spectrum values

        if xu not in ["nm", "Hz"]:
            raise ValueError(f'Waveunit {xu} not recognized. Use "nm" or "Hz".')

        self.x_nm = x
        self.y_nm = y
        self.x_hz = x.copy()
        self.y_hz = y.copy()

        if xu == "Hz":
            self.x_nm, self.y_nm = to_wavelength(self.x_nm, self.y_nm)
        else:
            self.x_hz, self.y_hz = to_frequency(self.x_hz, self.y_hz)

        self.xu = xu

    def data(self, xu: str = "nm") -> "tuple[ndarray, ndarray]":
        """
        Get the x and y values in the specified units.

        Parameters
        ----------
        xu : str, optional
            Units of the x values. Defaults to 'nm'. Can be 'nm' or 'Hz'.

        Returns
        -------
        x : ndarray
            X values of the spectrum in `xu` units.
        y : ndarray
            Transmission values.
        """
        xu_lower = xu.lower()
        if xu_lower not in ["nm", "hz"]:
            raise ValueError(f'Waveunit {xu} not recognized. Use "nm" or "Hz.')

        if xu_lower == "hz":
            return self.x_hz, self.y_hz
        return self.x_nm, self.y_nm

    def generate_plot(self, xu: str = "nm") -> "plt":
        """
        Plot the spectrum in the specified units.

        Parameters
        ----------
        xu : str, optional
            Units of the x values. Defaults to 'nm'. Can be 'nm' or 'Hz'.
        """
        from lib.plots import spectrum_plot

        x, y = self.data(xu)

        if xu == "Hz":
            xlabel = "Frequency [Hz]"
        else:
            xlabel = "Wavelength [nm]"

        return spectrum_plot(x, y, "Spectrum", xlabel=xlabel)
    
    def show_plot(self, xu: str = "nm") -> None:
        """
        Show the spectrum plot in the specified units.

        Parameters
        ----------
        xu : str, optional
            Units of the x values. Defaults to 'nm'. Can be 'nm' or 'Hz'.
        """
        return self.generate_plot(xu).show()

    def scale_by(self, factor: float) -> None:
        """
        Scale the y values of the spectrum by a factor.

        Parameters
        ----------
        factor : float
            Factor to scale the y values by.
        """
        self.y_nm *= factor
        self.y_hz *= factor


class MeasuredSpectrum(Spectrum):
    """
    A measured spectrum with additional metadata. The x values are always stored in nm,
    but can be given in nm or Hz.

    Parameters
    ----------
    x : ndarray
        X values of the spectrum in `xu` units.
    y : ndarray
        Transmission values.
    xu : str, optional
        Unit of the x values. Defaults to 'nm'. Can be 'nm' or 'Hz'.
    meas_name : Optional[str]
        Name of the measurement.
    center_freq : Optional[float]
        Center frequency of the radio frequency comb in Hz.
    freq_spacing : Optional[float]
        Spacing of the radio frequency comb in Hz.
    number_of_teeth : Optional[int]
        Number of teeth to consider.
    laser_wavelength : Optional[float]
        Approximate value of the laser wavelength in m.
    optical_comb_spacing : Optional[float]
        Spacing of the optical optical comb in Hz.
    acq_freq : Optional[float]
        Acquisition frequency used in the measurement in Hz.
    y_sdv : Optional[ndarray]
        Standard deviation of the y values in the same units as `y`.

    Other Parameters
    ----------------
    molecule : Optional[str]
        Name of the molecule.
    pressure : Optional[float]
        Pressure in Pa.
    temperature : Optional[float]
        Temperature in K.
    concentration : Optional[float]
        Concentration of the molecule.
    length : Optional[float]
        Length of the absorption path in m.
    """

    def __init__(
        self,
        x: "ndarray",
        y: "ndarray",
        xu: "str" = "nm",
        **kwargs: dict[str, str | float | int],
    ) -> None:
        super().__init__(x, y, xu, **kwargs)
        self.meas_name: "Optional[str]" = kwargs.get("meas_name", None)
        self.center_freq: "Optional[float]" = kwargs.get("center_freq", None)
        self.freq_spacing: "Optional[float]" = kwargs.get("freq_spacing", None)
        self.number_of_teeth: "Optional[int]" = kwargs.get("number_of_teeth", None)
        self.laser_wavelength: "Optional[float]" = kwargs.get("laser_wavelength", None)
        self.optical_comb_spacing: "Optional[float]" = kwargs.get(
            "optical_comb_spacing", None
        )
        self.acq_freq: "Optional[float]" = kwargs.get("acq_freq", None)
        
        self.y_sdv_nm: "Optional[ndarray]" = None
        self.y_sdv_hz: "Optional[ndarray]" = None

        y_sdv: "Optional[ndarray]" = kwargs.get("y_sdv", None)
        if y_sdv is not None:
            if xu == "Hz":
                self.y_sdv_hz = y_sdv
                _, self.y_sdv_nm = to_wavelength(x, y_sdv)
            else:
                self.y_sdv_nm = y_sdv
                _, self.y_sdv_hz = to_frequency(x, y_sdv)



class SimulatedSpectrum(Spectrum):
    """
    A simulated spectrum with additional metadata. The x values are always stored in nm,
    but can be given in nm or Hz.

    Parameters
    ----------
    x : ndarray
        X values of the spectrum in `xu` units.
    y : ndarray
        Transmission values.
    xu : str, optional
        Unit of the given values. Defaults to 'nm'. Can be 'nm' or 'Hz'.
    wl_min : Optional[float]
        Minimum wavelength for the simulation in nm.
    wl_max : Optional[float]
        Maximum wavelength for the simulation in nm.

    Other Parameters
    ----------------
    molecule : Optional[str]
        Name of the molecule.
    pressure : Optional[float]
        Pressure in Pa.
    temperature : Optional[float]
        Temperature in K.
    concentration : Optional[float]
        Concentration of the molecule.
    length : Optional[float]
        Length of the absorption path in m.
    """

    def __init__(
        self,
        x: "ndarray",
        y: "ndarray",
        xu: "str" = "nm",
        **kwargs: dict[str, str | float],
    ) -> None:
        super().__init__(x, y, xu, **kwargs)
        self.wl_min: "Optional[float]" = kwargs.get("wl_min", None)
        self.wl_max: "Optional[float]" = kwargs.get("wl_max", None)
