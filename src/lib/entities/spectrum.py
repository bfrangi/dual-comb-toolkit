from typing import TYPE_CHECKING

from lib.combs import to_wavelength, to_frequency

if TYPE_CHECKING:
    from typing import Optional

    from numpy import ndarray


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

    def __init__(self, x: 'ndarray', y: 'ndarray', xu: str = 'nm', **kwargs) -> None:

        # Properties

        self.molecule: 'Optional[str]' = kwargs.get('molecule', None)
        self.pressure: 'Optional[float]' = kwargs.get('pressure', None)
        self.temperature: 'Optional[float]' = kwargs.get('temperature', None)
        self.concentration: 'Optional[float]' = kwargs.get('concentration', None)
        self.length: 'Optional[float]' = kwargs.get('length', None)

        # Spectrum values

        if xu not in ['nm', 'Hz']:
            raise ValueError(f'Waveunit {xu} not recognized. Use "nm" or "Hz".')

        self.x_nm = x
        self.y_nm = y
        self.x_hz = x.copy()
        self.y_hz = y.copy()

        if xu == 'Hz':
            self.x_nm, self.y_nm = to_wavelength(self.x_nm, self.y_nm)
        else:
            self.x_hz, self.y_hz = to_frequency(self.x_hz, self.y_hz)

        self.xu = xu

    def plot(self, xu: str = 'nm') -> None:
        """
        Plot the spectrum in the specified units.

        Parameters
        ----------
        xu : str, optional
            Units of the x values. Defaults to 'nm'. Can be 'nm' or 'Hz'.
        """
        if xu not in ['nm', 'Hz']:
            raise ValueError(f'Waveunit {xu} not recognized. Use "nm" or "Hz.')
        
        from lib.plots import spectrum_plot

        if xu == 'Hz':
            x, y = self.x_hz, self.y_hz
            xlabel = 'Frequency [Hz]'
        else:
            x, y = self.x_nm, self.y_nm
            xlabel = 'Wavelength [nm]'
        
        spectrum_plot(x, y, 'Spectrum', xlabel=xlabel).show()


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
        Modulation frequency of the radio frequency comb in Hz.
    number_of_teeth : Optional[int]
        Number of teeth to consider.
    laser_wavelength : Optional[float]
        Approximate value of the laser wavelength in m.
    high_freq_modulation : Optional[float]
        Modulation frequency of the optical comb in Hz.
    acq_freq : Optional[float]  
        Acquisition frequency used in the measurement in Hz.

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
    def __init__(self, x: 'ndarray', y: 'ndarray', xu: 'str' = 'nm',
                 **kwargs: dict[str, str | float | int]) -> None:
        super().__init__(x, y, xu, **kwargs)
        self.meas_name: 'Optional[str]' = kwargs.get('meas_name', None)
        self.center_freq: 'Optional[float]' = kwargs.get('center_freq', None)
        self.freq_spacing: 'Optional[float]' = kwargs.get('freq_spacing', None)
        self.number_of_teeth: 'Optional[int]' = kwargs.get('number_of_teeth', None)
        self.laser_wavelength: 'Optional[float]' = kwargs.get('laser_wavelength', None)
        self.high_freq_modulation: 'Optional[float]' = kwargs.get('high_freq_modulation', None)
        self.acq_freq: 'Optional[float]' = kwargs.get('acq_freq', None)


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
    def __init__(self, x: 'ndarray', y: 'ndarray', xu: 'str' = 'nm',
                 **kwargs: dict[str, str | float]) -> None:
        super().__init__(x, y, xu, **kwargs)
        self.wl_min: 'Optional[float]' = kwargs.get('wl_min', None)
        self.wl_max: 'Optional[float]' = kwargs.get('wl_max', None)
