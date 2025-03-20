from typing import TYPE_CHECKING

from lib.combs import to_wavelength

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
        Unit of the x values. Defaults to 'nm'. Can be 'nm' or 'Hz'.
    """

    def __init__(self, x: 'ndarray', y: 'ndarray', xu: str = 'nm') -> None:
        if xu not in ['nm', 'Hz']:
            raise ValueError(f'Waveunit {xu} not recognized. Use "nm" or "Hz".')

        if xu == 'Hz':
            x, y = to_wavelength(x, y)

        self.x = x
        self.y = y
        self.xu = xu

    def plot(self):
        """
        Plot the spectrum.
        """
        from lib.plots import spectrum_plot
        spectrum_plot(self.x, self.y, 'Spectrum', xlabel='Wavelength [nm]').show()


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
    """
    def __init__(self,
                 x: 'ndarray',
                 y: 'ndarray',
                 xu: 'str' = 'nm',
                 meas_name: 'Optional[str]' = None,
                 center_freq: 'Optional[float]' = None,
                 freq_spacing: 'Optional[float]' = None,
                 number_of_teeth: 'Optional[int]' = None,
                 laser_wavelength: 'Optional[float]' = None,
                 high_freq_modulation: 'Optional[float]' = None,
                 acq_freq: 'Optional[float]' = None) -> None:
        super().__init__(x, y, xu)
        self.meas_name = meas_name
        self.center_freq = center_freq
        self.freq_spacing = freq_spacing
        self.number_of_teeth = number_of_teeth
        self.laser_wavelength = laser_wavelength
        self.high_freq_modulation = high_freq_modulation
        self.acq_freq = acq_freq


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
        Unit of the x values. Defaults to 'nm'. Can be 'nm' or 'Hz'.
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
    wl_min : Optional[float]
        Minimum wavelength for the simulation in nm.
    wl_max : Optional[float]
        Maximum wavelength for the simulation in nm.
    """
    def __init__(self,
                 x: 'ndarray',
                 y: 'ndarray',
                 xu: 'str' = 'nm',
                 molecule: 'Optional[str]' = None,
                 pressure: 'Optional[float]' = None,
                 temperature: 'Optional[float]' = None,
                 concentration: 'Optional[float]' = None,
                 length: 'Optional[float]' = None,
                 wl_min: 'Optional[float]' = None,
                 wl_max: 'Optional[float]' = None) -> None:
        super().__init__(x, y, xu)
        self.molecule = molecule
        self.pressure = pressure
        self.temperature = temperature
        self.concentration = concentration
        self.length = length
        self.wl_min = wl_min
        self.wl_max = wl_max
