from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


def simulate_line(molecule: str, wl_min: float, wl_max: float,
                  **conditions: dict[str, float]) -> 'tuple[np.ndarray, np.ndarray]':
    """
    Simulate a line spectrum.

    Parameters
    ----------
    molecule : str
        The name of the molecule.
    wl_min : float
        The minimum wavelength in nm.
    wl_max : float
        The maximum wavelength in nm.

    Keyword Arguments
    -----------------
    vmr : float
        The volume mixing ratio of the molecule in air.
    pressure : float
        The pressure in Pa.
    temperature : float
        The temperature in K.
    length : float
        The length of the absorption path in m.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The wavelength and transmission spectrum.
    """
    from lib.simulations import Simulator

    for key in ['vmr', 'pressure', 'temperature', 'length']:
        if key not in conditions:
            raise ValueError(f"Missing condition: {key}.")

    s = Simulator(
        molecule=molecule,
        vmr=conditions['vmr'],
        pressure=conditions['pressure'],
        temperature=conditions['temperature'],
        length=conditions['length'],
    )
    s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max)
    wl, transmission = s.get_transmission_spectrum(wl_min, wl_max)

    return wl, transmission


def get_measurement_transmission(meas_name: str, **specifications: dict[str, float | int]) -> 'tuple[np.ndarray, np.ndarray]':
    """
    Get the transmission spectrum of a measurement.

    Parameters
    ----------
    meas_name : str
        The name of the measurement.

    Keyword Arguments
    -----------------
    center_freq : float
        The center frequency of the measurement.
    freq_spacing : float
        The frequency spacing of the measurement.
    number_of_teeth : int
        The number of teeth in the measurement.
    laser_wavelength : float
        The laser wavelength of the measurement.
    high_freq_modulation : float
        The high frequency modulation of the measurement.
    acq_freq : float
        The acquisition frequency of the measurement.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The wavelength and transmission spectrum.
    """
    from lib.combs import to_wavelength
    from lib.measurements import Measurement

    specification_names = ['center_freq', 'freq_spacing', 'number_of_teeth', 'laser_wavelength',
                            'high_freq_modulation', 'acq_freq']

    for key in specification_names:
        if key not in specifications:
            raise ValueError(f"Missing specification: {key}.")

    m = Measurement(
        meas_name,
        center_freq=specifications['center_freq'],
        freq_spacing=specifications['freq_spacing'],
        number_of_teet=specifications['number_of_teeth'],
        laser_wavelength=specifications['laser_wavelength'],
        high_freq_modulation=specifications['high_freq_modulation'],
        acq_freq=specifications['acq_freq'],
    )
    x_meas, y_meas = m.transmission_freq, m.transmission_amp
    x_meas, y_meas = to_wavelength(x_meas, y_meas)

    return x_meas, y_meas
