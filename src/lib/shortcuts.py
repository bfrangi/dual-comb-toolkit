from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

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


def fit_measurement_concentration(meas_name: str, center_freq: float, freq_spacing: float,
                                  number_of_teeth: int, laser_wavelength: float, high_freq_modulation: float,
                                  acq_freq: float, molecule: str, pressure: float, temperature: float,
                                  length: float, wl_min: float, wl_max: float, **kwargs
                                  ) -> tuple[float, 'np.ndarray', 'np.ndarray', 'np.ndarray', 'np.ndarray']:
    """
    Fit the concentration of a measurement to a simulation.

    Parameters
    ----------
    meas_name : str
        Name of the measurement.
    center_freq : float
        Center frequency of the radio frequency comb in Hz.
    freq_spacing : float
        Modulation frequency of the radio frequency comb in Hz.
    number_of_teeth : int
        Number of teeth to consider.
    laser_wavelength : float
        Approximate value of the laser wavelength in nm.
    high_freq_modulation : float
        Modulation frequency of the optical comb in Hz.
    acq_freq : float
        Acquisition frequency used in the measurement in Hz.
    molecule : str
        Name of the molecule measured and simulated.
    pressure : float
        Pressure in Pa.
    temperature : float
        Temperature of the gas in K.
    length : float
        Length of the absorption path in m.
    wl_min : float
        Minimum wavelength for the simulation in nm.
    wl_max : float
        Maximum wavelength for the simulation in nm.

    Other Parameters
    ----------------
    initial_guess : float, optional
        Initial guess for the concentration. Defaults to 0.5.

    Returns
    -------
    tuple[float, np.ndarray, np.ndarray]
        The concentration of the molecule (in VMR) and the wavelength and transmission spectrum
        for both the measurement and the simulation.
    """
    from lib.combs import to_wavelength
    from lib.fitting import ConcentrationFitter

    f = ConcentrationFitter(
        meas_name=meas_name,
        center_freq=center_freq,
        freq_spacing=freq_spacing,
        number_of_teeth=number_of_teeth,
        laser_wavelength=laser_wavelength,
        high_freq_modulation=high_freq_modulation,
        acq_freq=acq_freq,
        molecule=molecule,
        pressure=pressure,
        temperature=temperature,
        length=length,
        wl_min=wl_min,
        wl_max=wl_max,
        **kwargs,
    )
    x_meas, y_meas = to_wavelength(f.meas_freq, f.meas_amp)
    x_sim, y_sim = to_wavelength(f.sim_freq, f.sim_amp)

    return f.concentration, x_sim, y_sim, x_meas, y_meas


def get_measurement_spectrum(meas_name: str, acq_freq: 'Optional[float]' = None) -> None:
    """
    Compute the spectrum of a measurement.

    Parameters
    ----------
    meas_name : str
        The name of the measurement.
    acq_freq : float, optional
        The acquisition frequency of the measurement. If not provided, it will be read from the
        measurement file.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The frequency and amplitude of the spectrum.
    """
    from lib.analysis.time_domain import FFTCalculator
    from lib.measurements import read_lvm

    t, a = read_lvm(meas_name, acq_freq=acq_freq)
    fftc = FFTCalculator(t, a)

    return fftc.fft_x, fftc.fft_y
