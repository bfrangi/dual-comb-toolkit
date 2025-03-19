from typing import TYPE_CHECKING

from lib.fitting import overlap_transmission

if TYPE_CHECKING:
    from typing import Optional

    from numpy import ndarray


def fit_concentration(meas_freq: 'ndarray', meas_amp: 'ndarray', molecule: str, wl_min: float,
                      wl_max: float, conditions: dict[str, float],
                      initial_guess: float = 0.001) -> 'tuple[float, ndarray, ndarray]':
    condition_names = ['pressure', 'temperature', 'length']
    for name in condition_names:
        if name not in conditions:
            raise ValueError(f'Missing condition: {name}.')

    from numpy import sum
    from scipy.optimize import minimize

    from lib.combs import to_frequency
    from lib.fitting import intersect_onto_common_grid
    from lib.simulations import Simulator

    # Fit concentration

    def f(conc, f_sample, a_sample):
        # Get the simulated curve
        s = Simulator(
            molecule=molecule,
            vmr=conc,
            pressure=conditions['pressure'],
            temperature=conditions['temperature'],
            length=conditions['length'],
        )
        s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max)
        wl_ref, a_ref = s.get_transmission_spectrum(wl_min, wl_max)

        # Transform the simulated data to frequency
        f_ref, a_ref = to_frequency(wl_ref, a_ref)

        # Shift the sample spectrum to overlap with the simulated data as much as possible
        f_sample = overlap_transmission(f_sample, a_sample, f_ref, a_ref)

        # Interpolate the sample data onto the common grid
        f_sample_com, a_sample_com, f_ref_com, a_ref_com = intersect_onto_common_grid(
            f_sample, a_sample, f_ref, a_ref)
        return sum((a_sample_com - a_ref_com)**2)

    result = minimize(f, initial_guess, args=(meas_freq, meas_amp), tol=1e-4, bounds=[(0, 1)])
    concentration = result.x[0]

    # Get the simulated curve
    s = Simulator(
        molecule=molecule,
        vmr=concentration,
        pressure=conditions['pressure'],
        temperature=conditions['temperature'],
        length=conditions['length'],
    )
    s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max)
    wl_ref, a_ref = s.get_transmission_spectrum(wl_min, wl_max)

    # Transform simulation data to frequency
    f_ref, a_ref = to_frequency(wl_ref, a_ref)

    return concentration, f_ref, a_ref


class ConcentrationFitter:
    """
    Fit the concentration of a measured spectrum to a simulated spectrum.

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
        Approximate value of the laser wavelength in m.
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
    """

    def __init__(self, meas_name: str, center_freq: float, freq_spacing: float,
                 number_of_teeth: int, laser_wavelength: float, high_freq_modulation: float,
                 acq_freq: float, molecule: str, pressure: float, temperature: float, length: float,
                 wl_min: float, wl_max: float, **kwargs) -> None:
        # Measurement parameters
        self.meas_name = meas_name
        self.center_freq = center_freq
        self.freq_spacing = freq_spacing
        self.number_of_teeth = number_of_teeth
        self.laser_wavelength = laser_wavelength
        self.high_freq_modulation = high_freq_modulation
        self.acq_freq = acq_freq

        self._meas_pre_freq: 'Optional[ndarray]' = None
        self._meas_pre_amp: 'Optional[ndarray]' = None
        self._meas_freq: 'Optional[ndarray]' = None
        self._meas_amp: 'Optional[ndarray]' = None

        # Simulation parameters
        self.molecule = molecule
        self.pressure = pressure
        self.temperature = temperature
        self.length = length
        self.wl_min = wl_min
        self.wl_max = wl_max
        self.initial_guess: float = kwargs.get('initial_guess', 0.5)

        self._sim_freq: 'Optional[ndarray]' = None
        self._sim_amp: 'Optional[ndarray]' = None

        # Results
        self._concentration: 'Optional[float]' = None

    # Properties

    @property
    def conditions(self) -> dict[str, float]:
        return {
            'pressure': self.pressure,
            'temperature': self.temperature,
            'length': self.length,
        }
    
    @property
    def meas_pre_freq(self) -> 'ndarray':
        if self._meas_pre_freq is None:
            self._pre_process_measurement()
        return self._meas_pre_freq

    @property
    def meas_pre_amp(self) -> 'ndarray':
        if self._meas_pre_amp is None:
            self._pre_process_measurement()
        return self._meas_pre_amp
        
    @property
    def concentration(self) -> float:
        if self._concentration is None:
            self._process_simulation()
        return self._concentration
    
    @property
    def sim_freq(self) -> 'ndarray':
        if self._sim_freq is None:
            self._process_simulation()
        return self._sim_freq
    
    @property
    def sim_amp(self) -> 'ndarray':
        if self._sim_amp is None:
            self._process_simulation()
        return self._sim_amp

    @property
    def meas_freq(self) -> 'ndarray':
        if self._meas_freq is None:
            self._process_measurement()
        return self._meas_freq
    
    @property
    def meas_amp(self) -> 'ndarray':
        if self._meas_amp is None:
            self._process_measurement()
        return self._meas_amp

    def _pre_process_measurement(self) -> None:
        from lib.measurements import Measurement

        m = Measurement(
            self.meas_name,
            center_freq=self.center_freq,
            freq_spacing=self.freq_spacing,
            number_of_teet=self.number_of_teeth,
            laser_wavelength=self.laser_wavelength,
            high_freq_modulation=self.high_freq_modulation,
            acq_freq=self.acq_freq,
        )
        self._meas_pre_freq, self._meas_pre_amp = m.transmission_freq, m.transmission_amp

    def _process_simulation(self) -> None:
        self._concentration, self._sim_freq, self._sim_amp = fit_concentration(
            self.meas_pre_freq,
            self.meas_pre_amp,
            self.molecule,
            self.wl_min,
            self.wl_max,
            self.conditions,
            initial_guess=self.initial_guess,
        )

    def _process_measurement(self) -> None:
        self._meas_amp = self.meas_pre_amp
        self._meas_freq = self.meas_pre_freq
        self._meas_freq = overlap_transmission(
            self.meas_pre_freq,
            self.meas_pre_amp,
            self.sim_freq,
            self.sim_amp,
        )
