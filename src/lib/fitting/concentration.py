from typing import TYPE_CHECKING

from lib.fitting import overlap_transmission

if TYPE_CHECKING:
    from typing import Optional

    from numpy import ndarray

    from lib.entities import MeasuredSpectrum, Result, SimulatedSpectrum


def fit_interpolated_concentration(meas_freq: 'ndarray', meas_amp: 'ndarray', molecule: str, 
                                   wl_min: float, wl_max: float, conditions: dict[str, float],
                                   initial_guess: float = 0.001) -> 'tuple[float, ndarray, ndarray]':
    condition_names = ['pressure', 'temperature', 'length']
    for name in condition_names:
        if name not in conditions:
            raise ValueError(f'Missing condition: {name}.')

    from numpy import sum
    from scipy.optimize import minimize

    from lib.combs import to_frequency
    from lib.fitting import intersect_onto_common_grid
    from lib.simulations import curry_interpolated_transmission_curve

    transmission_curve = curry_interpolated_transmission_curve(
        wl_min, wl_max, molecule, **conditions)

    # Fit concentration

    def f(conc, f_sample, a_sample):
        # Get the simulated curve
        wl_ref, a_ref = transmission_curve(conc)

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
    wl_ref, a_ref = transmission_curve(concentration)

    # Transform simulation data to frequency
    f_ref, a_ref = to_frequency(wl_ref, a_ref)

    return concentration, f_ref, a_ref


def fit_concentration_gpu(meas_freq: 'ndarray', meas_amp: 'ndarray', molecule: str, wl_min: float,
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

    s = Simulator(
        molecule=molecule,
        vmr=initial_guess,
        pressure=conditions['pressure'],
        temperature=conditions['temperature'],
        length=conditions['length'],
        use_gpu=True,
    )

    def f(conc, f_sample, a_sample):
        # Get the simulated curve
        s.vmr = conc
        s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max, exit_gpu=False)
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
    s.vmr = concentration
    s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max, exit_gpu=True)
    wl_ref, a_ref = s.get_transmission_spectrum(wl_min, wl_max)

    # Transform simulation data to frequency
    f_ref, a_ref = to_frequency(wl_ref, a_ref)

    return concentration, f_ref, a_ref


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
    Fit the concentration of a measured spectrum to a simulated spectrum. The measured spectrum
    object `meas_transmission` must have the properties `molecule`, `pressure`, `temperature`, and
    `length` set.

    Parameters
    ----------
    meas_transmission : MeasuredSpectrum
        Measured transmission spectrum to fit.
    wl_min : float
        Minimum wavelength for the simulation in nm.
    wl_max : float
        Maximum wavelength for the simulation in nm.

    Other Parameters
    ----------------
    initial_guess : float, optional
        Initial guess for the concentration. Defaults to 0.5.
    fitter : str, optional
        Fitter to use. Defaults to 'normal'. Possible values are 'normal', 'interp' and
        'normal_gpu'.
    """

    def __init__(self, meas_transmission: 'MeasuredSpectrum', wl_min: float,
                 wl_max: float, **kwargs: dict[str, float]) -> None:

        # Measurement parameters

        self._pre_meas_trasmission = meas_transmission

        self._meas_transmission: 'Optional[MeasuredSpectrum]' = None

        # Simulation parameters

        self.wl_min = wl_min
        self.wl_max = wl_max
        self.initial_guess: float = kwargs.get('initial_guess', 0.5)
        self.fitter: str = kwargs.get('fitter', 'normal')

        self._check_fitter()

        self._sim_transmission: 'Optional[SimulatedSpectrum]' = None

        # Conditions

        required_conditions = ['molecule', 'pressure', 'temperature', 'length']

        for condition in required_conditions:
            if getattr(meas_transmission, condition) is None:
                raise ValueError(f'Missing property: {condition}. This property must be set' +
                                 ' in the measured spectrum `meas_transmission`.')

        self.molecule = meas_transmission.molecule
        self.pressure = meas_transmission.pressure
        self.temperature = meas_transmission.temperature
        self.length = meas_transmission.length

    # Properties

    @property
    def result(self) -> 'Result':
        from lib.entities import Result

        return Result(self.measured_transmission, self.simulated_transmission)

    @property
    def measured_transmission(self) -> 'MeasuredSpectrum':
        if self._meas_transmission is None:
            self._process_measurement()
        return self._meas_transmission

    @property
    def simulated_transmission(self) -> 'SimulatedSpectrum':
        if self._sim_transmission is None:
            self._process_simulation()
        return self._sim_transmission

    @property
    def concentration(self) -> float:
        return self.simulated_transmission.concentration

    @property
    def conditions(self) -> dict[str, float]:
        return {
            'pressure': self.pressure,
            'temperature': self.temperature,
            'length': self.length,
        }

    def _check_fitter(self) -> None:
        if self.fitter not in ['normal', 'interp', 'normal_gpu']:
            raise ValueError(f'Invalid fitter: {self.fitter}. Possible values are "normal", ' +
                             '"interp" and "normal_gpu".')

    def _process_simulation(self) -> None:
        from lib.entities import SimulatedSpectrum

        self._check_fitter()

        fitter = fit_concentration

        if self.fitter == 'normal_gpu':
            fitter = fit_concentration_gpu
        elif self.fitter == 'interp':
            fitter = fit_interpolated_concentration

        concentration, sim_freq, sim_amp = fitter(
            self._pre_meas_trasmission.x_hz,
            self._pre_meas_trasmission.y_hz,
            self.molecule,
            self.wl_min,
            self.wl_max,
            self.conditions,
            self.initial_guess,
        )
        self._sim_transmission = SimulatedSpectrum(
            sim_freq,
            sim_amp,
            xu='Hz',
            molecule=self.molecule,
            pressure=self.pressure,
            temperature=self.temperature,
            concentration=concentration,
            length=self.length,
            wl_min=self.wl_min,
            wl_max=self.wl_max,
        )

    def _process_measurement(self) -> None:
        from lib.entities import MeasuredSpectrum

        meas_freq, meas_amp = self._pre_meas_trasmission.x_hz, self._pre_meas_trasmission.y_hz
        meas_freq = overlap_transmission(
            meas_freq,
            meas_amp,
            self.simulated_transmission.x_hz,
            self.simulated_transmission.y_hz,
        )
        self._meas_transmission = MeasuredSpectrum(
            meas_freq,
            meas_amp,
            xu='Hz',
            meas_name=self._pre_meas_trasmission.meas_name,
            center_freq=self._pre_meas_trasmission.center_freq,
            freq_spacing=self._pre_meas_trasmission.freq_spacing,
            number_of_teeth=self._pre_meas_trasmission.number_of_teeth,
            laser_wavelength=self._pre_meas_trasmission.laser_wavelength,
            optical_comb_spacing=self._pre_meas_trasmission.optical_comb_spacing,
            acq_freq=self._pre_meas_trasmission.acq_freq,
        )
