import re
from typing import TYPE_CHECKING

from lib.datatypes import MeasuredSpectrum, Result, SimulatedSpectrum
from lib.fitting import ConcentrationFitter

if TYPE_CHECKING:
    from typing import Optional

    from matplotlib import pyplot as plt
    from numpy import ndarray


class Mapper:
    """
    Fit a spatial sweep of measured spectra to simulated spectra.

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
    initial_guess : dict[str, float], optional
        Initial guess for the fitting. Defaults to {'concentration': 0.5}.
    verbose : bool, optional
        Print the fitting results. Defaults to False.
    """

    def __init__(self, meas_names: list[str], center_freq: float, freq_spacing: float,
                 number_of_teeth: int, laser_wavelength: float, high_freq_modulation: float,
                 acq_freq: float, molecule: str, pressure: float, temperature: float, length: float,
                 wl_min: float, wl_max: float, **kwargs) -> None:
        # Measurement parameters
        self.meas_names = meas_names
        self.center_freq = center_freq
        self.freq_spacing = freq_spacing
        self.number_of_teeth = number_of_teeth
        self.laser_wavelength = laser_wavelength
        self.high_freq_modulation = high_freq_modulation
        self.acq_freq = acq_freq

        # Simulation parameters
        self.molecule = molecule
        self.pressure = pressure
        self.temperature = temperature
        self.length = length
        self.wl_min = wl_min
        self.wl_max = wl_max
        self.initial_guess: dict[str, float] = kwargs.get('initial_guess', {'concentration': 0.5})

        # Other parameters
        self.verbose: bool = kwargs.get('verbose', False)

        # Results
        self._concentrations: 'Optional[list[float]]' = None
        self._results: 'Optional[list[Result]]' = None
        self._positions: 'Optional[list[float]]' = None

    @property
    def concentrations(self) -> list[float]:
        """
        List of concentrations of the measured spectra.
        """
        if self._concentrations is None:
            self.map_concentration()
        return self._concentrations

    @property
    def results(self) -> list[Result]:
        """
        List of results of the fitting.
        """
        if self._results is None:
            self.map_concentration()
        return self._results

    @property
    def positions(self) -> list[tuple[float, float]]:
        """
        List of positions of the measured spectra.
        """
        if self._positions is None:
            self.map_concentration()
        return self._positions

    @property
    def concentration_map(self) -> 'ndarray':
        """
        Concentration map of the measured spectra. Format is a 2D array with the concentration
        at each position (x, y).
        """
        from numpy import array, full, nan

        if any(item is None for item in [self._concentrations, self._positions]):
            self.map_concentration()

        x = array([int(t[0]) for t in self.positions])
        y = array([int(t[1]) for t in self.positions])

        conc_map = full((x.max() + 1, y.max() + 1), nan)

        conc_map[x, y] = self.concentrations

        return conc_map

    def map_concentration(self):
        """
        Fit the concentration of a spatial sweep of measured spectra to simulated spectra.
        """
        self._concentrations = []
        self._results = []
        self._positions = []

        for meas_name in self.meas_names:
            fitter = ConcentrationFitter(
                meas_name=meas_name,
                center_freq=self.center_freq,
                freq_spacing=self.freq_spacing,
                number_of_teeth=self.number_of_teeth,
                laser_wavelength=self.laser_wavelength,
                high_freq_modulation=self.high_freq_modulation,
                acq_freq=self.acq_freq,
                molecule=self.molecule,
                pressure=self.pressure,
                temperature=self.temperature,
                length=self.length,
                wl_min=self.wl_min,
                wl_max=self.wl_max,
                initial_guess=self.initial_guess['concentration'],
            )
            meas = MeasuredSpectrum(fitter.meas_freq, fitter.meas_amp, 'Hz', meas_name,
                                    self.center_freq, self.freq_spacing, self.number_of_teeth,
                                    self.laser_wavelength, self.high_freq_modulation, self.acq_freq)
            sim = SimulatedSpectrum(fitter.sim_freq, fitter.sim_amp, 'Hz', self.molecule,
                                    self.pressure, self.temperature, fitter.concentration,
                                    self.length, self.wl_min, self.wl_max)
            self._results.append(Result(meas, sim))
            self._concentrations.append(fitter.concentration)
            self._positions.append(
                tuple(float(pos)
                      for pos in re.findall(r'\-X([0-9\.]+)\-Y([0-9\.]+)$', meas_name)[0])
            )
            if self.verbose:
                print(f'Position: {self._positions[-1]}, Concentration:' +
                      f' {self._concentrations[-1]}', end='\n\n')

    def generate_concentration_heatmap(self) -> 'plt':
        """
        Generate a heatmap of the concentration as a function of position.
        """
        import numpy as np
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        conc = self.concentration_map.T
        non_nan_rows = np.where(~np.isnan(conc).all(axis=1))[0]
        non_nan_cols = np.where(~np.isnan(conc).all(axis=0))[0]

        df = pd.DataFrame(conc)
        ax = sns.heatmap(df, cmap="crest")
        ax.invert_yaxis()
        ax.set_title("Concentration as a function of position")
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim(non_nan_cols[0], non_nan_cols[-1] + 1)
        ax.set_ylim(non_nan_rows[0], non_nan_rows[-1] + 1)

        return plt

    def show_concentration_heatmap(self) -> None:
        """
        Show the heatmap of the concentration as a function of position.
        """
        return self.generate_concentration_heatmap().show()

    def generate_concentration_plot(self, x: 'Optional[int]' = None,
                                    y: 'Optional[int]' = None) -> 'plt':
        """
        Generate a heatmap of the concentration as a function of position.
        """
        import numpy as np

        rulers = [x, y]
        if all(r is None for r in rulers):
            non_nan_rows = np.where(~np.isnan(self.concentration_map.T).all(axis=1))[0]
            if not non_nan_rows:
                raise ValueError('No non-NaN values in the concentration map.')
            y = non_nan_rows[0]
        if all(r is not None for r in rulers):
            raise ValueError('Only one of x or y can be specified.')

        from matplotlib import pyplot as plt

        if y is not None:
            conc = self.concentration_map[:, y]
            ruler = 'y'
            other = 'x'
        else:
            conc = self.concentration_map[x, :]
            ruler = 'x'
            other = 'y'

        positions = np.arange(len(conc))
        non_nan_indices = np.where(~np.isnan(conc))[0]
        

        plt.scatter(positions, conc)
        plt.xlim(non_nan_indices[0] - 1, non_nan_indices[-1] + 1)
        plt.xlabel('Position')
        plt.ylabel('Concentration [VMR]')
        plt.title(f'Concentration as a function of the {other} position ' + 
                  f'for {ruler} = {x if y is None else y}')
        return plt

    def show_concentration_plot(self, x: 'Optional[int]' = None,
                                y: 'Optional[int]' = None) -> None:
        """
        Show the heatmap of the concentration as a function of position.
        """
        return self.generate_concentration_plot(x=x, y=y).show()
