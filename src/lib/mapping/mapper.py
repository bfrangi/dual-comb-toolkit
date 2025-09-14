import os
import re
from typing import TYPE_CHECKING

from lib.entities import MeasuredSpectrum, Result
from lib.files import get_figures_path
from lib.fitting import ConcentrationFitter
from lib.math import bounded

if TYPE_CHECKING:
    from typing import Optional

    from matplotlib import pyplot as plt
    from numpy import ndarray


class Mapper:
    """
    Fit a spatial sweep of measured spectra to simulated spectra.

    Parameters
    ----------
    meas_transmissions : list[MeasuredSpectrum]
        List of measured transmission spectra to fit.
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
    database : str, optional
        Database to use for the simulation. Defaults to 'hitran' (defined in `lib.defaults`).
    spectrum_plot_folder : str, optional
        Folder to save the spectrum plots. Defaults to None, so no plots are saved.
    lower_bound : float, optional
        Lower bound for the concentration. Defaults to 0.
    upper_bound : float, optional
        Upper bound for the concentration. Defaults to 1.
    verbose : bool, optional
        Whether to print the results of the fitting. Defaults to False.
    """

    def __init__(
        self,
        meas_transmissions: "list[MeasuredSpectrum]",
        wl_min: float,
        wl_max: float,
        **kwargs: dict[str, float],
    ) -> None:
        from lib.defaults import DATABASE

        # Measurement parameters
        self.meas_transmissions = meas_transmissions

        # Simulation parameters
        self.wl_min = wl_min
        self.wl_max = wl_max
        self.initial_guess: dict[str, float] = kwargs.get("initial_guess", 0.5)
        self.lower_bound: float = bounded(kwargs.get("lower_bound", 0.0), 0, 1)
        self.upper_bound: float = bounded(kwargs.get("upper_bound", 1.0), 0, 1)
        self.database: str = kwargs.get("database", DATABASE)

        if self.lower_bound >= self.upper_bound:
            raise ValueError("Lower bound must be less than upper bound.")

        # Other parameters
        self.verbose: bool = kwargs.get("verbose", False)
        self.fitter: str = kwargs.get("fitter", "normal")

        # Results
        self._concentrations: "Optional[list[float]]" = None
        self._results: "Optional[list[Result]]" = None
        self._positions: "Optional[list[float]]" = None
        self._spectrum_plot_folder: "Optional[str]" = kwargs.get(
            "spectrum_plot_folder", None
        )

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
    def concentration_map(self) -> "ndarray":
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

    def map_concentration(self) -> None:
        """
        Fit the concentration of a spatial sweep of measured spectra to simulated spectra.
        """
        self._concentrations = []
        self._results = []
        self._positions = []

        for meas_transmission in self.meas_transmissions:
            fitter = ConcentrationFitter(
                meas_transmission,
                self.wl_min,
                self.wl_max,
                initial_guess=self.initial_guess,
                lower_bound=self.lower_bound,
                upper_bound=self.upper_bound,
                verbose=self.verbose,
                fitter=self.fitter,
                database=self.database,
            )

            if self._spectrum_plot_folder is not None:
                folder_path = os.path.join(
                    get_figures_path(), self._spectrum_plot_folder
                )
                file_name = f"{fitter.result.measured_spectrum.meas_name.split('/')[-1].split('.')[0]}.svg"

                if not os.path.exists(folder_path):
                    os.makedirs(folder_path)

                plt = fitter.result.generate_plot()
                plt.savefig(os.path.join(folder_path, file_name))
                plt.close()

            self._results.append(fitter.result)
            self._concentrations.append(fitter.concentration)

            meas_name = fitter.result.measured_spectrum.meas_name
            self._positions.append(
                tuple(
                    float(pos)
                    for pos in re.findall(r"\-X([0-9\.]+)\-Y([0-9\.]+)$", meas_name)[0]
                )
            )

    def generate_concentration_heatmap(self) -> "plt":
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

        if conc.shape[1] < conc.shape[0] / 5:
            plt.figure(figsize=(4, 12), dpi=100)
            tight = {"pad": 0.1, "rect": (0.02, 0.02, 0.96, 0.99)}
        elif conc.shape[0] < conc.shape[1] / 5:
            plt.figure(figsize=(12, 4), dpi=100)
            tight = {"pad": 0.1, "rect": (0.02, 0.02, 1.08, 0.99)}
        else:
            tight = {"pad": 0.1, "rect": (0.02, 0.02, 0.99, 0.99)}

        df = pd.DataFrame(conc)
        ax = sns.heatmap(df, cmap="crest")
        ax.invert_yaxis()
        ax.set_title("Concentration as a function of position")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_xlim(non_nan_cols[0], non_nan_cols[-1] + 1)
        ax.set_ylim(non_nan_rows[0], non_nan_rows[-1] + 1)
        plt.tight_layout(**tight)

        return plt

    def show_concentration_heatmap(self) -> None:
        """
        Show the heatmap of the concentration as a function of position.
        """
        return self.generate_concentration_heatmap().show()

    def generate_concentration_plot(
        self, x: "Optional[int]" = None, y: "Optional[int]" = None
    ) -> "plt":
        """
        Generate a heatmap of the concentration as a function of position.
        """
        import numpy as np

        rulers = [x, y]
        if all(r is None for r in rulers):
            non_nan_rows = np.where(~np.isnan(self.concentration_map.T).all(axis=1))[0]
            if not any(non_nan_rows):
                raise ValueError("No non-NaN values in the concentration map.")
            y = non_nan_rows[0]
        if all(r is not None for r in rulers):
            raise ValueError("Only one of x or y can be specified.")

        from matplotlib import pyplot as plt

        if y is not None:
            conc = self.concentration_map[:, y]
            ruler = "y"
            other = "x"
        else:
            conc = self.concentration_map[x, :]
            ruler = "x"
            other = "y"

        positions = np.arange(len(conc))
        non_nan_indices = np.where(~np.isnan(conc))[0]

        plt.figure(figsize=(9, 6), dpi=100)
        plt.scatter(positions, conc, c="b")
        plt.xlim(non_nan_indices[0] - 1, non_nan_indices[-1] + 1)
        plt.xlabel("Position")
        plt.ylabel("Concentration [VMR]")
        plt.title(
            f"Concentration as a function of the {other} position "
            + f"for ${ruler} = {x if y is None else y}$"
        )
        tight = {"pad": 0.1, "rect": (0.02, 0.02, 0.98, 0.97)}
        plt.tight_layout(**tight)
        return plt

    def show_concentration_plot(
        self, x: "Optional[int]" = None, y: "Optional[int]" = None
    ) -> None:
        """
        Show the heatmap of the concentration as a function of position.
        """
        return self.generate_concentration_plot(x=x, y=y).show()
