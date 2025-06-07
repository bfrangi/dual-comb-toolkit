from typing import TYPE_CHECKING
from matplotlib import pyplot as plt
from lib.plots import tight

if TYPE_CHECKING:
    from lib.entities.spectrum import MeasuredSpectrum, SimulatedSpectrum


class Result:
    def __init__(self, measured_spectrum: 'MeasuredSpectrum',
                 simulated_spectrum: 'SimulatedSpectrum') -> None:
        self.measured_spectrum = measured_spectrum
        self.simulated_spectrum = simulated_spectrum

    @property
    def concentration(self) -> float:
        """
        Get the concentration of the molecule.
        """
        return self.simulated_spectrum.concentration

    vmr = concentration

    @property
    def molecule(self) -> str:
        """
        Get the name of the molecule.
        """
        return self.simulated_spectrum.molecule

    @property
    def pressure(self) -> float:
        """
        Get the pressure in Pa.
        """
        return self.simulated_spectrum.pressure

    @property
    def temperature(self) -> float:
        """
        Get the temperature in K.
        """
        return self.simulated_spectrum.temperature

    @property
    def length(self) -> float:
        """
        Get the length of the absorption path in m.
        """
        return self.simulated_spectrum.length

    def generate_plot(self) -> 'plt':
        """
        Generate a plot of the measured and simulated transmission spectra.
        
        Returns
        -------
        plt : matplotlib.pyplot
            The plot object containing the measured and simulated spectra.
        """
        plt.plot(
            self.simulated_spectrum.x_nm,
            self.simulated_spectrum.y_nm,
            label=f"Simulation for {self.concentration:.3f} VMR",
            color="blue",
            zorder=0,
        )
        plt.scatter(
            self.measured_spectrum.x_nm,
            self.measured_spectrum.y_nm,
            label="Measurement",
            color="red",
            zorder=1,
        )
        plt.legend()
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Transmission")
        plt.title(
            f"Measured and Simulated Transmission spectrum of {self.molecule}\nat {self.pressure:.2f} Pa "
            + f"and {self.temperature:.2f} K. {self.length:.3f} m path length, {self.concentration:.3f} VMR "
        )
        plt.tight_layout(**tight)
        return plt
    
    def show_plot(self) -> None:
        """
        Show the plot of the measured and simulated transmission spectra.
        """
        return self.generate_plot().show()

