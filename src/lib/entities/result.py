from typing import TYPE_CHECKING

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
