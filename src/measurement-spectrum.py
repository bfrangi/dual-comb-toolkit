from matplotlib import pyplot as plt

from lib.plots import article_tight, use_latex
from lib.shortcuts import get_measurement_spectrum

# Define the molecule, VMR, pressure, temperature, and length.
measurement_name = "cell-sweep-10-34-17-03-2025/Position-X1-Y1-sample.lvm"
acq_freq = 400000.0  # Hz
latex = True # Whether to use LaTeX for plotting.

if latex:
    use_latex()

wl, transmission = get_measurement_spectrum(measurement_name, acq_freq=acq_freq)

# Plot the measured transmission spectrum.
plt.plot(wl, transmission, color="red")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmission [-]")
plt.title(f'Measured spectrum of\n{measurement_name}.')
plt.tight_layout(**article_tight)
plt.show()
