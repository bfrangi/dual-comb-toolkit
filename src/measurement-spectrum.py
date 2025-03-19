from matplotlib import pyplot as plt

from lib.plots import tight
from lib.shortcuts import get_measurement_spectrum

# Define the molecule, VMR, pressure, temperature, and length.
meas_name = 'cell-sweep-1/Position-X1-Y1-sample.lvm'
acq_freq = 400000.0  # Hz

wl, transmission = get_measurement_spectrum(meas_name, acq_freq=acq_freq)

# Plot the measured transmission spectrum.
plt.plot(wl, transmission, color='red')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission')
plt.title(f'Measured spectrum of {meas_name}.')
plt.tight_layout(**tight)
plt.show()