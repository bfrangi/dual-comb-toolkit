from matplotlib import pyplot as plt

from lib.plots import tight
from lib.shortcuts import fit_measurement_concentration

# Define the molecule, VMR, pressure, temperature, and length.
molecule = 'CH4'
pressure = 53328.94736842  # Pa
temperature = 298  # K
length = 0.055  # m
wl_min = 3427.1  # nm
wl_max = 3427.8  # nm

# Specify the name and specifications of the measurement.
meas_name = 'cell-sweep-10-36-17-03-2025/Position-X7-Y1'
center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
number_of_teeth = 38
laser_wavelength = 3427.437e-9  # m
high_freq_modulation = 500e6  # Hz
acq_freq = 400000.0  # Hz

# Fit the concentration.

vmr, x_sim, y_sim, x_meas, y_meas = fit_measurement_concentration(
    meas_name,
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
    initial_guess=0.5
)

# Plot the simulated and measured transmission spectra.

plt.plot(x_sim, y_sim, label='Simulated', color='blue', zorder=0)
plt.scatter(x_meas, y_meas, label='Measured', color='red', zorder=1)
plt.legend()
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission')
plt.title(f'Transmission spectrum of {molecule} at {pressure:.2f} Pa ' +
          f'and {temperature:.2f} K.\n{length:.3f} m path length, {vmr:.3f} VMR.')
plt.tight_layout(**tight)
plt.show()

print(f'The concentration of {molecule} is {vmr:.3f} VMR.')