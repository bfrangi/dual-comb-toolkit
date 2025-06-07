from matplotlib import pyplot as plt

from lib.plots import tight
from lib.shortcuts import get_measurement_transmission, simulate_line

# Define the molecule, VMR, pressure, temperature, and length.
molecule = "CH4"
vmr = 0.005  # volume mixing ratio
pressure = 101325  # Pa
temperature = 298  # K
length = 0.07  # m
wl_min = 3427.0  # nm
wl_max = 3427.9  # nm

# Specify the name and specifications of the measurement.
meas_name = "4a/Position-X10-Y1"
baseline_names = []
center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
number_of_teeth = 12
laser_wavelength = 3427.54e-9  # m
optical_comb_spacing = 1250e6  # Hz
acq_freq = 400000.0  # Hz

# Simulate the transmission spectrum.
x_sim, y_sim = simulate_line(
    molecule=molecule,
    wl_min=wl_min,
    wl_max=wl_max,
    vmr=vmr,
    pressure=pressure,
    temperature=temperature,
    length=length,
)

# Get the measured transmission spectrum.
x_meas, y_meas = get_measurement_transmission(
    meas_name=meas_name,
    center_freq=center_freq,
    freq_spacing=freq_spacing,
    number_of_teeth=number_of_teeth,
    laser_wavelength=laser_wavelength,
    optical_comb_spacing=optical_comb_spacing,
    acq_freq=acq_freq,
    baseline_names=baseline_names,
)


# Plot the simulated and measured transmission spectra.

plt.plot(x_sim, y_sim, label="Simulated", color="blue", zorder=0)
plt.scatter(x_meas, y_meas, label="Measured", color="red", zorder=1)
plt.legend()
plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmission")
plt.title(
    f"Transmission spectrum of {molecule} at {pressure:.2f} Pa "
    + f"and {temperature:.2f} K.\n{length:.3f} m path length, {vmr:.3f} VMR."
)
plt.tight_layout(**tight)
plt.show()
