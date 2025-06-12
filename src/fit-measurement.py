from matplotlib import pyplot as plt

from lib.plots import tight
from lib.shortcuts import fit_measurement_concentration

# Define the molecule, VMR, pressure, temperature, and length.
molecule = "CH4"
pressure = 101325  # Pa
temperature = 298  # K
length = 0.07  # m
wl_min = 3426.8  # nm
wl_max = 3428.1  # nm

# Specify the name and specifications of the measurement.
meas_name = "1a/Position-X1-Y10"
baseline_names = []
center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
number_of_teeth = 12
laser_wavelength = (wl_max + wl_min) / 2 * 1e-9  # m
optical_comb_spacing = 1250e6  # Hz
fitter = "normal_gpu"
initial_guess = 0.0001  # VMR
lower_bound = 0.0  # VMR
upper_bound = 0.1  # VMR
verbose = True
acq_freq = 400000.0  # Hz

# Fit the concentration.

vmr, x_sim, y_sim, x_meas, y_meas = fit_measurement_concentration(
    meas_name,
    center_freq=center_freq,
    freq_spacing=freq_spacing,
    number_of_teeth=number_of_teeth,
    laser_wavelength=laser_wavelength,
    optical_comb_spacing=optical_comb_spacing,
    acq_freq=acq_freq,
    wl_min=wl_min,
    wl_max=wl_max,
    molecule=molecule,
    pressure=pressure,
    temperature=temperature,
    length=length,
    initial_guess=initial_guess,
    fitter=fitter,
    baseline_names=baseline_names,
    lower_bound=lower_bound,
    upper_bound=upper_bound,
    verbose=verbose,
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
