from lib.files import get_measurement_names
from lib.shortcuts import map_measurement_concentration
from lib.plots import use_latex

# Define the molecule, pressure, temperature, and length.

molecule = "CH4"
pressure = 101325  # Pa
temperature = 298  # K
length = 0.07  # m
wl_min = 3427.0  # nm
wl_max = 3427.9  # nm

# Name and specifications of the measurement.

sweep_name = "1a"
meas_names = get_measurement_names(sweep_name)
baseline_names = []
center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
number_of_teeth = 12
laser_wavelength = 3427.437e-9  # m
optical_comb_spacing = 1250e6  # Hz
acq_freq = 400000.0  # Hz
fitter = "normal"
initial_guess = 0.001  # VMR
spectrum_plot_folder = sweep_name

# Use LaTeX for plotting.

use_latex()

# Map the measurements.

mapper = map_measurement_concentration(
    meas_names=meas_names,
    center_freq=center_freq,
    freq_spacing=freq_spacing,
    number_of_teeth=number_of_teeth,
    laser_wavelength=laser_wavelength,
    optical_comb_spacing=optical_comb_spacing,
    acq_freq=acq_freq,
    molecule=molecule,
    pressure=pressure,
    temperature=temperature,
    length=length,
    wl_min=wl_min,
    wl_max=wl_max,
    baseline_names=baseline_names,
    fitter=fitter,
    initial_guess=initial_guess,
    spectrum_plot_folder=spectrum_plot_folder,
)

mapper.show_concentration_heatmap()

mapper.show_concentration_plot(x=1)
