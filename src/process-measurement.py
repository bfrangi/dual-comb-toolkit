import os

from matplotlib import pyplot as plt

from lib.entities import Result, SimulatedSpectrum
from lib.files import get_figures_path, initialize_figures_folder
from lib.plots import use_latex
from lib.shortcuts import get_measurement, simulate_line

####################################################################################################
#  Simulation parameters                                                                           #
####################################################################################################

# Molecule and database.

molecule = "CH4"
"""Molecule to simulate."""
database = "hitemp"
"""Database to use for the simulation. Can be 'hitran', 'hitemp', 'exomol' or 'geisa'. Some may not
be available for all molecules."""

# Physical conditions.

vmr = 0.02616  # VMR
"""Volume mixing ratio of the molecule in the gas mixture."""
pressure = 101325  # Pa
"""Pressure of the gas mixture."""
temperature = 1277.46  # K
"""Temperature of the gas mixture."""
length = 0.07  # m
"""Path length of the gas mixture."""

# Simulation range.

wl_min = 3427.0  # nm
"""Minimum wavelength of the simulation range."""
wl_max = 3427.9  # nm
"""Maximum wavelength of the simulation range."""

####################################################################################################
# Measurement parameters                                                                           #
####################################################################################################

# Measurement name.

measurement_name = "7b/Position-X1-Y1"
"""Name of the measurement to process."""
baseline_names = []
"""List of baseline measurement names to use for processing."""

# Radio frequency comb specifications.

center_freq = 40000.0  # Hz
"""Center frequency of the radio frequency comb."""
freq_spacing = 200.0  # Hz
"""Frequency spacing of the radio frequency comb."""
acq_freq = 400000.0  # Hz
"""Acquisition frequency of the radio frequency comb."""
flip = False
"""If True, the measured transmission spectrum will be flipped with respect to the center frequency."""

# Optical comb specifications.

number_of_teeth = 12
"""Number of teeth in the optical frequency comb."""
optical_comb_spacing = 1250e6  # Hz
"""Optical frequency comb spacing."""
laser_wavelength = 3427.41e-9 + 0.07e-9  # m
"""Wavelength of the laser used to probe the gas mixture."""


####################################################################################################
# Processing parameters                                                                            #
####################################################################################################

# Noise filtering.

tooth_std_threshold = 10
"""Threshold for the standard deviation of the comb teeth. If the standard deviation of a tooth is 
larger than `tooth_std_threshold` times the mean standard deviation of all teeth, the tooth will be 
discarded. Note this could give unexpected results if combined with `remove_teeth_indices`."""
sub_measurements = 10
"""Number of sub-measurements to use for calculating the standard deviation of the teeth."""

# Removing noisy teeth
remove_teeth_indices = [11, 12]
"""List of tooth indices to be removed from the fitting."""

# Measurement scaling factor.

scaling_factor = 0.993  # Adjust the measured spectrum by this factor.
"""Scaling factor to apply to the measured transmission spectrum."""

# Plotting parameters.

spectrum_plot_folder = "process-measurement-output"
"""Folder to save the spectrum plots. If None, plots will not be saved."""

# Use LaTeX for plotting.

latex = False
"""If True, use LaTeX for plotting."""

if latex:
    use_latex()


####################################################################################################
# Simulation                                                                                       #
####################################################################################################

# Simulate the transmission spectrum.

x_sim, y_sim = simulate_line(
    molecule=molecule,
    wl_min=wl_min,
    wl_max=wl_max,
    vmr=vmr,
    pressure=pressure,
    temperature=temperature,
    length=length,
    database=database,
)

simulated_spectrum = SimulatedSpectrum(
    x=x_sim,
    y=y_sim,
    xu="nm",
    molecule=molecule,
    pressure=pressure,
    temperature=temperature,
    length=length,
    concentration=vmr,
)


####################################################################################################
# Measurement                                                                                      #
####################################################################################################

# Get the measured transmission spectrum.

measurement = get_measurement(
    meas_name=measurement_name,
    center_freq=center_freq,
    freq_spacing=freq_spacing,
    number_of_teeth=number_of_teeth,
    laser_wavelength=laser_wavelength,
    optical_comb_spacing=optical_comb_spacing,
    acq_freq=acq_freq,
    baseline_names=baseline_names,
    sub_measurements=sub_measurements,
    tooth_std_threshold=tooth_std_threshold,
    flip=flip,
)
measurement._compute_transmission()
measurement.remove_teeth(remove_teeth_indices)

measured_spectrum = measurement.transmission_spectrum
measured_spectrum.scale_by(scaling_factor)


####################################################################################################
# Plots                                                                                            #
####################################################################################################


result = Result(
    measured_spectrum=measured_spectrum, simulated_spectrum=simulated_spectrum
)

result.generate_plot()

if spectrum_plot_folder:
    initialize_figures_folder(spectrum_plot_folder)

    file_name = measurement_name.split("/")[-1] + ".svg"
    folder_path = os.path.join(get_figures_path(), spectrum_plot_folder)
    file_path = os.path.join(folder_path, file_name)

    plt.savefig(file_path)

    print(f"Plot saved to {file_path}.")

plt.show()
