import os

from matplotlib import pyplot as plt

from lib.entities import Result, SimulatedSpectrum
from lib.files import get_figures_path, initialize_figures_folder
from lib.plots import use_latex
from lib.shortcuts import get_measurement, simulate_line

####################################################################################################
# Simulation parameters                                                                            #
####################################################################################################

# Molecule and physical conditions.

molecule = "CH4"
vmr = 0.00065  # VMR
pressure = 101325  # Pa
temperature = 298  # K
length = 0.07  # m

# Simulation range.

wl_min = 3426.8  # nm
wl_max = 3428.1  # nm


####################################################################################################
# Measurement parameters                                                                           #
####################################################################################################

# Measurement name.

measurement_name = "11a/Position-X1-Y1"
baseline_names = []

# Radio frequency comb specifications.

center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
acq_freq = 400000.0  # Hz

# Optical comb specifications.

number_of_teeth = 12
laser_wavelength = (wl_max + wl_min) / 2 * 1e-9 + 0.11e-9  # m
optical_comb_spacing = 1250e6  # Hz


####################################################################################################
# Processing parameters                                                                            #
####################################################################################################

# Noise filtering.

tooth_std_threshold = 0.4 # Teeth with a standard deviation above this threshold will be discarded.
sub_measurements = 10 # Number of sub-measurements used to obtain the standard deviation of the teeth.

# Measurement scaling factor.

scaling_factor = 1.001  # Adjust the measured spectrum by this factor.

# Plotting parameters.

spectrum_plot_folder = 'process-measurement-output'
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
)

measured_spectrum = measurement.transmission_spectrum
measured_spectrum.scale_by(scaling_factor)


####################################################################################################
# Plots                                                                                            #
####################################################################################################


result = Result(measured_spectrum=measured_spectrum, simulated_spectrum=simulated_spectrum)

result.generate_plot()

if spectrum_plot_folder:
    initialize_figures_folder(spectrum_plot_folder)

    file_name = measurement_name.split("/")[-1] + ".svg"
    folder_path = os.path.join(get_figures_path(), spectrum_plot_folder)
    file_path = os.path.join(folder_path, file_name)

    plt.savefig(file_path)

    print(f"Plot saved to {file_path}.")

plt.show()