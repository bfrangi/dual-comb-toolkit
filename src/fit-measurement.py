import os

from lib.files import get_figures_path, initialize_figures_folder
from lib.plots import use_latex
from lib.shortcuts import fit_measurement_concentration

####################################################################################################
# Measurement parameters                                                                           #
####################################################################################################

# Molecule and physical conditions.

molecule = "CH4"
pressure = 101325  # Pa
temperature = 298  # K
length = 0.07  # m

# Simulation range.

wl_min = 3426.8  # nm
wl_max = 3428.1  # nm

# Measurement name.

measurement_name = "2a/Position-X1-Y19"
baseline_names = []

# Radio frequency comb specifications.

center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
acq_freq = 400000.0  # Hz

# Optical comb specifications.

number_of_teeth = 26
laser_wavelength = (wl_max + wl_min) / 2 * 1e-9  # m
optical_comb_spacing = 700e6  # Hz


####################################################################################################
# Fitting parameters                                                                               #
####################################################################################################

# Fitter, initial guess and allowed concentration bounds.

fitter = "normal_gpu"
initial_guess = 0.0001  # VMR
lower_bound = 0.0  # VMR
upper_bound = 0.15  # VMR

# Noise filtering.

tooth_std_threshold = 1.5  # Teeth with a standard deviation above `tooth_std_threshold * mean_std` will be discarded.
sub_measurements = (
    10  # Number of sub-measurements used to obtain the standard deviation of the teeth.
)

# Output and plotting parameters.

verbose = True
spectrum_plot_folder = "fit-measurement-output"

# Use LaTeX for plotting.

use_latex()


####################################################################################################
# Fitting                                                                                          #
####################################################################################################

f = fit_measurement_concentration(
    measurement_name,
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
    sub_measurements=sub_measurements,
    tooth_std_threshold=tooth_std_threshold,
    verbose=verbose,
)


####################################################################################################
# Results plots                                                                                    #
####################################################################################################

plt = f.result.generate_plot()

if spectrum_plot_folder:
    initialize_figures_folder(spectrum_plot_folder)

    file_name = measurement_name.split("/")[-1] + ".svg"
    folder_path = os.path.join(get_figures_path(), spectrum_plot_folder)
    file_path = os.path.join(folder_path, file_name)

    plt.savefig(file_path)

    print(f"Plot saved to {file_path}.")

plt.show()
