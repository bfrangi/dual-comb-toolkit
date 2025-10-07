from datetime import datetime

from lib.files import get_figures_path, get_measurement_names, save_mapping_report
from lib.plots import use_latex
from lib.shortcuts import map_measurement_concentration

####################################################################################################
# Measurement parameters                                                                           #
####################################################################################################

# Molecule and physical conditions.

molecule = "CH4"
pressure = 101325  # Pa
temperature = 2100  # K
length = 0.07  # m

# Simulation range.

wl_min = 3426.8  # nm
wl_max = 3428.1  # nm

# Measurement names.

mapping_name = "7a"
measurement_names = get_measurement_names(mapping_name)
baseline_names = []

# Radio frequency comb specifications.

center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
acq_freq = 400000.0  # Hz

# Optical comb specifications.

number_of_teeth = 12
laser_wavelength = (wl_max + wl_min) / 2 * 1e-9  # m
optical_comb_spacing = 1250e6  # Hz


####################################################################################################
# Fitting parameters                                                                               #
####################################################################################################

# Fitter, database, initial guess and allowed concentration bounds.

fitter = "normal_gpu"
database = "hitemp"
initial_guess = 0.0001  # VMR
lower_bound = 0.0  # VMR
upper_bound = 1.0  # VMR

# Noise filtering.

tooth_std_threshold = 10  # Teeth with a standard deviation above `tooth_std_threshold * mean_std` will be discarded.
sub_measurements = (
    10  # Number of sub-measurements used to obtain the standard deviation of the teeth.
)

# Output and plotting parameters.

verbose = True
spectrum_plot_folder = mapping_name

# Use LaTeX for plotting.

use_latex()


####################################################################################################
# Mapping                                                                                          #
####################################################################################################

mapper = map_measurement_concentration(
    meas_names=measurement_names,
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
    database=database,
    initial_guess=initial_guess,
    lower_bound=lower_bound,
    upper_bound=upper_bound,
    verbose=verbose,
    spectrum_plot_folder=spectrum_plot_folder,
    tooth_std_threshold=tooth_std_threshold,
    sub_measurements=sub_measurements,
)


####################################################################################################
# Results plots                                                                                    #
####################################################################################################

if spectrum_plot_folder:
    folder_path = f"{get_figures_path()}{spectrum_plot_folder}"
    heatmap_path = f"{folder_path}/heatmap.pdf"
    plot_path = f"{folder_path}/plot.pdf"

    plt = mapper.generate_concentration_heatmap()
    plt.savefig(heatmap_path)
    plt.close()

    plt = mapper.generate_concentration_plot(x=1)
    plt.savefig(plot_path)
    plt.close()

    if verbose:
        print(f"Concentration heatmap saved to {heatmap_path}.")
        print(f"Concentration plot saved to {plot_path}.")
else:
    mapper.show_concentration_heatmap()
    mapper.show_concentration_plot(x=1)


####################################################################################################
# Mapping report                                                                                   #
####################################################################################################

report_filename = (
    f"{mapping_name.split('/')[-1]} @ {datetime.now().strftime('%Y-%m-%d %H-%M-%S')}"
)

save_mapping_report(
    filename=report_filename,
    data={
        "molecule": molecule,
        "pressure": pressure,
        "temperature": temperature,
        "path_length": length,
        "optical_comb_spacing": optical_comb_spacing / 1e9,  # Convert to GHz
        "number_of_teeth": number_of_teeth,
        "laser_wavelength": laser_wavelength * 1e9,  # Convert to nm
        "wl_min": wl_min,
        "wl_max": wl_max,
        "comb_spacing": freq_spacing,
        "center_frequency": center_freq,
        "acquisition_frequency": acq_freq,
        "fitter": fitter,
        "initial_guess": initial_guess,
        "lower_bound": lower_bound,
        "upper_bound": upper_bound,
        "sub_measurements": sub_measurements,
        "tooth_std_threshold": tooth_std_threshold,
        "results": mapper.results,
    },
    verbose=verbose,
)
