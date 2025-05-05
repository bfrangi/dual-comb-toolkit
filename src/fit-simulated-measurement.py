import csv
import os
import time
import warnings

from matplotlib import pyplot as plt
from numpy import array

from lib.benchmarking import Time
from lib.constants import c
from lib.conversions import delta_frequency_to_delta_wavelength
from lib.files import get_figures_path, get_reports_path
from lib.plots import tight
from lib.shortcuts import fit_simulated_measurement_concentration

warnings.filterwarnings("ignore", category=FutureWarning)

# Define the parameters for simulating the measurement.

molecule = "CH4"
database = "hitran"

vmr = 0.01  # volume mixing ratio
pressure = 53328.94736842  # Pa
temperature = 298  # K
length = 0.055  # m
laser_wavelength = 3427.41  # nm

wl_min = 3427.0  # nm
wl_max = 3427.9  # nm
min_comb_span = 0.15  # nm

std_dev = 0.014  # unitless
number_of_teeth_for_std_dev = 30  # teeth
x_shift_std_dev = 0.02  # nm
scaling_std_dev = 1  # unitless


n_simulations_per_config = 100

comb_spacings = [
    2.6, 2.7, 2.8, 2.9, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    2.8, 2.9, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4,
    2.5, 2.6, 2.7, 2.8, 2.9, 1.3, 1.4, 1.5, 1.6, 1.7,
    1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 1.2,
    1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2,
    2.3, 2.4, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.1, 2.2, 2.3, 2.4, 1.0, 1.1, 1.2, 1.3, 1.4,
    1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 1.0, 1.1,
    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 0.6, 0.7, 0.8,
    0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 0.6, 0.7,
    0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 0.6, 0.7,
    0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 0.5, 0.6, 0.7,
    0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 0.5, 0.6, 0.7, 0.8,
    0.9, 1.0, 1.1, 1.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.1, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 0.4,
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.5, 0.6, 0.7,
    0.8, 0.9, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.4, 0.5, 0.6, 0.7, 0.8, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.3, 0.4, 0.5, 0.6, 0.7, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.3, 0.4, 0.5, 0.6, 0.7, 0.3, 0.4,
    0.5, 0.6
]  # GHz
comb_spacings = [i * 1e9 for i in comb_spacings]  # Hz
numbers_of_teeth = [
    5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
    6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 10, 11, 11, 11, 11, 11,
    11, 11, 11, 11, 11, 11, 11, 11, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 14, 14, 14,
    14, 14, 14, 14, 14, 14, 14, 14, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 17, 17, 17,
    17, 17, 17, 17, 17, 17, 18, 18, 18, 18,
    18, 18, 18, 18, 19, 19, 19, 19, 19, 19,
    19, 20, 20, 20, 20, 20, 20, 20, 20, 21,
    21, 21, 21, 21, 21, 21, 22, 22, 22, 22,
    22, 22, 23, 23, 23, 23, 23, 23, 24, 24,
    24, 24, 24, 25, 25, 25, 25, 25, 26, 26,
    26, 26, 26, 27, 27, 27, 27, 27, 28, 28,
    28, 28, 28, 29, 29, 29, 29, 29, 30, 30,
    30, 30
]  # teeth

# Define the fitting parameters.

normalize = True
initial_guess = 0.001
fitter = "normal_gpu"

# Info

print(f"Simulating {min(len(numbers_of_teeth), len(comb_spacings))} configurations.")

# Simulate for every combination of number of teeth and comb spacing.

results = []

for n_teeth, spacing in zip(numbers_of_teeth, comb_spacings):
    # Make sure the comb span is within the limits
    comb_span = delta_frequency_to_delta_wavelength(
        spacing * (n_teeth - 1), c / laser_wavelength * 1e9
    )
    comb_spacing = comb_span / (n_teeth - 1)  # nm

    # Calculate the laser wavelength slack range
    laser_wavelength_slack = (-comb_spacing, comb_spacing)  # nm

    if comb_span + 2 * comb_spacing > wl_max - wl_min or comb_span < min_comb_span:
        continue

    # Print the simulation parameters

    print(
        f"Number of teeth: {n_teeth}, High frequency "
        + f"modulation: {spacing / 1e9:.2f} GHz"
    )

    # Create a folder in `figures` with name like `8 x 1.0 GHz` to store the plots.

    folder_name = f"{n_teeth} x {spacing / 1e9:.1f} GHz"
    folder_path = os.path.join(get_figures_path(), folder_name)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Simulate the measurement and fit the concentration `n_simulations_per_config` times

    config_results = []

    for i in range(n_simulations_per_config):
        indent = " " * len(f"{i + 1}. ")

        # Simulate the transmission spectrum

        with Time(
            f"{i + 1}. Simulating for {n_teeth} teeth and {spacing / 1e9:.2f} GHz"
        ):
            x_meas, y_meas, f = fit_simulated_measurement_concentration(
                molecule=molecule,
                wl_min=wl_min,
                wl_max=wl_max,
                vmr=vmr,
                pressure=pressure,
                temperature=temperature,
                length=length,
                laser_wavelength=laser_wavelength,
                high_freq_modulation=spacing,
                number_of_teeth=n_teeth,
                database=database,
                std_dev=std_dev,
                number_of_teeth_for_std_dev=number_of_teeth_for_std_dev,
                x_shift_std_dev=x_shift_std_dev,
                scaling_std_dev=scaling_std_dev,
                laser_wavelength_slack=laser_wavelength_slack,
                normalize=normalize,
                initial_guess=initial_guess,
                fitter=fitter,
            )

            meas = f.measured_transmission
            sim = f.simulated_transmission

            x_sim, y_sim, x_meas_fitted, y_meas_fitted = (
                sim.x_nm,
                sim.y_nm,
                meas.x_nm,
                meas.y_nm,
            )

        config_results.append(f.concentration)
        print(f"{indent}The fitted concentration is {f.concentration:.6f} VMR.")

        # Plot the simulated and measured transmission spectra

        with Time(f"{indent}Plotting"):
            plt.plot(
                x_sim,
                y_sim,
                label=f"Simulation for {f.concentration:.3f} VMR",
                color="blue",
                zorder=0,
            )
            plt.scatter(
                x_meas,
                y_meas,
                label="Simulated Measurement (pre-fitting)",
                color="red",
                zorder=1,
            )
            plt.scatter(
                x_meas_fitted,
                y_meas_fitted,
                label="Simulated Measurement (fitted)",
                color="green",
                zorder=2,
            )
            plt.legend()
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Transmission")
            plt.title(
                f"Transmission spectrum of {molecule} at {pressure:.2f} Pa "
                + f"and {temperature:.2f} K.\n{length:.3f} m path length, {vmr:.3f} VMR "
                + f"(fitted as {f.concentration:.3f} VMR)."
            )
            plt.tight_layout(**tight)
            plt.savefig(f"{folder_path}/fit-simulated-measurement-{i}.svg")
            plt.clf()

    config_results = array(config_results)
    result = (n_teeth, spacing, config_results.mean(), config_results.std())

    results.append(result)

    # Print the result
    print(
        f"Number of teeth: {n_teeth}, Comb spacing: {spacing / 1e9:.2f} GHz, "
        + f"Mean concentration: {result[2]:.6f} VMR, Standard deviation: {result[3]:.6f} VMR"
    )
    print("--------------------")


timestr = time.strftime("%Y%m%d-%H%M%S")
with open(f"{get_reports_path()}report-{timestr}.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter=",")
    writer.writerow(
        (
            "Number of teeth",
            "Comb spacing (Hz)",
            "Mean concentration (VMR)",
            "Standard deviation (VMR)",
        )
    )
    writer.writerows(results)
