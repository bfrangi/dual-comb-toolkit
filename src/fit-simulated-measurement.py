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

molecule = 'CH4'
database = 'hitran'

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
laser_wavelength_std_dev = 0.01  # nm

n_simulations_per_config = 10

comb_spacings = [(i+1)*100e6 for i in range(30)] * 26 # Hz
numbers_of_teeth = [i for i in range(5, 31) for _ in range(30)] # teeth

# Define the fitting parameters.

normalize = True
initial_guess = 0.001
fitter = 'normal_gpu'

# Info

print(f'Simulating {min(len(numbers_of_teeth), len(comb_spacings))} configurations.')

# Simulate for every combination of number of teeth and comb spacing.

results = []

for n_teeth, spacing in zip(numbers_of_teeth, comb_spacings):
    # Make sure the comb span is within the limits.

    comb_span = delta_frequency_to_delta_wavelength(
        spacing*(n_teeth - 1), c / laser_wavelength * 1e9)

    if comb_span > wl_max - wl_min or comb_span < min_comb_span:
        continue

    # Print the simulation parameters

    print(f'Number of teeth: {n_teeth}, High frequency ' +
          f'modulation: {spacing/1e9:.2f} GHz')

    # Create a folder in `figures` with name like `8 x 1.0 GHz` to store the plots.

    folder_name = f'{n_teeth} x {spacing / 1e9:.1f} GHz'
    folder_path = os.path.join(get_figures_path(), folder_name)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Simulate the measurement and fit the concentration `n_simulations_per_config` times.

    config_results = []

    for i in range(n_simulations_per_config):
        indent = ' '*len(f'{i+1}. ')

        # Simulate the transmission spectrum.

        with Time(f'{i+1}. Simulating for {n_teeth} teeth and {spacing/1e9:.2f} GHz'):
            x_meas, y_meas, f = fit_simulated_measurement_concentration(
                molecule=molecule, wl_min=wl_min, wl_max=wl_max, vmr=vmr, pressure=pressure,
                temperature=temperature, length=length, laser_wavelength=laser_wavelength,
                high_freq_modulation=spacing, number_of_teeth=n_teeth, database=database,
                std_dev=std_dev, number_of_teeth_for_std_dev=number_of_teeth_for_std_dev,
                x_shift_std_dev=x_shift_std_dev, scaling_std_dev=scaling_std_dev,
                laser_wavelength_std_dev=laser_wavelength_std_dev, normalize=normalize,
                initial_guess=initial_guess, fitter=fitter)

            meas = f.measured_transmission
            sim = f.simulated_transmission

            x_sim, y_sim, x_meas_fitted, y_meas_fitted = sim.x_nm, sim.y_nm, meas.x_nm, meas.y_nm

        config_results.append(f.concentration)
        print(f'{indent}The fitted concentration is {f.concentration:.6f} VMR.')

        # Plot the simulated and measured transmission spectra.

        with Time(f'{indent}Plotting'):
            plt.plot(x_sim, y_sim,
                     label=f'Simulation for {f.concentration:.3f} VMR', color='blue', zorder=0)
            plt.scatter(x_meas, y_meas, label='Simulated Measurement (pre-fitting)',
                        color='red', zorder=1)
            plt.scatter(x_meas_fitted, y_meas_fitted,
                        label='Simulated Measurement (fitted)', color='green', zorder=2)
            plt.legend()
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Transmission')
            plt.title(f'Transmission spectrum of {molecule} at {pressure:.2f} Pa ' +
                      f'and {temperature:.2f} K.\n{length:.3f} m path length, {vmr:.3f} VMR ' +
                      f'(fitted as {f.concentration:.3f} VMR).')
            plt.tight_layout(**tight)
            plt.savefig(f'{folder_path}/fit-simulated-measurement-{i}.svg')
            plt.clf()

    config_results = array(config_results)
    result = (n_teeth, spacing, config_results.mean(), config_results.std())

    results.append(result)

    # Print the result
    print(f'Number of teeth: {n_teeth}, Comb spacing: {spacing/1e9:.2f} GHz, ' +
          f'Mean concentration: {result[2]:.6f} VMR, Standard deviation: {result[3]:.6f} VMR')
    print('--------------------')


timestr = time.strftime("%Y%m%d-%H%M%S")
with open(f'{get_reports_path()}report-{timestr}.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(('Number of teeth', 'Comb spacing (Hz)',
                    'Mean concentration (VMR)', 'Standard deviation (VMR)'))
    writer.writerows(results)
