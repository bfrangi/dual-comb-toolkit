import csv
import os
import time

from matplotlib import pyplot as plt
from numpy import array

from lib.constants import c
from lib.conversions import delta_frequency_to_delta_wavelength
from lib.files import get_figures_path, get_reports_path
from lib.plots import tight
from lib.shortcuts import fit_simulated_measurement_concentration

# Define the parameters for simulating the measurement.
molecule = 'CH4'
wl_min = 3427.0  # nm
wl_max = 3427.9  # nm
vmr = 0.01  # volume mixing ratio
pressure = 53328.94736842  # Pa
temperature = 298  # K
length = 0.055  # m
laser_wavelength = 3427.41  # nm
database = 'hitran'
std_dev = 0.014  # unitless
x_shift_std_dev = 0.02  # nm
scaling_std_dev = 1  # unitless
laser_wavelength_std_dev = 0.01  # nm
normalize = True
min_comb_span = 0.15  # nm
high_freq_modulations = [(i+1)*100e6 for i in range(30)]  # Hz
numbers_of_teeth = range(5,31) # teeth

# Define the fitting parameters.
initial_guess = 0.5

results = []

for number_of_teeth in numbers_of_teeth:
    for high_freq_modulation in high_freq_modulations:
        comb_span = delta_frequency_to_delta_wavelength(
            high_freq_modulation*(number_of_teeth-1), c / laser_wavelength * 1e9)

        if comb_span > wl_max - wl_min or comb_span < min_comb_span:
            continue

        # Print the parameters
        print(
            f'Number of teeth: {number_of_teeth}, High frequency modulation: {high_freq_modulation/1e9:.2f} GHz')

        config_results = []

        # Create a folder in `figures` with name like: 8 x 1.0 GHz
        folder_name = f'{number_of_teeth} x {high_freq_modulation / 1e9:.1f} GHz'
        folder_path = os.path.join(get_figures_path(), folder_name)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        for i in range(10):
            # Simulate the transmission spectrum.
            x_meas, y_meas, f = fit_simulated_measurement_concentration(
                molecule=molecule, wl_min=wl_min, wl_max=wl_max, vmr=vmr, pressure=pressure, temperature=temperature,
                length=length, laser_wavelength=laser_wavelength, high_freq_modulation=high_freq_modulation,
                number_of_teeth=number_of_teeth, database=database, std_dev=std_dev, x_shift_std_dev=x_shift_std_dev,
                scaling_std_dev=scaling_std_dev, laser_wavelength_std_dev=laser_wavelength_std_dev, normalize=normalize,
                initial_guess=initial_guess)

            meas = f.measured_transmission
            sim = f.simulated_transmission

            x_sim, y_sim, x_meas_fitted, y_meas_fitted = sim.x_nm, sim.y_nm, meas.x_nm, meas.y_nm

            config_results.append(f.concentration)

            # Plot the simulated and measured transmission spectra.

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

            print(f'The fitted concentration of {molecule} is {f.concentration:.6f} VMR.')

        config_results = array(config_results)
        result = (number_of_teeth, high_freq_modulation,
                  config_results.mean(), config_results.std())

        results.append(result)

        # Print the result
        print(f'Number of teeth: {number_of_teeth}, High frequency modulation: {high_freq_modulation/1e9:.2f} GHz, ' +
              f'Mean concentration: {result[2]:.6f} VMR, Standard deviation: {result[3]:.6f} VMR')
        print('--------------------')


timestr = time.strftime("%Y%m%d-%H%M%S")
with open(f'{get_reports_path()}report-{timestr}.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(('Number of teeth', 'High frequency modulation (Hz)',
                    'Mean concentration (VMR)', 'Standard deviation (VMR)'))
    writer.writerows(results)
