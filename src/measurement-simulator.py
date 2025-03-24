from matplotlib import pyplot as plt

from lib.plots import tight
from lib.simulations import simulate_measurement

# Define the molecule, VMR, pressure, temperature, and length.
molecule = 'CH4'
vmr = 0.01  # volume mixing ratio
pressure = 53328.94736842  # Pa
temperature = 298  # K
length = 0.055  # m
wl_min = 3427.1  # nm
wl_max = 3427.8  # nm
database = 'hitran'

# Specify the name and specifications of the measurement.
number_of_teeth = 38
laser_wavelength = 3427.45  # nm
high_freq_modulation = 500e6  # Hz
std_dev = 0.001

# Simulate the transmission spectrum.
x_sim, y_sim = simulate_measurement(molecule=molecule, wl_min=wl_min,
                                    wl_max=wl_max, laser_wavelength=laser_wavelength,
                                    high_freq_modulation=high_freq_modulation,
                                    number_of_teeth=number_of_teeth,
                                    vmr=vmr,
                                    pressure=pressure,
                                    temperature=temperature,
                                    length=length,
                                    database=database,
                                    std_dev=std_dev)


# Plot the simulated transmission spectrum.
plt.scatter(x_sim, y_sim, label='Simulated Measurement', color='blue', zorder=0)
plt.legend()
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission')
plt.title(f'Transmission spectrum of {molecule} at {pressure:.2f} Pa ' +
          f'and {temperature:.2f} K.\n{length:.3f} m path length, {vmr:.3f} VMR.')
plt.tight_layout(**tight)
plt.show()
