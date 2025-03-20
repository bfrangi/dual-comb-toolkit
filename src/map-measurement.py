from lib.files import get_measurement_names
from lib.mapping import Mapper

# Define the molecule, VMR, pressure, temperature, and length.
molecule = 'CH4'
pressure = 53328.94736842  # Pa
temperature = 298  # K
length = 0.055  # m
wl_min = 3427.1  # nm
wl_max = 3427.8  # nm

# Specify the name and specifications of the measurement.
sweep_name = 'cell-sweep-10-34-17-03-2025'
meas_names = get_measurement_names(sweep_name)
center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
number_of_teeth = 38
laser_wavelength = 3427.437e-9  # m
high_freq_modulation = 500e6  # Hz
acq_freq = 400000.0  # Hz

# Map the measurements.

mapper = Mapper(
    meas_names,
    center_freq=center_freq,
    freq_spacing=freq_spacing,
    number_of_teeth=number_of_teeth,
    laser_wavelength=laser_wavelength,
    high_freq_modulation=high_freq_modulation,
    acq_freq=acq_freq,
    molecule=molecule,
    pressure=pressure,
    temperature=temperature,
    length=length,
    wl_min=wl_min,
    wl_max=wl_max,
    initial_guess={'concentration': 0.5}
)

mapper.show_concentration_heatmap()

mapper.show_concentration_plot()
