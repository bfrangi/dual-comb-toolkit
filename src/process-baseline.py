from lib.entities import Baseline

# Define the molecule, VMR, pressure, temperature, and length.
molecule = 'CH4'
vmr = 0.01  # volume mixing ratio
pressure = 53328.94736842  # Pa
temperature = 298  # K
length = 0.055  # m
wl_min = 3427.1  # nm
wl_max = 3427.8  # nm

# Specify the name and specifications of the measurement.
meas_baseline = [
    'cell-sweep-10-36-17-03-2025/Position-X12-Y1',
    'cell-sweep-10-36-17-03-2025/Position-X13-Y1',
    'cell-sweep-10-36-17-03-2025/Position-X14-Y1',
    'cell-sweep-10-36-17-03-2025/Position-X15-Y1',
    'cell-sweep-10-36-17-03-2025/Position-X16-Y1',
]
center_freq = 40000.0  # Hz
freq_spacing = 200.0  # Hz
number_of_teeth = 38
laser_wavelength = 3427.45e-9  # m
high_freq_modulation = 500e6  # Hz
acq_freq = 400000.0  # Hz

# Compute the baseline.
bl = Baseline(measurement_names=meas_baseline, center_freq=center_freq,
              freq_spacing=freq_spacing, number_of_teeth=number_of_teeth,
              laser_wavelength=laser_wavelength, high_freq_modulation=high_freq_modulation,
              acq_freq=acq_freq)

bl.show_baseline_plot()

bl.show_baselines_plot()

bl.show_standard_deviation_plot()
