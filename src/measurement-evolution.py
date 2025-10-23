import os

from matplotlib import animation
from matplotlib import pyplot as plt

from lib.entities import MeasuredSpectrum
from lib.files import get_animations_path
from lib.fitting import ConcentrationFitter
from lib.measurements import Measurement
from lib.plots import use_latex

####################################################################################################
# Measurement parameters                                                                           #
####################################################################################################

# Molecule and physical conditions.

molecule = "CH4"
database = "hitemp"
pressure = 101325  # Pa
temperature = 1180.235  # K
length = 0.07  # m

# Simulation range.

wl_min = 3426.8  # nm
wl_max = 3428.1  # nm

# Measurement name.

measurement_name = "8a/Position-X1-Y10"
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

# Fitter, initial guess and allowed concentration bounds.

fitter = "normal_gpu"
initial_guess = 0.0001  # VMR
lower_bound = 0.0  # VMR
upper_bound = 0.15  # VMR

# Noise filtering.

tooth_std_threshold = 10  # Teeth with a standard deviation above `tooth_std_threshold * mean_std` will be discarded.
sub_measurements = (
    0  # Number of sub-measurements used to obtain the standard deviation of the teeth.
)

# Removing noisy teeth

remove_teeth_indices = [1, 12]  # List of tooth indices to be removed from the fitting.

# Output and plotting parameters.

verbose = True
spectrum_plot_folder = "fit-measurement-output"

# Use LaTeX for plotting.

# use_latex()

####################################################################################################
# Fitting                                                                                          #
####################################################################################################

integration_time = 0.1  # s
time_between_measurements = 0.05  # s


m = Measurement(
    measurement_name=measurement_name,
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
)

total_time = m.measurement_time

intervals = [
    (i * time_between_measurements, i * time_between_measurements + integration_time)
    for i in range(int(total_time / time_between_measurements))
]

print(f"Total measurement time: {total_time:.2f} s")

measured_spectra: list[MeasuredSpectrum] = []
concentrations = []
fitted_spectra = []
min_abs = 10
for start, end in intervals:
    m.compute_transmission(start=start, end=end, save=True)
    m.remove_teeth(remove_teeth_indices)
    measured_spectrum = m.transmission_spectrum
    f = ConcentrationFitter(
        meas_transmission=measured_spectrum,
        wl_min=wl_min,
        wl_max=wl_max,
        initial_guess=initial_guess,
        fitter=fitter,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        verbose=verbose,
        database=database,
    )
    concentrations.append(f.concentration)
    fitted_spectra.append(f.simulated_transmission)
    measured_spectra.append(f.measured_transmission)
    min_abs = min(min_abs, measured_spectra[-1].y_nm.min())


fig = plt.figure(figsize=(6, 4.5), dpi=(1920 / 6))
ax = plt.gca()  # Get current axes
sim = ax.plot(
    fitted_spectra[0].x_nm,
    fitted_spectra[0].y_nm,
    color="blue",
    label="Simulated Transmission",
)[0]
meas = ax.plot(
    measured_spectra[0].x_nm,
    measured_spectra[0].y_nm,
    "o-",
    color="red",
    label="Measured Transmission",
)[0]
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Transmission [-]")
ax.set_title(
    f"{measurement_name}\n$t = 0$ s\t Concentration = {concentrations[0]:.6f} VMR"
)
ax.set_xlim(wl_min, wl_max)
ax.set_ylim(min_abs * 0.96, 1.03)
ax.legend(loc="upper right")


def update(frame):
    # for each frame, update the data stored on each artist.
    measured_spectrum = measured_spectra[frame]
    x_meas = measured_spectrum.x_nm
    y_meas = measured_spectrum.y_nm
    x_sim = fitted_spectra[frame].x_nm
    y_sim = fitted_spectra[frame].y_nm
    meas.set_xdata(x_meas)
    meas.set_ydata(y_meas)
    sim.set_xdata(x_sim)
    sim.set_ydata(y_sim)
    ax.set_title(
        f"{measurement_name}\n$t = {frame * time_between_measurements:.2f}$ s\t Concentration = {concentrations[frame]:.6f} VMR"
    )
    return meas


animation_name = f"{measurement_name.replace('/', '-')} - it{integration_time}s tbm{time_between_measurements}s.gif"
animation_path = os.path.join(get_animations_path(), animation_name)

ani = animation.FuncAnimation(
    fig=fig,
    func=update,
    frames=len(measured_spectra),
    interval=time_between_measurements * 1000,
)
ani.save(filename=animation_path, writer="pillow")
plt.show()
