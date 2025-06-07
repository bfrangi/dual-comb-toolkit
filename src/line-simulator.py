from lib.simulations import Simulator

# Define the molecule, VMR, pressure, temperature, and length.

molecule = "CH4"
vmr = 0.01  # volume mixing ratio
pressure = 53328.94736842  # Pa
temperature = 298  # K
length = 0.055  # m

# Simulation configuration.

wl_min = 3427.1  # nm
wl_max = 3427.8  # nm
use_gpu = False

# Simulate the transmission spectrum.

s = Simulator(
    molecule=molecule,
    vmr=vmr,
    pressure=pressure,
    temperature=temperature,
    length=length,
    use_gpu=use_gpu,
)
s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max)
s.show_transmission_spectrum()
