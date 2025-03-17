from lib.simulations import Simulator

# Define the molecule, VMR, pressure, temperature, and length.
molecule = 'CH4'
vmr = 0.01 # volume mixing ratio
pressure = 53328.94736842 # Pa
temperature = 298 # K
length = 0.055 # m
wl_min = 3426.0 # nm
wl_max = 3428.8 # nm


s = Simulator(
    molecule=molecule,
    vmr=vmr,
    pressure=pressure,
    temperature=temperature,
    length=length,
)
s.compute_transmission_spectrum(wl_min=wl_min, wl_max=wl_max)
s.show_transmission_spectrum()