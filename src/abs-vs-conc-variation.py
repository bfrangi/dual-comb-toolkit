from labellines import labelLines
from matplotlib import pyplot as plt

from lib.files import get_figures_path
from lib.math import get_x_of_min_y
from lib.plots import article_figsize, article_tight, get_cmap, use_latex
from lib.simulations import Simulator

####################################################################################################
# Simulation parameters                                                                            #
####################################################################################################

# Molecule and database.

molecule = "CH4"
"""Molecule to simulate."""
database = "hitemp"
"""Database to use for the simulation. Can be 'hitran', 'hitemp', 'exomol' or 'geisa'. Some may not
be available for all molecules."""

# Physical conditions.

pressure = 101325  # Pa
"""Pressure of the gas mixture."""
temperature = 1200  # K
"""Temperature of the gas mixture."""
length = 0.07  # m
"""Path length of the gas mixture."""
vmrs = [0.005*(i+1) for i in range(10)] # VMR
"""List of volume mixing ratios to simulate."""

# Absorption feature ranges.
ranges = [
    # (3427.0, 3427.9),  # nm
    (3391, 3391.83),
    (3391.83, 3392.14),
    (3392.14, 3392.44),
    (3392.44, 3393),
    (3402, 3403.75),
    (3403.75, 3404.1),
    (3404.1, 3404.38),
    (3404.38, 3405.3),
    (3414.5, 3415.68),
    (3415.68, 3416.2),
    (3416.2, 3416.5),
    (3416.5, 3417.6),
    (3426.8, 3427.82),
    (3427.82, 3428.35),
    (3428.35, 3428.68),
    (3428.68, 3429.8),
    (3439.2, 3440),
]
"""List of (minimum wavelength, maximum wavelength) tuples for the absorption feature ranges."""

use_gpu = True
"""Whether to use the GPU for the simulation."""

####################################################################################################
# Output parameters                                                                                #
####################################################################################################

# Output folder for plots.
output_folder = "abs-vs-conc-variation-output"
"""Folder to save the output plots."""

# Latex for plotting.

latex = True
if latex:
    use_latex()

####################################################################################################
# Simulate                                                                                         #
####################################################################################################

results: dict[tuple[float, float], dict[str, list[float] | float]] = {}
sims: list[Simulator] = []

figures_path = f"{get_figures_path()}{output_folder}"

wl_min = min(r[0] for r in ranges)
wl_max = max(r[1] for r in ranges)

for i, vmr in enumerate(vmrs):
    print(
        f"Preparing simulator for {vmr} VMR ({i + 1}/{len(vmrs)})..."
    )

    sim = Simulator(
        molecule=molecule,
        database=database,
        vmr=vmr,
        pressure=pressure,
        temperature=temperature,
        length=length,
        use_gpu=use_gpu,
    )
    sim.compute_transmission_spectrum(wl_min, wl_max)
    sims.append(sim)

for i, r in enumerate(ranges):
    print(f"Processing range {r[0]}-{r[1]} nm ({i + 1}/{len(ranges)})...")

    wl_min, wl_max = r
    results[r] = {
        "peak_abs": [],
        "central_wl": None,
    }

    # Iterate over each concentration and temperature
    for j, sim in enumerate(sims):
        wl, tr = sim.get_transmission_spectrum(wl_min, wl_max)
        peak_absorption = 1 - tr.min()
        results[r]["peak_abs"].append(peak_absorption)

        if results[r]["central_wl"] is None:
            results[r]["central_wl"] = get_x_of_min_y(wl, tr)

ranges.sort(key=lambda ran: results[ran]["central_wl"])
plt.figure(figsize=article_figsize)

cmap = get_cmap("CMRmap", len(ranges))

for i, result in enumerate(results.values()):
    plt.plot(
        vmrs,
        result["peak_abs"],
        "o-",
        label=f"{result['central_wl']:.2f} nm",
        color=cmap(i),
    )

labelLines(plt.gca().get_lines(), zorder=2.5)
plt.tight_layout(**article_tight)
plt.xlabel("Concentration [VMR]")
plt.ylabel("Peak Absorption [-]")
# plt.title(
#     f"Peak Absorption vs Concentration for {molecule}\n{temperature} K, {pressure} Pa, {length} m path length"
# )
plt.savefig(f"{figures_path}/abs_vs_conc.pdf", dpi=300)
plt.savefig(f"{figures_path}/abs_vs_conc.svg", dpi=300)
plt.clf()
