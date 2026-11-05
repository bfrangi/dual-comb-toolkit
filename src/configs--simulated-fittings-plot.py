"""Paper-style configuration plots for a simulated-fittings run, with selected configs highlighted.

Reuses the per-condition summary plotter (``lib.simulated_fittings.summary_plots``) and just adds
markers on the chosen ``(number_of_teeth, comb_spacing_GHz)`` configurations. For every condition
in the run it writes ``paper-conc.{svg,pdf,png}`` and ``paper-sdv.{svg,pdf,png}`` next to that
condition's ``results.csv`` (alongside, not overwriting, the auto-generated ``summary-*`` figures).

Edit ``RUN_ID`` / ``SELECTED_CONFIGS`` below and run from the ``src`` directory:

    python configs--simulated-fittings-plot.py
"""

from pathlib import Path

from lib.files import get_simulations_path
from lib.simulated_fittings.config import load_config
from lib.simulated_fittings.summary_plots import make_summary_plots

####################################################################################################
# What to plot                                                                                     #
####################################################################################################

RUN_ID = "run-20260616-151858"
SELECTED_CONFIGS = [
    (19, 0.8),
    # (23, 0.8),
]  # (number of teeth, comb spacing in GHz) to highlight


def main() -> None:
    run_dir = Path(get_simulations_path()) / RUN_ID
    if not run_dir.is_dir():
        raise SystemExit(f"Run not found: {run_dir}")

    # Molecule / true VMR per condition come from the run's config snapshot.
    cfg = load_config(run_dir / "config.snapshot.py")
    by_id = {c["id"]: c for c in cfg.conditions}

    for results in sorted(run_dir.glob("*/results.csv")):
        cond_dir = results.parent
        cond = by_id.get(cond_dir.name)
        if cond is None:
            print(f"  ! {cond_dir.name}: not in snapshot, skipping")
            continue

        stems = make_summary_plots(
            str(results),
            str(cond_dir),
            molecule=cond["molecule"],
            vmr_true=cond["vmr"],
            highlight=SELECTED_CONFIGS,
            stem_prefix="paper",
        )
        names = ", ".join(Path(s).name for s in stems) or "(no data)"
        print(f"  {cond_dir.name}: {names}")


if __name__ == "__main__":
    main()
