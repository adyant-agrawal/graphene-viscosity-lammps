"""Viscosity of the suspension from the shear force on the walls.

Run in 3D-discs/ or 2D-plates/ after input.shear.lammps:
    python ../shared/analysis.py                    viscosity eta
    python ../shared/analysis.py --water ../water   also eta/eta0 and (eta/eta0 - 1)/c
"""
import argparse
import json
import numpy as np

KCAL_MOL_A_TO_N = 4184 / 6.02214076e23 / 1e-10  # LAMMPS force unit (kcal/mol/A) in N
NBLOCKS = 5  # blocks for the standard error


def viscosity(folder):
    """Viscosity eta (mPa s) and its standard error, from the output of input.shear.lammps in folder."""
    # x force on the top and bottom walls (kcal/mol/A), box size (A), shear rate (1/s)
    _, f_top, f_bot = np.loadtxt(f"{folder}/force.dat", unpack=True)
    _, _, _, lx, ly = np.loadtxt(f"{folder}/box.dat", unpack=True)
    gamma = float(open(f"{folder}/gamma.dat").read())

    # shear stress (Pa): mean x force of the two walls per wall area
    stress = (f_bot - f_top) / 2 * KCAL_MOL_A_TO_N / (lx * ly * 1e-20)

    # skip the first 20 % (start-up of the flow)
    stress = stress[len(stress) // 5:]
    if len(stress) < NBLOCKS:
        raise SystemExit(f"{folder}: only {len(stress)} force samples, run longer")

    # eta = stress / shear rate, in mPa s, averaged over blocks
    eta = np.array([block.mean() for block in np.array_split(stress, NBLOCKS)]) / gamma * 1e3
    return eta.mean(), eta.std(ddof=1) / np.sqrt(NBLOCKS)


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("folder", nargs="?", default=".", help="run folder (default: here)")
parser.add_argument("--water", help="folder of the pure-water run (built with -n 0)")
args = parser.parse_args()

eta, err = viscosity(args.folder)
print(f"eta = {eta:.4f} +- {err:.4f} mPa s")

# relative and reduced viscosity, using the pure-water run and the volume fraction c from system.json
if args.water:
    eta0, err0 = viscosity(args.water)
    with open(f"{args.folder}/system.json") as f:
        c = json.load(f)["c"]
    print(f"eta0 = {eta0:.4f} +- {err0:.4f} mPa s")
    print(f"eta/eta0 = {eta / eta0:.4f}")
    print(f"(eta/eta0 - 1)/c = {(eta / eta0 - 1) / c:.2f}")
