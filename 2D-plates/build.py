"""Build graphene plates (periodic along y) in water between two FeO walls, as a LAMMPS data file.

Example: python build.py -n 3 -d 50
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import solvate

# building blocks in ../shared/system.py
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "shared"))
from system import LY, make_plate, make_walls, tip4p2005, insert_mol, write_data

# options; H and L default to multiples of d
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-n", type=int, default=1, help="number of plates, 0 for pure water (default 1)")
parser.add_argument("-d", type=float, default=50, help="plate length 2a in A (default 50)")
parser.add_argument("-H", type=float, help="channel height, wall to wall, in A (default 3*d)")
parser.add_argument("-L", type=float, help="box size along x in A (default 2*d); along y it is fixed")
parser.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
args = parser.parse_args()
H = args.H or 3 * args.d
L = args.L or 2 * args.d
np.random.seed(args.seed)

# 1. two FeO walls at zbot and ztop, H apart; box length LY along y
u, zbot, ztop = make_walls(L, LY, H)

# 2. n plates, randomly placed and rotated about y, with 12 A clearance from walls and other plates
u = insert_mol(u, make_plate(args.d), args.n, axis=[0, 1, 0], zmin=zbot + 12, zmax=ztop - 12)

# 3. water at 1 g/cm3 (0.0334 per A^3)
u = solvate.SolvatePlanar(u, tip4p2005(), density=0.0334, zmin=zbot, zmax=ztop, distance=2.0)

# 4. wrap into the box and write the LAMMPS data file
u.atoms.wrap()
write_data(u, "data.lammps")

# 5. box size and plate volume fraction c (thickness 2b = 5 A), saved to system.json
lx, ly = u.dimensions[:2]
c = float(args.n * args.d * 5 * ly / (lx * ly * H))
info = {"n": args.n, "d": args.d, "H": H, "lx": float(lx), "ly": float(ly), "c": c}
with open("system.json", "w") as f:
    json.dump(info, f, indent=1)

nwater = u.select_atoms("type 1").n_atoms
print(f"data.lammps: {args.n} plates of d = {args.d} A, box {lx:.2f} x {ly:.2f} A, H = {H} A, {nwater} waters, c = {c:.4f}")
