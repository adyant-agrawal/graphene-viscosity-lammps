# Viscosity of graphene suspensions in shear flow

These scripts compute the viscosity of graphene nanoparticle suspensions in water with molecular dynamics, following A. Agrawal, C. Kamal, S. Gravelle and L. Botto, [Negative Intrinsic Viscosity in Graphene Nanoparticle Suspensions Induced by Hydrodynamic Slip](https://doi.org/10.1021/acsnano.5c10415), ACS Nano 19, 38350 (2025). The particles are either 3D discs, free to rotate, or quasi-2D plates that span the periodic y direction and rotate only about y. They are placed in TIP4P/2005 water between two rock-salt FeO walls. After equilibration at $10^5\,\text{Pa}$, the walls move at opposite velocities along x to impose a Couette flow, and the viscosity follows from the shear force on the walls.

Slip at the graphene surface is set by the LJ energy between carbon and the water oxygen, `eps_CO`: the default 0 keeps the mixing rule (pristine graphene, large slip), and 1.4 kcal/mol gives no slip.

## Layout

```
3D-discs/build.py     builds discs in water between the walls
2D-plates/build.py    builds plates in water between the walls
shared/system.py      molecule and wall builders used by both build.py
shared/PARM.lammps    force field
shared/input.equilibration.lammps
shared/input.shear.lammps
shared/analysis.py    viscosity from the wall force
```

## Requirements

- Python with numpy, ase, MDAnalysis and solvate: `pip install numpy ase MDAnalysis solvate` (tested with numpy 2.5.3, ase 3.29.0, MDAnalysis 2.10.0 and solvate 0.1.0)
- LAMMPS with the MOLECULE, KSPACE, RIGID and EXTRA-DUMP packages (tested with the 28 Mar 2023 version)

## Usage

A run happens inside `3D-discs/` or `2D-plates/`, or inside a copy of one of them:

```
python build.py -n 3 -d 50
lmp -in ../shared/input.equilibration.lammps
lmp -in ../shared/input.shear.lammps -var gamma 5e10 -var eps_CO 1.4
python ../shared/analysis.py --water ../water
```

The last command needs a pure-water reference run with the same box. Since `-d` also sets the default `-H` and `-L`, build it with the same `-d` (or the same `-H` and `-L`) as the particle run:

```
cp -r 3D-discs water && cd water && python build.py -n 0 -d 50
```

Use `2D-plates` instead of `3D-discs` for a plate run, then run the two `lmp` commands in `water/` with the same `-var gamma` as the particle run. From the particle folder, `analysis.py --water ../water` prints the suspension viscosity eta, the water viscosity eta0, eta/eta0, and (eta/eta0 - 1)/c, with c the particle concentration. Without `--water`, it prints eta only.

## Options

`build.py` (both folders):

| option | default | meaning |
|--------|---------|---------|
| `-n` | 1 | number of particles; 0 gives pure water |
| `-d` | 50 | disc diameter or plate length 2a, in Å |
| `-H` | 3*d | channel height, wall to wall, in Å |
| `-L` | 2*d | box size in Å: Lx = Ly = L for discs; Lx = L for plates, whose Ly is fixed at 21.08 Å |
| `--seed` | 1 | random seed for the particle placement |

Particles are kept at least 12 Å from the walls and from each other.

`input.shear.lammps`, set with `-var name value`: `gamma` is the shear rate in 1/s (default 5e10), `eps_CO` the C-O(water) LJ energy in kcal/mol (default 0), and `nsteps` the number of production steps (default 1000000, with a 1 fs timestep).

## Output files

`build.py` writes `data.lammps` and `system.json`. The equilibration writes `data.equilibrium`, which the shear run reads, and the diagnostics `box.equilibration.dat`, `temperature.equilibration.dat` and `dump.equilibration.xtc`. The shear run writes `force.dat`, `box.dat`, `temperature.dat`, `gamma.dat`, the trajectory `dump.shear.xtc` and the final state `data.shear`; `analysis.py` reads `force.dat`, `box.dat`, `gamma.dat` and `system.json`.

## Citation

If you use these scripts, please cite:

A. Agrawal, C. Kamal, S. Gravelle, L. Botto, "Negative Intrinsic Viscosity in Graphene Nanoparticle Suspensions Induced by Hydrodynamic Slip", ACS Nano 19, 38350 (2025), https://doi.org/10.1021/acsnano.5c10415

## Licence

GPL-3.0, see `LICENSE`.
