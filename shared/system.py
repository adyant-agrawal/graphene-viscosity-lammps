"""Building blocks for graphene particles in water between two FeO walls.

Types, numbered as in PARM.lammps:

    atom types   1 water O, 2 water H, 3 wall Fe, 4 wall O, 5 graphene C, 6 graphene H
    bond types   1 water O-H, 2 graphene C-C, 3 graphene C-H
    angle types  1 water H-O-H, 2 graphene C-C-C, 3 graphene C-C-H
    dihedral 1 and improper 1: graphene
"""
import re
from itertools import combinations
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, minimize_vectors
from ase.build import graphene

ACC = 1.40      # C-C distance (A), equal to the OPLS bond length
ACH = 1.08      # C-H distance (A)
DFEO = 2.085    # Fe-O distance in rock-salt FeO (A)
NLAYER = 6      # atomic layers per wall
VACUUM = 30.0   # empty gap between the two walls across the periodic z boundary (A)
LY = 21.08      # box length along y for plates (A): 5 armchair periods of graphene, 5 FeO cells
NTYPES = {"atom": 6, "bond": 3, "angle": 3, "dihedral": 1, "improper": 1}


def universe(types, charges, positions, resindex, resname, box):
    """MDAnalysis Universe of atoms with types, charges and positions; resindex numbers the molecules from 0."""
    nres = max(resindex) + 1
    u = mda.Universe.empty(len(types), nres, atom_resindex=resindex, trajectory=True)
    u.add_TopologyAttr("types", types)
    u.add_TopologyAttr("names", types)
    u.add_TopologyAttr("charges", charges)
    u.add_TopologyAttr("resids", np.arange(1, nres + 1))
    u.add_TopologyAttr("resnames", [resname] * nres)
    for attr in ("bonds", "angles", "dihedrals", "impropers"):  # empty, so that mda.Merge keeps them
        u.add_TopologyAttr(attr, [])
    u.atoms.positions = positions
    u.dimensions = box
    return u


def graphene_sheet(size):
    """Carbon positions of a graphene sheet in the x-y plane, larger than size (A), centred at the origin."""
    n = int(size / ACC)
    sheet = graphene(a=ACC * 3**0.5, size=(2 * n, 2 * n, 1)).positions
    return sheet - sheet.mean(axis=0)


def make_disc(d):
    """Graphene disc of diameter d (A), including 2.5 A for the edge atoms."""
    sheet = graphene_sheet(d)
    distance_from_centre = np.linalg.norm(sheet[:, :2], axis=1)
    carbons = sheet[distance_from_centre <= (d - 5) / 2]
    no_periodicity = np.array([1000, 1000, 1000, 90, 90, 90])  # box much larger than the disc
    return add_topology(add_hydrogens(carbon_flake(carbons, no_periodicity), no_periodicity))


def make_plate(d):
    """Graphene plate of length d (A) along x, including 2.5 A for the edge atoms, periodic along y."""
    sheet = graphene_sheet(max(d, LY))
    # place an atom at y = 0.2 A, so that no atom lies on the cut at y = 0
    closest_to_centre = sheet[np.argmin(np.linalg.norm(sheet, axis=1))]
    sheet = sheet - closest_to_centre + [0, 0.2, 0]
    period = 15 * ACC  # 5 periods of graphene along y, 21.0 A, stretched below to LY
    keep = (np.abs(sheet[:, 0]) <= (d - 5) / 2) & (sheet[:, 1] >= 0) & (sheet[:, 1] < period)
    carbons = sheet[keep]
    carbons[:, 1] *= LY / period
    periodic_in_y = np.array([1000, LY, 1000, 90, 90, 90])
    return add_topology(add_hydrogens(carbon_flake(carbons, periodic_in_y), periodic_in_y))


def carbon_flake(carbons, box):
    """Carbon atoms (type 5) with C-C bonds (type 2) when closer than 1.1 ACC; dangling carbons are removed."""
    n = len(carbons)
    u = universe(["5"] * n, np.zeros(n), carbons, [0] * n, "GRA", box)
    pairs = capped_distance(carbons, carbons, 1.1 * ACC, box=box, return_distances=False)
    bonds = pairs[pairs[:, 0] < pairs[:, 1]]  # each pair once
    u.add_bonds(bonds, types=["2"] * len(bonds))

    kept = u.atoms[[len(atom.bonds) >= 2 for atom in u.atoms]]  # dangling carbons have one neighbour
    flake = mda.Merge(kept)  # keeps the bonds between the kept atoms
    flake.dimensions = box
    return flake


def add_hydrogens(flake, box):
    """Add an H (type 6) to every edge carbon, i.e. every carbon with two neighbours."""
    edge = [atom for atom in flake.atoms if len(atom.bonds) == 2]
    hydrogens = []
    for atom in edge:  # the H points away from the two neighbours
        to_neighbours = minimize_vectors(atom.bonded_atoms.positions - atom.position, box)
        outward = -to_neighbours.sum(axis=0)
        hydrogens.append(atom.position + ACH * outward / np.linalg.norm(outward))

    nc, nh = flake.atoms.n_atoms, len(hydrogens)  # carbons first, then hydrogens
    charges = np.r_[np.zeros(nc), np.full(nh, 0.115)]
    charges[[atom.index for atom in edge]] = -0.115
    u = universe(["5"] * nc + ["6"] * nh, charges, np.vstack([flake.atoms.positions] + hydrogens), [0] * (nc + nh), "GRA", box)
    u.add_bonds(flake.bonds.indices, types=["2"] * len(flake.bonds))
    u.add_bonds([(atom.index, nc + n) for n, atom in enumerate(edge)], types=["3"] * nh)
    return u


def add_topology(u):
    """Add the angles, dihedrals and impropers of a graphene flake from its bonds."""
    angles, angle_types, dihedrals, impropers = [], [], [], []

    for atom in u.atoms:
        neighbours = atom.bonded_atoms
        for a, b in combinations(neighbours, 2):
            angles.append((a.index, atom.index, b.index))
            if a.type == "6" or b.type == "6":
                angle_types.append("3")  # C-C-H
            else:
                angle_types.append("2")  # C-C-C
        if len(neighbours) == 3:
            impropers.append((atom.index, *neighbours.indices))  # central atom first

    for bond in u.bonds:
        b, c = bond.atoms
        for a in b.bonded_atoms:
            for d in c.bonded_atoms:
                if a != c and d != b:
                    dihedrals.append((a.index, b.index, c.index, d.index))

    u.add_angles(angles, types=angle_types)
    u.add_dihedrals(dihedrals, types=["1"] * len(dihedrals))
    u.add_impropers(impropers, types=["1"] * len(impropers))
    return u


def make_walls(lx, ly, H):
    """Two FeO walls with inner planes H apart; returns the Universe and the z of both inner planes."""
    nx, ny = 2 * round(lx / DFEO / 2), 2 * round(ly / DFEO / 2)
    dx, dy = lx / nx, ly / ny
    zbot = VACUUM / 2 + (NLAYER - 1) * DFEO
    ztop = zbot + H

    # one atom per lattice site (i, j, k); layers k = 0..5 form the bottom wall, k = 6..11 the top wall
    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(2 * NLAYER), indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()
    top = k >= NLAYER
    z = VACUUM / 2 + k * DFEO + top * (H - DFEO)
    iron = (i + j + k) % 2 == 0  # Fe and O alternate on the cubic lattice

    box = [nx * dx, ny * dy, ztop + (NLAYER - 1) * DFEO + VACUUM / 2, 90, 90, 90]
    u = universe(np.where(iron, "3", "4"), np.where(iron, 2.0, -2.0), np.c_[i * dx, j * dy, z],
                 top.astype(int), "WAL", box)
    return u, zbot, ztop


def tip4p2005():
    """One TIP4P/2005 water molecule (O, H, H; LAMMPS adds the M site)."""
    theta = np.deg2rad(104.52)
    positions = [[0, 0, 0], [0.9572, 0, 0], [0.9572 * np.cos(theta), 0.9572 * np.sin(theta), 0]]
    u = universe(["1", "2", "2"], [-1.1128, 0.5564, 0.5564], positions, [0, 0, 0], "SOL", [10, 10, 10, 90, 90, 90])
    u.add_bonds([(0, 1), (0, 2)], types=["1", "1"])
    u.add_angles([(1, 0, 2)], types=["1"])
    return u


def insert_mol(u, mol, n, axis=None, zmin=0, zmax=None, distance=12, tries=1000):
    """Insert n copies of mol at random positions and orientations, at least distance (A) from other atoms.

    axis fixes the rotation axis (None: random axis).
    """
    box = u.dimensions
    for k in range(n):
        mol.residues.resids = u.residues.resids.max() + 1  # each copy is a new molecule
        new = mda.Merge(u.atoms, mol.atoms)
        new.dimensions = box
        copy = new.atoms[u.atoms.n_atoms:]
        existing = mda.lib.NeighborSearch.AtomNeighborSearch(new.atoms[:u.atoms.n_atoms], box=box)
        for _ in range(tries):
            centre = np.random.uniform([0, 0, zmin], [box[0], box[1], zmax or box[2]])
            copy.translate(centre - copy.center_of_geometry())
            rotation_axis = np.random.randn(3) if axis is None else axis
            copy.rotateby(np.random.uniform(0, 360), rotation_axis, point=copy.center_of_geometry())
            if not existing.search(copy, distance):  # no atom of u within distance
                break
        else:
            raise RuntimeError(f"no room for particle {k + 1}: increase -H or -L, or decrease -n")
        u = new
    return u


def write_data(u, filename):
    """Write u as a LAMMPS data file for PARM.lammps."""
    u.trajectory.ts.data["molecule_tag"] = u.atoms.resids  # molecule IDs
    u.atoms.write(filename, file_format="DATA")
    with open(filename) as f:
        text = f.read()
    for key, n in NTYPES.items():  # all types of PARM.lammps, also when absent (e.g. pure water)
        text = re.sub(rf"\d+\s+{key} types", f"{n} {key} types", text)
    text = re.sub(r"\nMasses\n\n.*?\n\n", "\n", text, flags=re.S)  # masses are in PARM.lammps
    with open(filename, "w") as f:
        f.write(text)
