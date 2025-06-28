from ase.io import read
from ase.geometry import find_mic
import numpy as np
import argparse

def find_atoms_within_radius(filename, R0, Rcut):
    atoms = read(filename, format='espresso-in')

    positions = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.get_pbc()

    # Compute distance with minimum image convention (ASE function)
    deltas, dists = find_mic(positions - R0, cell=cell, pbc=pbc)

    indices_within_cutoff = np.where(dists <= Rcut)[0]

    return atoms, indices_within_cutoff, positions[indices_within_cutoff], dists[indices_within_cutoff]

def main():
    parser = argparse.ArgumentParser(description="Find atom indices within radius R from point R0 using ASE and PBC.")
    parser.add_argument("--input", required=True, help="QE input file")
    parser.add_argument("--R0", nargs=3, type=float, required=True, help="Reference point (in Å)")
    parser.add_argument("--R", type=float, required=True, help="Cutoff radius (in Å)")
    args = parser.parse_args()

    atoms, indices, pos_list, dists = find_atoms_within_radius(args.input, np.array(args.R0), args.R)

    print(f"\nAtoms within {args.R:.2f} Å of point {args.R0}:")
    
    indexes_list = []
    
    for i, pos, dist in zip(indices, pos_list, dists):
        symbol = atoms[i].symbol
        print(f"Index {i+1:3d}: {symbol:>2} at {pos} (distance = {dist:.3f} angstroms)")
        indexes_list.append(i+1)
        
    print('Total of atoms in this file: ', len(atoms))
    print('Total of atoms within radius: ', len(indexes_list))
        
    print(f"\nIndexes: {indexes_list}")

if __name__ == "__main__":
    main()
