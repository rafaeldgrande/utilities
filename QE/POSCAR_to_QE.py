
from ase.io import read, write
import argparse

'''Usage:
python POSCAR_to_QE.py --vasp_file POSCAR --qe_input qe_input.in --E_cutoff 40.0 --conv_thr 1e-10 --nbnds 1
'''

parser = argparse.ArgumentParser()
parser.add_argument('--vasp_file', type=str, default='POSCAR', help='POSCAR file from VASP')
parser.add_argument('--qe_input', type=str, default='qe_input.in', help='Quantum ESPRESSO input file to be generated using data from POSCAR file')
parser.add_argument('--E_cutoff', type=float, default=40.0, help='Cutoff energy in Ry')
parser.add_argument('--conv_thr', type=float, default=1e-10, help='Convergence threshold in Ry')
parser.add_argument('--nbnds', type=int, default=1, help='Number of bands')

args = parser.parse_args()
vasp_file = args.vasp_file
qe_input_file = args.qe_input
E_cutoff = args.E_cutoff
conv_thr = args.conv_thr
nbnds = args.nbnds

# Read POSCAR (VASP format)
atoms = read(vasp_file)

# print("Elements:", atoms.get_chemical_symbols())
# print("Lattice vectors:\n", atoms.get_cell())
# print("Atomic positions:\n", atoms.get_positions())

# Define pseudopotentials mapping (element → filename without `.UPF` or `.pbe-n-kjpaw_psl.1.0.0.UPF`)
pseudopotentials = {
   "Si": "Si.upf"
}

input_data_dict = {
        "control": {
            "calculation": "scf",
            "prefix": "qe_calc",
            "pseudo_dir": "./",
            "verbosity": "high",
            "tstress": True,
            "tprnfor": True,
            "outdir": "./out"
        },
        "system": {
            "ecutwfc": E_cutoff,
            "nbnd": nbnds,
            "occupations": "smearing",
            "smearing": "gaussian",
            "degauss": 0.01
        },
        "electrons": {
            "conv_thr": conv_thr,
            'diagonalization': 'david',
            'diago_david_ndim': 4,
            'diago_full_acc': True
        }
}


# Write QE input
write(qe_input_file,
    atoms,
    format="espresso-in",
    pseudopotentials=pseudopotentials,
    input_data=input_data_dict)

