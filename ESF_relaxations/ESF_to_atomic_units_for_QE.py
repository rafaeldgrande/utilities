
import numpy as np
import sys

print("Usage: python ESF_to_atomic_units_for_QE.py forces_file x")
print("Read the forces_file from ESF code and convert them to atomic units for Quantum Espresso input files")
print("x is the exciton concentration in exciton / unit cell")

# Check if the user provided the forces_file as a command-line argument
if len(sys.argv) != 3:
    print("Error: Usage: python ESF_to_atomic_units_for_QE.py forces_file x")
    sys.exit(1)

# Get forces_file from command-line argument
forces_file = sys.argv[1]

# exciton concentration
x = float(sys.argv[2])

flavor = 2

bohr2A = 0.529177
ry2eV = 13.6057039763

eV2ry = 1 / ry2eV
A2bohr = 1 / bohr2A

def read_excited_forces(excited_state_forces_file, flavor):
    # flavor = 1 -> RPA_diag
    # flavor = 2 -> RPA_diag_offiag
    # flavor = 3 -> RPA_diag_Kernel

    # data = np.loadtxt(excited_state_forces_file, usecols=flavor+1)
    data = np.loadtxt(excited_state_forces_file, usecols=flavor+1, dtype=str)  # Read as string
    data = np.array([complex(x).real for x in data])  # Convert to complex and extract real part

    return data

# Read data and convert units
data = x * read_excited_forces(forces_file, flavor) * eV2ry / A2bohr

# Number of atoms (each atom has 3 components: x, y, z)
Natoms = int(data.shape[0] / 3)

# Output atomic forces
print("ATOMIC_FORCES")
for iatom in range(Natoms):
    print(f"{iatom+1}   {data[3*iatom]:.8f}   {data[3*iatom+1]:.8f}   {data[3*iatom+2]:.8f} ")

