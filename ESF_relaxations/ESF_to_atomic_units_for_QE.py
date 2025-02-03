
import numpy as np
import sys

# Check if the user provided the forces_file as a command-line argument
if len(sys.argv) < 2:
    print("Error: Please provide the forces file path as a command-line argument.")
    sys.exit(1)

# Get forces_file from command-line argument
forces_file = sys.argv[1]

flavor = 2

bohr2A = 0.529177
ry2eV = 13.6057039763

eV2ry = 1 / ry2eV
A2bohr = 1 / bohr2A

def read_excited_forces(excited_state_forces_file, flavor):
    # flavor = 1 -> RPA_diag
    # flavor = 2 -> RPA_diag_offiag
    # flavor = 3 -> RPA_diag_Kernel

    data = np.loadtxt(excited_state_forces_file, usecols=flavor+1)
    return data

# Read data and convert units
data = read_excited_forces(forces_file, flavor) * eV2ry / A2bohr

# Number of atoms (each atom has 3 components: x, y, z)
Natoms = int(data.shape[0] / 3)

# Output atomic forces
print("ATOMIC_FORCES")
for iatom in range(Natoms):
    print(f"{iatom+1}   {data[3*iatom]:.8f}   {data[3*iatom+1]:.8f}   {data[3*iatom+2]:.8f} ")

