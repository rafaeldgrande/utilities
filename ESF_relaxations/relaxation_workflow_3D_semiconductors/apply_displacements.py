
import numpy as np

atoms_orig_file = 'ATOMS'
atoms_new_file = 'DISP_ATOMS'
displacement_file = '10-displacements/displacements_Newton_method.dat'

# read atoms original file
print(f'Reading atoms_orig_file = {atoms_orig_file}')

symbols, positions0 = [], [] 
arq = open(atoms_orig_file)
for line in arq:
    line_split = line.split()
    if len(line_split) >= 4:
        symbols.append(line_split[0])
        x, y, z = float(line_split[1]), float(line_split[2]), float(line_split[3])
        positions0.append([x, y, z])
arq.close()

positions0 = np.array(positions0)
print('positions0.shape = ', positions0.shape)

# read displacements
print(f"Reading displacements file {displacement_file}")
disp = np.loadtxt(displacement_file)[:, 1:]
print('disp.shape = ', disp.shape)
# print(disp[:10, :])

# apply displacements and write new file
positions_f = positions0 + disp

print(f'Writing final positions file {atoms_new_file}')

arq = open(atoms_new_file, 'w')
arq.write('ATOMIC_POSITIONS angstrom\n')
for iatom, pos in enumerate(positions_f):
    arq.write(f"{symbols[iatom]}   {pos[0]:.8f}   {pos[1]:.8f}   {pos[2]:.8f}\n")
arq.close()

print('Finished!')





