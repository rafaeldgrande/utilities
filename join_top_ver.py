
# ==============================================================================
#
#
#         This program join two lammps topologies in one!
#         You must edit the editable part of the code above. You must
#         inform the two topo files that you want to merge and if you want
#         to use the box boundaires looking for the highest and lowest
#         x, y and z coordinates of atoms in both files or use the limits
#         box boundaires of the input topo files (recommended if you are
#         working with crystal or any molecule that for some reason is
#         going through some of the box walls)
#
#         In lammps it's possible to read two different topo files
#         to perform one simulation. Although I've experienced some
#         problems using this feature, so this script is still useful :)
#
#         Ps.: -Be aware that this scrip doesn't check if there is
#         superposition among atoms from each file. I recommend to
#         minimize each part of your system before merging files
#              -I've wrote this script in a linux machine and I use
#         some linux commands in this routine. Maybe you should comment
#         the lines of the code that uses the "os" module (Use CTRL+F!)
#
#         Author : Rafael Del Grande - rafaeldgrande@gmail.com
#         Last update : 07/Jan/2016
#
#         7jan2016 - Bug in dihedrals section fixed (wrong list's indexes)
#        25Feb2016 - If there is not some dihedral, angle or bond don't
#        write the respective section
#        10Mar2016 - PEP8 style!
#
# ==============================================================================

# EDITABLE PART OF THE CODE

top1 = '2SiNP.top'
top2 = 'water.top'
output = '2SiNP_water.top'

new_box_booundaries = True   # If true we define the box boundaries of
                             # the output topo file by the atomic
                             # coordinates. If else then use the box
                             # boundaries of the input topo files.

# Modules

import string as string
import os

# Files

atoms_data = 'atoms.data'    # Auxiliary files
bonds_data = 'bonds.data'
angles_data = 'angles.data'
dihed_data = 'dihed.data'

# Variables

n_atoms1, n_atoms2 = 0, 0
n_mol1, n_mol2 = 0, 0
n_bonds1, n_bonds2 = 0, 0
n_angles1, n_angles2 = 0, 0
n_dihedrals1, n_dihedrals2 = 0, 0

atom_types1, atom_types2 = 0, 0
bond_types1, bond_types2 = 0, 0
anlge_types1, angle_types2 = 0, 0
dihed_types1, dihed_types2 = 0, 0

xlo, xhi, ylo, yhi, zlo, zhi = 10000, -10000, 10000, -10000, 10000, -10000

pair_coeffs = []
bond_coeffs = []
angle_coeffs = []
dihedral_coeffs = []
masses = []

# ==============================================================================
# ==============================================================================
#
# Reading first file
#
# ==============================================================================
# ==============================================================================

print 'Reading ', top1, 'file \n'

arq = open(top1)

j_pair_coeffs = 0
j_bond_coeffs = 0
j_angle_coeffs = 0
j_dihed_coeffs = 0
j_masses = 0

i_pair_coeffs = False
i_bond_coeffs = False
i_angle_coeffs = False
i_dihed_coeffs = False
i_masses = False

for line in arq:

    linha = string.split(line)

    if len(linha) >= 1:
        if linha[0] == 'Masses':
            i_masses = True

    if len(linha) >= 2:
        if linha[1] == 'atoms':            # ... atoms
            n_atoms1 = int(linha[0])
        elif linha[1] == 'bonds':          # ... bonds
            n_bonds1 = int(linha[0])
        elif linha[1] == 'angles':         # ... angles
            n_angles1 = int(linha[0])
        elif linha[1] == 'dihedrals':      # ... dihedrals
            n_dihedrals1 = int(linha[0])
        elif linha[0] == 'Pair':           # Pair Coeffs
            i_pair_coeffs = True
        elif linha[0] == 'Bond':           # Bond Coeffs
            i_bond_coeffs = True
        elif linha[0] == 'Angle':          # Angle Coeffs
            i_angle_coeffs = True
        elif linha[0] == 'Dihedral':       # Dihedral Coeffs
            i_dihed_coeffs = True

    if len(linha) >= 3:
        if linha[1] == 'atom':             # ... atom types
            atom_types1 = int(linha[0])
        elif linha[1] == 'bond':           # ... bond types
            bond_types1 = int(linha[0])
        elif linha[1] == 'angle':          # ... angle types
            angle_types1 = int(linha[0])
        elif linha[1] == 'dihedral':       # ... dihedral types
            dihed_types1 = int(linha[0])

    if new_box_booundaries is False:
        if len(linha) >= 4:
            if linha[3] == 'xhi':
                if float(linha[0]) < xlo:
                    xlo = float(linha[0])
                if float(linha[1]) > xhi:
                    xhi = float(linha[1])
            if linha[3] == 'yhi':
                if float(linha[0]) < ylo:
                    ylo = float(linha[0])
                if float(linha[1]) > yhi:
                    yhi = float(linha[1])
            if linha[3] == 'zhi':
                if float(linha[0]) < zlo:
                    zlo = float(linha[0])
                if float(linha[1]) > zhi:
                    zhi = float(linha[1])

    if len(linha) > 1 and i_masses is True:
        if linha[0] != 'Masses':
            j_masses += 1
            masses.append(linha)
            if j_masses == atom_types1:
                i_masses = False

    if len(linha) > 1 and i_pair_coeffs is True:
        if linha[1] != 'Coeffs':
            j_pair_coeffs += 1
            pair_coeffs.append(linha)
            if j_pair_coeffs == atom_types1:
                i_pair_coeffs = False

    if len(linha) > 1 and i_bond_coeffs is True:
        if linha[1] != 'Coeffs':
            j_bond_coeffs += 1
            bond_coeffs.append(linha)
            if j_bond_coeffs == bond_types1:
                i_bond_coeffs = False

    if len(linha) > 1 and i_angle_coeffs is True:
        if linha[1] != 'Coeffs':
            j_angle_coeffs += 1
            angle_coeffs.append(linha)
            if j_angle_coeffs == angle_types1:
                i_angle_coeffs = False

    if len(linha) > 1 and i_dihed_coeffs is True:
        if linha[1] != 'Coeffs':
            j_dihed_coeffs += 1
            dihedral_coeffs.append(linha)
            if j_dihed_coeffs == dihed_types1:
                i_dihed_coeffs = False

arq.close()

print n_atoms1, ' atoms'
print n_bonds1, ' bonds'
print n_angles1, ' angles'
print n_dihedrals1, ' dihedrals\n'

print atom_types1, ' atom types'
print bond_types1, ' bond types'
print angle_types1, ' angle types'
print dihed_types1, ' dihedral types\n\n'

# ==============================================================================
# ==============================================================================
#
# Reading second file
#
# ==============================================================================
# ==============================================================================

print 'Reading ', top2, 'file \n'

arq = open(top2)

i_pair_coeffs = False
i_bond_coeffs = False
i_angle_coeffs = False
i_dihed_coeffs = False
i_masses = False

j_pair_coeffs = 0
j_bond_coeffs = 0
j_angle_coeffs = 0
j_dihed_coeffs = 0
j_masses = 0

for line in arq:

    linha = string.split(line)

    if len(linha) >= 1:
        if linha[0] == 'Masses':
            i_masses = True

    if len(linha) >= 2:
        if linha[1] == 'atoms':
            n_atoms2 = int(linha[0])
        elif linha[1] == 'bonds':
            n_bonds2 = int(linha[0])
        elif linha[1] == 'angles':
            n_angles2 = int(linha[0])
        elif linha[1] == 'dihedrals':
            n_dihedrals2 = int(linha[0])
        elif linha[0] == 'Pair':
            i_pair_coeffs = True
        elif linha[0] == 'Bond':
            i_bond_coeffs = True
        elif linha[0] == 'Angle':
            i_angle_coeffs = True
        elif linha[0] == 'Dihedral':
            i_dihed_coeffs = True

    if len(linha) >= 3:
        if linha[1] == 'atom':
            atom_types2 = int(linha[0])
        elif linha[1] == 'bond':
            bond_types2 = int(linha[0])
        elif linha[1] == 'angle':
            angle_types2 = int(linha[0])
        elif linha[1] == 'dihedral':
            dihed_types2 = int(linha[0])

    if new_box_booundaries is False:
        if len(linha) >= 4:
            if linha[3] == 'xhi':
                if float(linha[0]) < xlo:
                    xlo = float(linha[0])
                if float(linha[1]) > xhi:
                    xhi = float(linha[1])
            if linha[3] == 'yhi':
                if float(linha[0]) < ylo:
                    ylo = float(linha[0])
                if float(linha[1]) > yhi:
                    yhi = float(linha[1])
            if linha[3] == 'zhi':
                if float(linha[0]) < zlo:
                    zlo = float(linha[0])
                if float(linha[1]) > zhi:
                    zhi = float(linha[1])

    if len(linha) > 1 and i_masses is True:
        if linha[0] != 'Masses':
            j_masses += 1
            linha[0] = str(int(linha[0]) + atom_types1)
            masses.append(linha)
            if j_masses == atom_types2:
                i_masses = False

    if len(linha) > 1 and i_pair_coeffs is True:
        if linha[-1] != 'Coeffs':
            j_pair_coeffs += 1
            linha[0] = str(int(linha[0]) + atom_types1)
            pair_coeffs.append(linha)
            if j_pair_coeffs == atom_types2:
                i_pair_coeffs = False

    if len(linha) > 1 and i_bond_coeffs is True:
        if linha[-1] != 'Coeffs':
            j_bond_coeffs += 1
            linha[0] = str(int(linha[0]) + bond_types1)
            bond_coeffs.append(linha)
            if j_bond_coeffs == bond_types2:
                i_bond_coeffs = False

    if len(linha) > 1 and i_angle_coeffs is True:
        if linha[-1] != 'Coeffs':
            j_angle_coeffs += 1
            linha[0] = str(int(linha[0]) + angle_types1)
            angle_coeffs.append(linha)
            if j_angle_coeffs == angle_types2:
                i_angle_coeffs = False

    if len(linha) > 1 and i_dihed_coeffs is True:
        if linha[-1] != 'Coeffs':
            j_dihed_coeffs += 1
            linha[0] = str(int(linha[0]) + dihed_types1)
            dihedral_coeffs.append(linha)
            if j_dihed_coeffs == dihed_types2:
                i_dihed_coeffs = False

arq.close()

print n_atoms2, ' atoms'
print n_bonds2, ' bonds'
print n_angles2, ' angles'
print n_dihedrals2, ' dihedrals\n'

print atom_types2, ' atom types'
print bond_types2, ' bond types'
print angle_types2, ' angle types'
print dihed_types2, ' dihedral types\n\n'


# ==============================================================================
# ==============================================================================
#
# Intermidiarie files
#
# ==============================================================================
# ==============================================================================


Atoms = open(atoms_data, 'w')
Bonds = open(bonds_data, 'w')
Angles = open(angles_data, 'w')
Dihedrals = open(dihed_data, 'w')


# ==============================================================================
# ==============================================================================
#
# Reading first file again
#
# ==============================================================================
# ==============================================================================

print 'Saving data from ', top1, ' in intermediarie files'

i_atoms = False
i_bonds = False
i_angles = False
i_dihed = False

j_atoms = 0
j_bonds = 0
j_angles = 0
j_dihed = 0

arq = open(top1)

for line in arq:

    linha = string.split(line)

    if len(linha) >= 1:
        if linha[0] == 'Atoms':
            i_atoms = True
        elif linha[0] == 'Bonds':
            i_bonds = True
        elif linha[0] == 'Angles':
            i_angles = True
        elif linha[0] == 'Dihedrals':
            i_dihed = True
        elif linha[0] == 'Masses':
            i_masses = True

    # atom_id mol_id atom_type charge x y z flagx flagy flagz #comment
    if len(linha) > 1 and i_atoms is True:

        if linha[0] != 'Atoms':

            j_atoms += 1

            if int(linha[1]) > n_mol1:
                n_mol1 = int(linha[1])

            for i in linha:  # Writing atoms data
                Atoms.write(i + '     ')
            Atoms.write('\n')

            # Use atoms coordinates to set new box boundaries?
            if new_box_booundaries is True:
                if xlo > float(linha[4]):
                    xlo = float(linha[4])
                elif xhi < float(linha[4]):
                    xhi = float(linha[4])
                if ylo > float(linha[5]):
                    ylo = float(linha[5])
                elif yhi < float(linha[5]):
                    yhi = float(linha[5])
                if zlo > float(linha[6]):
                    zlo = float(linha[6])
                elif zhi < float(linha[6]):
                    zhi = float(linha[6])

            if j_atoms == n_atoms1:
                i_atoms = False

    # bond_id bond_type atom1 atom2
    if len(linha) > 1 and i_bonds is True:

        if linha[0] != 'Bonds':

            j_bonds += 1

            for i in linha:
                Bonds.write(i + '     ')
            Bonds.write('\n')

            if j_bonds == n_bonds1:
                i_bonds = False

    # Angle_id angle_type atom1 atom2 atom3
    if len(linha) > 1 and i_angles is True:

        if linha[0] != 'Angles':

            j_angles += 1

            for i in linha:
                Angles.write(i + '        ')
            Angles.write('\n')

            if j_angles == n_angles1:
                i_angles = False

    # Dihed_id dihed_type atom1 atom2 atom3 atom4
    if len(linha) > 1 and i_dihed is True:

        if linha[0] != 'Dihedrals':

            j_dihed += 1

            for i in linha:
                Dihedrals.write(i + '     ')
            Dihedrals.write('\n')

            if j_dihed == n_dihedrals1:
                i_dihed = False

arq.close()

# ==============================================================================
# ==============================================================================
#
# Reading second file again
#
# ==============================================================================
# ==============================================================================

print 'Saving data from ', top2, ' in intermediarie files'

i_atoms = False
i_bonds = False
i_angles = False
i_dihed = False

j_atoms = 0
j_bonds = 0
j_angles = 0
j_dihed = 0

arq = open(top2)

for line in arq:

    linha = string.split(line)

    if len(linha) >= 1:
        if linha[0] == 'Atoms':
            i_atoms = True
        elif linha[0] == 'Bonds':
            i_bonds = True
        elif linha[0] == 'Angles':
            i_angles = True
        elif linha[0] == 'Dihedrals':
            i_dihed = True
        elif linha[0] == 'Masses':
            i_masses = True

    if len(linha) >= 2:
        if linha[0] == 'Pair':
            i_pair_coeffs = True
        elif linha[0] == 'Bond':
            i_bond_coeffs = True
        elif linha[0] == 'Angle':
            i_angle_coeffs = True
        elif linha[0] == 'Dihedral':
            i_dihed_coeffs = True

    # atom_id mol_id atom_type charge x y z flagx flagy flagz
    if len(linha) > 1 and i_atoms is True:

        if linha[0] != 'Atoms':

            j_atoms += 1

            linha[0] = str(int(linha[0]) + n_atoms1)
            linha[1] = str(int(linha[1]) + n_mol1)
            linha[2] = str(int(linha[2]) + atom_types1)

            for i in linha:
                Atoms.write(i + '    ')
            Atoms.write('\n')

            # Use atoms coordinates to set new box boundaries?
            if new_box_booundaries is True:

                if xlo > float(linha[4]):
                    xlo = float(linha[4])
                elif xhi < float(linha[4]):
                    xhi = float(linha[4])

                if ylo > float(linha[5]):
                    ylo = float(linha[5])
                elif yhi < float(linha[5]):
                    yhi = float(linha[5])

                if zlo > float(linha[6]):
                    zlo = float(linha[6])
                elif zhi < float(linha[6]):
                    zhi = float(linha[6])

            if j_atoms == n_atoms2:
                i_atoms = False

    # bond_id bond_type atom1 atom2
    if len(linha) > 1 and i_bonds is True:

        if linha[0] != 'Bonds':

            j_bonds += 1

            linha[0] = str(int(linha[0]) + n_bonds1)
            linha[1] = str(int(linha[1]) + bond_types1)
            linha[2] = str(int(linha[2]) + n_atoms1)
            linha[3] = str(int(linha[3]) + n_atoms1)

            for i in linha:
                Bonds.write(i + '     ')
            Bonds.write('\n')

            if j_bonds == n_bonds2:
                i_bonds = False

    # Angle_id angle_type atom1 atom2 atom3
    if len(linha) > 1 and i_angles is True:

        if linha[0] != 'Angles':

            j_angles += 1

            linha[0] = str(int(linha[0]) + n_angles1)
            linha[1] = str(int(linha[1]) + angle_types1)
            linha[2] = str(int(linha[2]) + n_atoms1)
            linha[3] = str(int(linha[3]) + n_atoms1)
            linha[4] = str(int(linha[4]) + n_atoms1)

            for i in linha:
                Angles.write(i + '        ')
            Angles.write('\n')

            if j_angles == n_angles2:
                i_angles = False

    # Dihed_id dihed_type atom1 atom2 atom3 atom4
    if len(linha) > 1 and i_dihed is True:

        if linha[0] != 'Dihedrals':

            j_dihed += 1

            linha[0] = str(int(linha[0]) + n_dihedrals1)
            linha[1] = str(int(linha[1]) + dihed_types1)
            linha[2] = str(int(linha[2]) + n_atoms1)
            linha[3] = str(int(linha[3]) + n_atoms1)
            linha[4] = str(int(linha[4]) + n_atoms1)
            linha[5] = str(int(linha[5]) + n_atoms1)

        for i in linha:
            Dihedrals.write(i + '     ')
        Dihedrals.write('\n')

        if j_dihed == n_dihedrals2:
            i_dihed = False

arq.close()


# ==============================================================================
# ==============================================================================
#
# Closing intermediaries files
#
# ==============================================================================
# ==============================================================================

Atoms.close()
Bonds.close()
Angles.close()
Dihedrals.close()

# ==============================================================================
# ==============================================================================
#
# Writing final top file!!!
#
# ==============================================================================
# ==============================================================================

N_atoms = n_atoms1 + n_atoms2
N_bonds = n_bonds1 + n_bonds2
N_angles = n_angles1 + n_angles2
N_dihedrals = n_dihedrals1 + n_dihedrals2

Atom_types = atom_types1 + atom_types2
Bond_types = bond_types1 + bond_types2
Angle_types = angle_types1 + angle_types2
Dihedral_types = dihed_types1 + dihed_types2

print 'New box boundaries :\n'
print 'xlo xhi = ', xlo, xhi
print 'ylo yhi = ', ylo, yhi
print 'zlo zhi = ', zlo, zhi

print '\nWriting final top final file\n'

top_final = open(output, 'w')

Atoms = open(atoms_data)
Bonds = open(bonds_data)
Angles = open(angles_data)
Dihedrals = open(dihed_data)

top_final.write('\n\n')

top_final.write(str(N_atoms) + '  atoms\n')
top_final.write(str(N_bonds) + '  bonds\n')
top_final.write(str(N_angles) + ' angles\n')
top_final.write(str(N_dihedrals) + '  dihedrals\n')

top_final.write('\n')

top_final.write(str(Atom_types) + '   atom types\n')
top_final.write(str(Bond_types) + '   bond types\n')
top_final.write(str(Angle_types) + '   angle types\n')
top_final.write(str(Dihedral_types) + '   dihedral types\n')

top_final.write('\n')

top_final.write(str(xlo) + '  ' + str(xhi) + '    xlo xhi\n')
top_final.write(str(ylo) + '  ' + str(yhi) + '    ylo yhi\n')
top_final.write(str(zlo) + '  ' + str(zhi) + '    zlo zhi\n')

top_final.write('\n')

top_final.write(' Masses\n\n')

for item in masses:
    for i in item:
        top_final.write(i + '	')
    top_final.write('\n')

top_final.write('\n\n')

top_final.write(' Pair Coeffs\n\n')

for item in pair_coeffs:
    for i in item:
        top_final.write(i + ' ')
    top_final.write('\n')

top_final.write('\n\n')

top_final.write(' Bond Coeffs\n\n')

if N_bonds > 0:
    for item in bond_coeffs:
        for i in item:
            top_final.write(i + ' ')
        top_final.write('\n')

    top_final.write('\n\n')

if N_angles > 0:
    top_final.write(' Angle Coeffs\n\n')

    for item in angle_coeffs:
        for i in item:
            top_final.write(i + ' ')
        top_final.write('\n')

    top_final.write('\n\n')

if N_dihedrals > 0:
    top_final.write(' Dihedral Coeffs\n\n')

    for item in dihedral_coeffs:
        for i in item:
            top_final.write(i + ' ')
        top_final.write('\n')

    top_final.write('\n\n')

print 'Writing atoms section\n'

top_final.write('\n\n Atoms\n\n')

for line in Atoms:
    top_final.write(line)

if N_bonds > 0:
    print 'Writing bonds section\n'

    top_final.write('\n\n Bonds\n\n')

    for line in Bonds:
        top_final.write(line)

if N_angles > 0:
    print 'Writing angles section\n'

    top_final.write('\n\n Angles\n\n')

    for line in Angles:
        top_final.write(line)

if N_dihedrals > 0:
    print 'Writing dihedrals section\n'

    top_final.write('\n\n Dihedrals\n\n')

    for line in Dihedrals:
        top_final.write(line)


top_final.close()

print 'Deleting auxiliary files'

os.system('rm atoms.data bonds.data angles.data dihed.data')

print '\n\nDone!'
