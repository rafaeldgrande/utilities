#==================================================================================
#
#
#          Script that cuts a box of silica and creates a spherical nanoparticle
#
#          Input : .xyz file
#
#          Output : nanopart_Si.xyz
#
#          updates: 11Fev15 - Now it's possible to define the OH/H ratio
#                             over the surface
#
#==================================================================================


# Modules

import string as string
import numpy as np
import time as time

Ti = time.time()   # Counting how much time the program gets to do stuff!

# Classes


class atom:
	x, y, z = float, float, float
	element = str
	num_coord = 0

# Files

input_xyz = 'amorphous_quartz.xyz'
output_xyz = 'nanoparticle_Si_10nm.xyz'

# Variables

old_atoms = []            # List containing atoms of the original file
new_atoms = []            # List containing atoms of the new file

OH_groups = 0             # How many OH bounded to Si
H_groups = 0              # How many H bounded to Si

x_min_old, x_max_old = float, float
y_min_old, y_max_old = float, float       # Old box boundaries
z_min_old, z_max_old = float, float

X0, Y0, Z0 = float, float, float          # Old box center

R = 20.0                                  # Nanoparticle radius (angstrons)

delta = 2.0                               # Atoms at a distance R < r < R - delta
                                          # probably will lose some neighbors during the nanoparticle
                                          # creation. The script will save those atoms' indexes, then
                                          # we can complete their valences

hidrogenation = 0.5                       # How much OH there are over the surface
                                          # 0 = no OH and full of H,
					  # 1 = full of OH and no H

Si_coord = 4                              # Atoms's coordination
O_coord = 2

Si_O = 1.8                                # Bond lenght
O_H = 1.0

# Functions

def dist (atom1, atom2) :                # Calculates the distance between two atoms
	return np.sqrt((atom1.x-atom2.x)**2 + (atom1.y-atom2.y)**2 + (atom1.z-atom2.z)**2)

def spherical_constrain (x, y, z, x0, y0, z0, R) :     # Nanoparticle creator!
	if ((x-x0)**2 + (y-y0)**2 + (z-z0)**2) < R**2:
		return True
	else :
		return False

#--------------------------Reading input file -----------------------------------------------------

print 'Reading input file\n'

input_file = open(input_xyz)

for line in input_file:
	lin = string.split(line)
	if len(lin) == 4:  # Element, x, y, z
		A = atom()
		A.element = lin[0]
		A.x, A.y, A.z = float(lin[1]), float(lin[2]), float(lin[3])
		old_atoms.append(A)

input_file.close()

x_min_old, x_max_old = old_atoms[0].x, old_atoms[0].x
y_min_old, y_max_old = old_atoms[0].y, old_atoms[0].y
z_min_old, z_max_old = old_atoms[0].z, old_atoms[0].z

for item in old_atoms:         # Looking for the old box boundaries
	if item.x > x_max_old:
		x_max_old = item.x
	elif item.x < x_min_old:
		x_min_old = item.x
	if item.y > y_max_old:
		y_max_old = item.y
	elif item.y < y_min_old:
		y_min_old = item.y
	if item.z > z_max_old:
		z_max_old = item.z
	elif item.z < z_min_old:
		z_min_old = item.z

print 'Old box boundaries'
print 'Xmin ', x_min_old, ' Xmax ', x_max_old
print 'Ymin ', y_min_old, ' Ymax ', y_max_old
print 'Zmin ', z_min_old, ' Zmax ', z_max_old, '\n'

x0 = (x_max_old + x_min_old	)/2
y0 = (y_max_old + y_min_old	)/2
z0 = (z_max_old + z_min_old)/2

print 'Box center:', x0, ', ', y0, ', ', z0

center = atom()    # Virtual atom representing the box center
center.x, center.y, center.z = x0, y0, z0

#--------------------------Saving oxigen atoms and cutting the box!!!-------------------------------

print 'Checking atoms valences that are near the surfaces\n'

Si_out_shell = []
O_out_shell = []

Si_in_shell = []
O_in_shell = []

for item in old_atoms :
	if item.element == 'O' :
		if spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R+delta) == True and spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R) == False:
			O_out_shell.append(item)  # R < r < R - delta
		elif spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R) == True and spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R - delta) == False:
			O_in_shell.append(item)   # R - delta < r < 0
	elif item.element == 'Si' :
		if spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R+delta) == True and spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R) == False:
			Si_out_shell.append(item) # R < r < R - delta
		elif spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R) == True and spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R - delta) == False:
			Si_in_shell.append(item)  # R - delta < r < 0

	if spherical_constrain (item.x, item.y, item.z, x0, y0, z0, R) == True :
		new_atoms.append(item)

#-------------------------Which atoms will be saved?-----------------------------------------------------------

print 'Completing the Si atoms valence by adding OH or H groups to them\n'

for O in O_out_shell :
	for Si in Si_in_shell :
		if dist(Si, O) <= Si_O :
			if new_atoms.count(O) == 0 :

				random_var = np.random.rand()

				if random_var <= hidrogenation:
					OH_groups += 1
					new_atoms.append(O)
					H = atom()                        #Adding H atom in the radial direction
					delx = O_H*(O.x - center.x)/dist(O, center)
					dely = O_H*(O.y - center.y)/dist(O, center)
					delz = O_H*(O.z - center.z)/dist(O, center)
					H.element, H.x, H.y, H.z = 'H', O.x + delx, O.y + dely, O.z + delz
					new_atoms.append(H)

				else:
					H_groups += 1
					O.element = 'H'
					new_atoms.append(O)


#---------------------------Adding missing H's-----------------------------------------------------------------

print 'Completing the O atoms valence by adding H atoms to them or transforming them in H atoms\n'

dump_O = []

for O in O_in_shell :
	for item in new_atoms :
		if item.element == 'Si' :
			if 0 < dist(item, O) < Si_O :
				O.num_coord += 1
	if O.num_coord == 1 :
		random_var = np.random.rand()
		if random_var <= hidrogenation:
			OH_groups += 1
			H = atom()                        #Adding H atom in the radial direction
			delx = O_H*(O.x - center.x)/dist(O, center)
			dely = O_H*(O.y - center.y)/dist(O, center)
			delz = O_H*(O.z - center.z)/dist(O, center)
			H.element, H.x, H.y, H.z = 'H', O.x + delx, O.y + dely, O.z + delz
			new_atoms.append(H)
		else:
			H_groups += 1
			O.element = 'H'
	elif O.num_coord == 0 :
		dump_O.append(O)

for atom in new_atoms :
	if dump_O.count(atom) != 0 :
		del new_atoms[new_atoms.index(atom)]

#--------------------------Translating all atoms!--------------------------------------------------------------

print 'Putting the nanoparticle center in the origin'

for item in new_atoms :
	item.x -= x0
	item.y -= y0
	item.z += z0

#--------------------------Writing output file----------------------------------------------------------------

print 'Writing the output file\n'

output_file = open(output_xyz, 'w')

output_file.write(str(len(new_atoms))+'\n\n')

for item in new_atoms :
	output_file.write(item.element+'	'+str(item.x)+'	'+str(item.y)+'	'+str(item.z)+'\n')

output_file.close()

print 'Finished!\n'

Tf = time.time()

print 'There are ', H_groups, ' H groups = ', 100.0*H_groups/(OH_groups+H_groups),' %'
print 'There are ', OH_groups, ' OH groups =', 100.0*OH_groups/(OH_groups+H_groups), '%'

print 'I got ', float(Tf-Ti)/60.0,' minutes to do all those stuff'
