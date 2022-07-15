#==================================================================================
#
#
#          Script that cuts a box of some solid material in some topology
#          and creates a slab of silica
#           
#          Input : .xyz file
#
#          Output : out.xyz
#
#
#==================================================================================

#Modules

import string as string
import numpy as np
import time as time

Ti = time.time()   #Counting how much time the program gets to do stuff!

#Classes

class atom :
    x, y, z = float, float, float
    element = str
    num_coord = 0
    OH_groups = 0   #Those two last properties are responsible for the analysis of the Si(OH)n proportions
    H_groups = 0

#Files

input_xyz = 'amorphous_quartz.xyz'
output_xyz = 'out.xyz'

#Variables

old_atoms = []            #List containing atoms of the original file
new_atoms = []            #List containing atoms of the new file

x_min_old, x_max_old = float, float     
y_min_old, y_max_old = float, float       #Old box boundaries 
z_min_old, z_max_old = float, float

x_min_new, x_max_new = float, float     
y_min_new, y_max_new = float, float       #New box boundaries 
z_min_new, z_max_new = float, float

nSiOH_1 = 0
nSiOH_2 = 0
nSiOH_3 = 0

lz = 40                                   #New box height

Si_coord = 4
O_coord = 2

Si_O = 1.8 
O_H = 1.0

#Functions

def dist (atom1, atom2) :
    return np.sqrt((atom1.x-atom2.x)**2 + (atom1.y-atom2.y)**2 + (atom1.z-atom2.z)**2)
    
#--------------------------Reading input file -----------------------------------------------------

print 'Reading input file\n'

input_file = open(input_xyz)

for line in input_file :
    lin = string.split(line)
    if len(lin) == 4 :  #Element, x, y, z
        A = atom()
        A.element = lin[0]
        A.x, A.y, A.z = float(lin[1]), float(lin[2]), float(lin[3])
        old_atoms.append(A)
        
input_file.close()

x_min_old, x_max_old = old_atoms[0].x, old_atoms[0].x
y_min_old, y_max_old = old_atoms[0].y, old_atoms[0].y
z_min_old, z_max_old = old_atoms[0].z, old_atoms[0].z

for item in old_atoms :
    if item.x > x_max_old :
        x_max_old = item.x
    elif item.x < x_min_old :
        x_min_old = item.x
    if item.y > y_max_old :
        y_max_old = item.y
    elif item.y < y_min_old :
        y_min_old = item.y
    if item.z > z_max_old :
        z_max_old = item.z
    elif item.z < z_min_old :
        z_min_old = item.z

print 'Old box boundaries'
print 'Xmin ', x_min_old, ' Xmax ', x_max_old 
print 'Ymin ', y_min_old, ' Ymax ', y_max_old
print 'Zmin ', z_min_old, ' Zmax ', z_max_old, '\n\n'

#--------------------------Saving oxigen atoms from death and cutting the box!!!-------------------------------

print 'Checking atoms valences that are near the surfaces\n'

O_upper_surf_out = []         #Oxigen atoms out of the box -> maybe they will be appended to complete the Si valence
O_lower_surf_out = []

Si_upper_surf_out = []        #Si atoms out of the box -> maybe they will be replaced by H atomsto complete the O valence
Si_lower_surf_out = []

Si_upper_surf_in = []         #Si atoms near the surface and inside the box. Their valence must be checked
Si_lower_surf_in = []

O_upper_surf_in = []          #O atoms near the surface and inside the box. Their valence must be checked
O_lower_surf_in = []

for item in old_atoms :
    if item.element == 'O' :
        if lz/2.0 < item.z <= lz/2.0 + 2.0 :
            O_upper_surf_out.append(item)
        elif lz/2.0 < -item.z <= lz/2.0 + 2.0 :
            O_lower_surf_out.append(item)
        elif lz/2.0 - 1.5 < item.z < lz/2.0 :
            O_upper_surf_in.append(item)
        elif (lz/2.0 - 1.5) <= -item.z <= lz/2.0  :
            O_lower_surf_in.append(item)
    if item.element == 'Si' :
        if lz/2.0 < item.z <= lz/2.0 + 2.0 :
            Si_upper_surf_out.append(item)
        elif lz/2.0 < -item.z <= lz/2.0 + 2.0 :
            Si_lower_surf_out.append(item)
        elif (lz/2.0 - 3.5) <= item.z <= lz/2.0 :
            Si_upper_surf_in.append(item)
        elif (lz/2.0 - 3.5) <= -item.z <= lz/2.0  :
            Si_lower_surf_in.append(item)
        
    if -lz/2 < item.z < lz/2 :
        new_atoms.append(item)
        
        
#-------------------------Which atoms will be saved?-----------------------------------------------------------

print 'Completing the Si atoms valence by adding OH groups to them\n'

for O in O_upper_surf_out :
    for Si in Si_upper_surf_in :
        if 0 < dist(Si, O) < Si_O :
            if new_atoms.count(O) == 0 :
                new_atoms.append(O)
                H = atom()
                H.element, H.x, H.y, H.z = 'H', O.x, O.y, O.z + O_H
                new_atoms.append(H)
                Si.OH_groups += 1
            
for O in O_lower_surf_out :
    for Si in Si_lower_surf_in :
        if 0 < dist(Si, O) < Si_O :
            if new_atoms.count(O) == 0 :
                new_atoms.append(O)
                H = atom()
                H.element, H.x, H.y, H.z = 'H', O.x, O.y, O.z - O_H
                new_atoms.append(H)
                Si.OH_groups += 1

#---------------------------Adding missing H's-----------------------------------------------------------------

print 'Completing the O atoms valence by adding H atoms to them\n'

dump_O = []

g = O_H/Si_O

for O in O_upper_surf_in :
    for Si in Si_upper_surf_in :
        if 0 < dist(Si, O) < Si_O :
            O.num_coord += 1
            temp = Si_upper_surf_in.index(Si)
    if O.num_coord == 1 :
        H = atom()
        H.element, H.x, H.y, H.z = 'H', O.x, O.y, O.z + O_H
        new_atoms.append(H)
        Si_upper_surf_in[temp].OH_groups += 1
    elif O.num_coord == 0 :
        dump_O.append(O)

for O in O_lower_surf_in :
    for Si in Si_lower_surf_in :
        if 0 < dist(Si, O) < Si_O :
            O.num_coord += 1   
    if O.num_coord == 1 :
        H = atom()
        H.element, H.x, H.y, H.z = 'H', O.x, O.y, O.z - O_H
        new_atoms.append(H)
    elif O.num_coord == 0 :
        dump_O.append(O)



for atom in new_atoms :
    if dump_O.count(atom) != 0 :
        del new_atoms[new_atoms.index(atom)]
        
#--------------------------Translating all atoms!--------------------------------------------------------------

for item in new_atoms :
    item.x -= x_min_old
    item.y -= y_min_old
    item.z -= lz/2 + 2.0
    
#--------------------------Writing output file----------------------------------------------------------------

print 'Writing the output file\n'

output_file = open(output_xyz, 'w')

output_file.write(str(len(new_atoms))+'\n\n')

for item in new_atoms :
    output_file.write(item.element+'    '+str(item.x)+'    '+str(item.y)+'    '+str(item.z)+'\n')

output_file.close()

print 'Finished!\n'

Tf = time.time()

print 'I got ', float(Tf-Ti)/60.0,' minutes to do all those stuff'

for i in range (len(Si_upper_surf_in)) :
    if Si_upper_surf_in[i].OH_groups == 1 :
        nSiOH_1 += 1
    elif Si_upper_surf_in[i].OH_groups == 2 :
        nSiOH_2 += 1
    elif Si_upper_surf_in[i].OH_groups == 3 :
        nSiOH_3 += 1
        
for i in range (len(Si_lower_surf_in)) :
    if Si_lower_surf_in[i].OH_groups == 1 :
        nSiOH_1 += 1
    elif Si_lower_surf_in[i].OH_groups == 2 :
        nSiOH_2 += 1
    elif Si_lower_surf_in[i].OH_groups == 3 :
        nSiOH_3 += 1
        
print 'Number of SiOH =', nSiOH_1
print 'Number of SiOH2 =', nSiOH_2
print 'Number of SiOH3 =', nSiOH_3




