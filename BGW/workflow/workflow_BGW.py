
'''

BGW workflow for semiconductors

-> following workflow for typical semiconductors

1-scf 
2-wfn
3-wfnq
4-wfn_co
5-wfn_fi
6-wfnq_fi
7-epsilon
8-sigma
9-kernel
10.1-absorption_momentum
10.2-absorption_velocity

'''

import os 
import numpy as np
from workflow_BGW_config import *
from workflow_BGW_mod import *

# First, assign variables to dictionary
prefix = configurations['prefix']
pseudo_dir = configurations['pseudo_dir']
qe_atoms_file = configurations['qe_atoms_file']
qe_cell_file = configurations['qe_cell_file']
kgrid_coarse = configurations['kgrid_coarse']
qshift_wfnq = configurations['qshift_wfnq']
kgrid_fine = configurations['kgrid_fine']
kpoints_file_scf = configurations['kpoints_file_scf']
kpoints_file_wfn = configurations['kpoints_file_wfn']                                                                                                                                                                                                                                                                              
kpoints_file_wfnq = configurations['kpoints_file_wfnq']
kpoints_file_wfn_fi = configurations['kpoints_file_wfn_fi']
kpoints_file_wfnq_fi = configurations['kpoints_file_wfnq_fi']
nval = configurations['nval']
ncond_wfn = configurations['ncond_wfn']
ncond_wfnq = configurations['ncond_wfnq']
ncond_wfnco = configurations['ncond_wfnco']
ncond_wfnfi = configurations['ncond_wfnfi']
nmin_gw = configurations['nmin_gw']
nmax_gw = configurations['nmax_gw']
nval_kernel = configurations['nval_kernel']
ncond_kernel = configurations['ncond_kernel']
nval_bse = configurations['nval_bse']
ncond_bse = configurations['ncond_bse']
dft_ekin_cutoff = configurations['dft_ekin_cutoff']
qshift_bse = configurations['qshift_bse']
epsilon_cutoff = configurations['epsilon_cutoff']
truncation_scheme = configurations['truncation_scheme']

# getting number of atoms
Ntypes, Natoms = 0, 0
arq = open(qe_atoms_file)
for line in arq:
    line_split = line.split()
    if len(line_split) == 3:
        Ntypes += 1
    elif len(line_split) == 4:
        Natoms += 1

                                        
# writing QE inputs

def GWBSE_write_files():
    
    # TODO: change for the option to use symmetry using kgrid.x from BGW
    print('Generating k grids - not using symmetry')
    generate_regular_kgrid(kpoints_file_scf, kgrid_fine, [0,0,0])
    generate_regular_kgrid(kpoints_file_wfn, kgrid_coarse, [0,0,0])
    generate_regular_kgrid(kpoints_file_wfnq, kgrid_coarse, qshift_wfnq)
    generate_regular_kgrid(kpoints_file_wfn_fi, kgrid_fine, [0,0,0])
    generate_regular_kgrid(kpoints_file_wfnq_fi, kgrid_fine, qshift_bse)
    print('')
    
    print('Creating directories for DFT calculations')
    os.system('mkdir -p 1-scf 2-wfn 3-wfnq 4-wfn_co 5-wfn_fi 6-wfnq_fi')

    for WFN_DIR in ['2-wfn', '3-wfnq', '4-wfn_co', '5-wfn_fi', '6-wfnq_fi']:
        print(f"Creating links for scf data. Now at dir {WFN_DIR}")
        os.system(f'mkdir -p {WFN_DIR}/{prefix}.save')
        os.system(f'ln -sf ../../1-scf/{prefix}.save/charge-density.hdf5 {WFN_DIR}/{prefix}.save/charge-density.hdf5')
        os.system(f'ln -sf ../../1-scf/{prefix}.save/charge-density.dat {WFN_DIR}/{prefix}.save/charge-density.dat')
        os.system(f'ln -sf ../../1-scf/{prefix}.save/data-file-schema.xml {WFN_DIR}/{prefix}.save/data-file-schema.xml')

    print('')
    print('Writing DFT inputs')

    # 1-scf
    nbnds_temp = nval + ncond_bse + 1
    qshift = [0,0,0]
    qgrid = kgrid_fine
    dir_temp = '1-scf'
    write_QE_input(f'{dir_temp}/scf.in', nbnds_temp, prefix, Natoms, Ntypes, pseudo_dir, 'scf', kpoints_file_scf, dft_ekin_cutoff, qe_cell_file, qe_atoms_file)
    write_pw2bgw_input(f'{dir_temp}/pw2bgw.in', nbnds_temp, prefix, qshift, qgrid)
    write_ph_input(f'{dir_temp}/ph.in', prefix)
    
    # 2-wfn
    nbnds_temp = nval + ncond_wfn
    qshift = [0,0,0]
    qgrid = kgrid_coarse
    dir_temp = '2-wfn'
    write_QE_input(f'{dir_temp}/bands.in', nbnds_temp, prefix, Natoms, Ntypes, pseudo_dir, 'bands', kpoints_file_wfn, dft_ekin_cutoff, qe_cell_file, qe_atoms_file)
    write_pw2bgw_input(f'{dir_temp}/pw2bgw.in', nbnds_temp, prefix, qshift, qgrid)
    
    # 3-wfnq
    nbnds_temp = nval + ncond_wfnq
    qshift = qshift_wfnq
    qgrid = kgrid_coarse
    dir_temp = '3-wfnq'
    write_QE_input(f'{dir_temp}/bands.in', nbnds_temp, prefix, Natoms, Ntypes, pseudo_dir, 'bands', kpoints_file_wfnq, dft_ekin_cutoff, qe_cell_file, qe_atoms_file)
    write_pw2bgw_input(f'{dir_temp}/pw2bgw.in', nbnds_temp, prefix, qshift, qgrid)
     
    # 4-wfn_co
    nbnds_temp = nval + ncond_wfnco
    qshift = [0,0,0]
    qgrid = kgrid_coarse
    dir_temp = '4-wfn_co'
    write_QE_input(f'{dir_temp}/bands.in', nbnds_temp, prefix, Natoms, Ntypes, pseudo_dir, 'bands', kpoints_file_wfn, dft_ekin_cutoff, qe_cell_file, qe_atoms_file)
    write_pw2bgw_input(f'{dir_temp}/pw2bgw.in', nbnds_temp, prefix, qshift, qgrid)
        
    # 5-wfn_fi
    nbnds_temp = nval + ncond_bse + 1
    qshift = [0,0,0]
    qgrid = kgrid_coarse
    dir_temp = '5-wfn_fi'
    write_QE_input(f'{dir_temp}/bands.in', nbnds_temp, prefix, Natoms, Ntypes, pseudo_dir, 'bands', kpoints_file_wfn_fi, dft_ekin_cutoff, qe_cell_file, qe_atoms_file)
    write_pw2bgw_input(f'{dir_temp}/pw2bgw.in', nbnds_temp, prefix, qshift, qgrid)        

    # 6-wfnq_fi
    nbnds_temp = nval,
    qshift = qshift_bse
    qgrid = kgrid_coarse
    dir_temp = '6-wfnq_fi'
    write_QE_input(f'{dir_temp}/bands.in', nbnds_temp, prefix, Natoms, Ntypes, pseudo_dir, 'bands', kpoints_file_wfnq_fi, dft_ekin_cutoff, qe_cell_file, qe_atoms_file)
    write_pw2bgw_input(f'{dir_temp}/pw2bgw.in', nbnds_temp, prefix, qshift, qgrid)  

    # writing BGW inputs
    print('Writing GW/BSE inputs')

    os.system('mkdir -p 7-epsilon 8-sigma 9-kernel 10.1-absorption_mom 10.2-absorption_vel')

    # epsilon
    print('Linking files for epsilon')
    os.system('ln -sf ../2-wfn/wfn.complex 7-epsilon/WFN')
    os.system('ln -sf ../3-wfnq/wfn.complex 7-epsilon/WFNq')

    # sigma
    print('Linking files for sigma')
    os.system('ln -sf ../2-wfn/wfn.complex 8-sigma/WFN_inner')
    os.system('ln -sf ../2-wfn/rho.complex 8-sigma/RHO')
    os.system('ln -sf ../4-wfn_co/vxc.dat 8-sigma/vxc.dat')
    os.system('ln -sf ../4-wfn_co/wfn.complex 8-sigma/WFN_outer')
    os.system('ln -sf ../7-epsilon/eps0mat.h5 8-sigma/eps0mat.h5')
    os.system('ln -sf ../7-epsilon/epsmat.h5 8-sigma/epsmat.h5')

    # kernel
    print('Linking files for kernel')
    os.system('ln -sf ../4-wfn_co/wfn.complex 9-kernel/WFN_co')
    os.system('ln -sf ../7-epsilon/eps0mat.h5 9-kernel/eps0mat.h5')
    os.system('ln -sf ../7-epsilon/epsmat.h5 9-kernel/epsmat.h5')

    # absorption momentum
    print('Linking files for absorption - momentum operator')
    os.system('ln -sf ../4-wfn_co/wfn.complex 10.1-absorption_mom/WFN_co')
    os.system('ln -sf ../1-scf/wfn.complex 10.1-absorption_mom/WFN_fi')
    os.system('ln -sf ../7-epsilon/eps0mat.h5 10.1-absorption_mom/eps0mat.h5')
    os.system('ln -sf ../7-epsilon/epsmat.h5 10.1-absorption_mom/epsmat.h5')
    os.system('ln -sf ../8-sigma/eqp1.dat 10.1-absorption_mom/eqp_co.dat')
    os.system('ln -sf ../9-kernel/bsemat.h5 10.1-absorption_mom/bsemat.h5')

    print('Linking files for absorption - velocity operator')
    os.system('ln -sf ../4-wfn_co/wfn.complex 10.2-absorption_vel/WFN_co')
    os.system('ln -sf ../5-wfn_fi/wfn.complex 10.2-absorption_vel/WFN_fi')
    os.system('ln -sf ../6-wfnq_fi/wfn.complex 10.2-absorption_vel/WFNq_fi')
    os.system('ln -sf ../7-epsilon/eps0mat.h5 10.2-absorption_vel/eps0mat.h5')
    os.system('ln -sf ../7-epsilon/epsmat.h5 10.2-absorption_vel/epsmat.h5')
    os.system('ln -sf ../8-sigma/eqp1.dat 10.2-absorption_vel/eqp_co.dat')
    os.system('ln -sf ../9-kernel/bsemat.h5 10.2-absorption_vel/bsemat.h5')
    
    write_epsilon('7-epsilon/epsilon.inp', kpoints_file_wfn, epsilon_cutoff, truncation_scheme, qshift_wfnq)
    write_sigma('8-sigma/sigma.inp', nmin_gw, nmax_gw, kpoints_file_wfn, truncation_scheme)
    write_kernel('9-kernel/kernel.inp', nval_kernel, ncond_kernel, truncation_scheme)
    write_absorption('10.1-absorption_mom/absorption.inp', nval_bse, ncond_bse, nval_kernel, ncond_kernel, truncation_scheme, use_momentum=True)
    write_absorption('10.2-absorption_vel/absorption.inp', nval_bse, ncond_bse, nval_kernel, ncond_kernel, truncation_scheme, use_momentum=False)

# print('Done!')

GWBSE_write_files()