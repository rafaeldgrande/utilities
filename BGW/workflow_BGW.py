
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



Ntypes, Natoms = 0, 0
arq = open(QE_ATOMS)
for line in arq:
    line_split = line.split()
    if len(line_split) == 3:  #C      12.0107   C.upf
        Ntypes += 1
    if len(line_split) == 4:  #C            0.5327850000       0.2500000000       0.9372530000
        Natoms += 1


# writing QE inputs

os.system('mkdir -p 1-scf 2-wfn 3-wfnq 4-wfn_co 5-wfn_fi 6-wfnq_fi')

for WFN_DIR in ['2-wfn', '3-wfnq', '4-wfn_co', '5-wfn_fi', '6-wfnq_fi']:
    os.system(f'mkdir -p {WFN_DIR}/{PREFIX}.save')
    os.system(f'ln -sf ../../1-scf/{PREFIX}.save/charge-density.hdf5 {WFN_DIR}/{PREFIX}.save/charge-density.hdf5')
    os.system(f'ln -sf ../../1-scf/{PREFIX}.save/data-file-schema.xml {WFN_DIR}/{PREFIX}.save/data-file-schema.xml')

write_QE_input('1-scf/scf.in', Nval, 'scf', KPOINTS_FILE_SCF)

write_QE_input('2-wfn/wfn.in', Nval + Ncond_WFN, 'bands', KPOINTS_FILE_WFN)
write_pw2bgw_input('2-wfn/pw2bgw.in', Nval + Ncond_WFN, [0,0,0], Kgrid_coarse)

write_QE_input('3-wfnq/wfn.in', Nval + Ncond_WFNq, 'bands', KPOINTS_FILE_WFNq)
write_pw2bgw_input('3-wfnq/pw2bgw.in', Nval + Ncond_WFNq, qshift_WFNq, Kgrid_coarse)

write_QE_input('4-wfn_co/wfn.in', Nval + Ncond_WFNco, 'bands', KPOINTS_FILE_WFN)
write_pw2bgw_input('4-wfn_co/pw2bgw.in', Nval + Ncond_WFNco, [0,0,0], Kgrid_coarse)

write_QE_input('5-wfn_fi/wfn.in', Nval + Ncond_WFNfi, 'bands', KPOINTS_FILE_WFN_fi)
write_pw2bgw_input('5-wfn_fi/pw2bgw.in', Nval + Ncond_WFNfi, [0,0,0], Kgrid_coarse)
write_ph_input('5-wfn_fi/ph.in')

write_QE_input('6-wfnq_fi/wfn.in', Nval + Ncond_WFNfi, 'bands', KPOINTS_FILE_WFNq_fi)
write_pw2bgw_input('6-wfnq_fi/pw2bgw.in', Nval + Ncond_WFNfi, qshift_WFNq, Kgrid_coarse)

# writing BGW inputs

os.system('mkdir -p 7-epsilon 8-sigma 9-kernel 10.1-absorption_mom 10.2-absorption_vel')

# epsilon
os.system('ln -sf ../2-wfn/wfn.complex 7-epsilon/WFN')
os.system('ln -sf ../3-wfnq/wfn.complex 7-epsilon/WFN')

# sigma
os.system('ln -sf ../2-wfn/wfn.complex 8-sigma/WFN_inner')
os.system('ln -sf ../2-wfn/rho.complex 8-sigma/RHO')
os.system('ln -sf ../2-wfn/vxc.dat 8-sigma/vxc.dat')
os.system('ln -sf ../7-epsilon/eps0mat.h5 8-sigma/eps0mat.h5')
os.system('ln -sf ../7-epsilon/epsmat.h5 8-sigma/epsmat.h5')

# kernel
os.system('ln -sf ../4-wfn_co/wfn.complex 9-kernel/WFN_co')
os.system('ln -sf ../7-epsilon/eps0mat.h5 9-kernel/eps0mat.h5')
os.system('ln -sf ../7-epsilon/epsmat.h5 9-kernel/epsmat.h5')

# absorption momentum
os.system('ln -sf ../4-wfn_co/wfn.complex 10.1-absorption_mom/WFN_co')
os.system('ln -sf ../5-wfn_fi/wfn.complex 10.1-absorption_mom/WFN_fi')
os.system('ln -sf ../7-epsilon/eps0mat.h5 10.1-absorption_mom/eps0mat.h5')
os.system('ln -sf ../7-epsilon/epsmat.h5 10.1-absorption_mom/epsmat.h5')
os.system('ln -sf ../8-sigma/eqp1.dat 10.1-absorption_mom/eqp_co.dat')
os.system('ln -sf ../9-kernel/bsemat.h5 10.1-absorption_mom/bsemat.h5')

# absorption velocity
os.system('ln -sf ../4-wfn_co/wfn.complex 10.2-absorption_vel/WFN_co')
os.system('ln -sf ../5-wfn_fi/wfn.complex 10.2-absorption_vel/WFN_fi')
os.system('ln -sf ../7-epsilon/eps0mat.h5 10.2-absorption_vel/eps0mat.h5')
os.system('ln -sf ../7-epsilon/epsmat.h5 10.2-absorption_vel/epsmat.h5')
os.system('ln -sf ../8-sigma/eqp1.dat 10.2-absorption_vel/eqp_co.dat')
os.system('ln -sf ../9-kernel/bsemat.h5 10.2-absorption_vel/bsemat.h5')

write_epsilon('7-epsilon/epsilon.inp')
write_sigma('8-sigma/sigma.inp', NminGW, NmaxGW)
write_kernel('9-kernel/kernel.inp', nvalKernel, ncondKernel)
write_absorption('10.1-absorption_mom/absorption.inp', nvalBSE, ncondBSE, nvalKernel, ncondKernel, use_momentum=True)
write_absorption('10.2-absorption_vel/absorption.inp', nvalBSE, ncondBSE, nvalKernel, ncondKernel, use_momentum=False)

print('Done!')