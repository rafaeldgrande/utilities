

# QE parameters
PREFIX = 'MAPI'
pseudo_dir = './'

QE_ATOMS = 'ATOMS_FILE'
QE_CELL = 'CELL_QE'

Kgrid_coarse = [6,6,6]
qshift_WFNq = [0,0, 1e-3]
Kgrid_fine = [6,6,6]
KPOINTS_FILE_SCF     = 'kpoints_scf'
KPOINTS_FILE_WFN     = 'WFN.out'
KPOINTS_FILE_WFNq    = 'WFNq.out'
KPOINTS_FILE_WFN_fi  = 'WFN_fi.out'
KPOINTS_FILE_WFNq_fi = 'WFNq_fi.out'

# Bands numbers
Nval = 25   # highest valence band
Ncond_WFN = 500 # number of cond bands for WFN file
Ncond_WFNq = 50 # number of cond bands for WFNq file
Ncond_WFNco = 10 # WFN_co
Ncond_WFNfi = 10

# BGW params 
NminGW, NmaxGW = 90, 110
nvalKernel, ncondKernel = 10, 10
nvalBSE, ncondBSE = 5, 5 

Ntypes, Natoms = 0, 0
arq = open(QE_ATOMS)
for line in arq:
    line_split = line.split()
    if len(line_split) == 3:  #C      12.0107   C.upf
        Ntypes += 1
    if len(line_split) == 4:  #C            0.5327850000       0.2500000000       0.9372530000
        Natoms += 1