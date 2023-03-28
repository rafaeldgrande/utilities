
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
        
# def read_test(input_file):
#     arq = open(input_file)
#     for line in arq:
#         print(line)
# read_test('workflow_input')


# reading input file
# arq = open('workflow_input')

# for line in arq:
#     line_split = line.split()
#     if len(line_split) > 1:
#         if line_split[0] == 'PREFIX':
#             PREFIX = line_split[1]
#         elif line_split[0] == 'pseudo_dir':
#             pseudo_dir = line_split[1]
#         elif line_split[0] == 'QE_ATOMS':
        
import configparser

config = configparser.ConfigParser()
config.read("workflow_input")

PREFIX               =  config.get("VARS", "PREFIX" )  
pseudo_dir           =  config.get("VARS", "pseudo_dir" )  
QE_ATOMS             =  config.get("VARS", "QE_ATOMS" )  
QE_CELL              =  config.get("VARS", "QE_CELL" )  
Kgrid_coarse         =  config.get("VARS", "Kgrid_coarse" )  
qshift_WFNq          =  config.get("VARS", "qshift_WFNq" )  
Kgrid_fine           =  config.get("VARS", "Kgrid_fine" )  
KPOINTS_FILE_SCF     =  config.get("VARS", "KPOINTS_FILE_SCF" )  
KPOINTS_FILE_WFN     =  config.get("VARS", "KPOINTS_FILE_WFN" ) 
KPOINTS_FILE_WFNq    =  config.get("VARS", "KPOINTS_FILE_WFNq" )  
KPOINTS_FILE_WFN_fi  =  config.get("VARS", "KPOINTS_FILE_WFN_fi" )  
KPOINTS_FILE_WFNq_fi =  config.get("VARS", "KPOINTS_FILE_WFNq_fi" )  
Nval                 =  config.get("VARS", "Nval" )  
Ncond_WFN            =  config.get("VARS", "Ncond_WFN" ) 
Ncond_WFNq           =  config.get("VARS", "Ncond_WFNq" )   
Ncond_WFNco          =  config.get("VARS", "Ncond_WFNco" )   
Ncond_WFNfi          =  config.get("VARS", "Ncond_WFNfi" )  
NminGW               =  config.get("VARS", "NminGW" )  
NmaxGW               =  config.get("VARS", "NmaxGW" )  
nvalKernel           =  config.get("VARS", "nvalKernel" )   
ncondKernel          =  config.get("VARS", "ncondKernel" )  
nvalBSE              =  config.get("VARS", "nvalBSE" )  
ncondBSE             =  config.get("VARS", "ncondBSE" )  

print('##########')
print(PREFIX)
print(QE_ATOMS)

# getting number of atoms and number of atomic types
Ntypes, Natoms = 0, 0
arq = open(QE_ATOMS)
for line in arq:
    line_split = line.split()
    if len(line_split) == 3:  #C      12.0107   C.upf
        Ntypes += 1
    if len(line_split) == 4:  #C            0.5327850000       0.2500000000       0.9372530000
        Natoms += 1
        
# 
temp = Kgrid_coarse.split()
kx, ky, kz = int(temp[0]), int(temp[1]), int(temp[2])
Kgrid_coarse = [kx, ky, kz]

temp = qshift_WFNq.split()
kx, ky, kz = float(temp[0]), float(temp[1]), float(temp[2])
qshift_WFNq = [kx, ky, kz]

temp = Kgrid_fine.split()
kx, ky, kz = int(temp[0]), int(temp[1]), int(temp[2])
Kgrid_fine = [kx, ky, kz]

Nval = int(Nval)
Ncond_WFN = int(Ncond_WFN)
Ncond_WFNq = int(Ncond_WFNq)
Ncond_WFNco = int(Ncond_WFNco)
Ncond_WFNfi = int(Ncond_WFNfi)
NminGW = int(NminGW)
NmaxGW = int(NmaxGW)
nvalKernel = int(nvalKernel)
ncondKernel = int(ncondKernel)
nvalBSE = int(nvalBSE)
ncondBSE = int(ncondBSE)
