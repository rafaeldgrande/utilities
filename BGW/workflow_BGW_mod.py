

from workflow_BGW_config import *

# QE functions

def write_QE_CONTROL(calc_type, extra_flags):

    text = '&CONTROL \n'
    text += f"    prefix = '{PREFIX}' \n"
    if calc_type == 'scf':
        text += '    calculation = \'scf\' \n'
        text += '    wf_collect = .true. \n'
    elif calc_type == 'bands':
        text += '    calculation = \'bands\' \n'
        text += '    wf_collect = .true. \n'        
    text += '    outdir = \'./\' \n'
    text += '    wfcdir = \'./\' \n'
    text += '    tprnfor =  .true. \n'
    text += '    pseudo_dir = '+pseudo_dir+'\n'
    for extra_line in extra_flags:
        text += extra_line + '\n'

    return text

def write_QE_SYSTEM(nbnds):

    text = ('/\n'+
    '&SYSTEM \n'+
    '    ibrav = 0 \n'+
    f'    nat = {Natoms} \n'+
    f'    ntyp = {Ntypes} \n'+
    '    nbnd = '+str(nbnds)+' \n'+
    '    ecutwfc = 100.0 \n'+
    '    occupations = \'smearing\' \n'+
    '    smearing = \'gaussian\' \n'+
    '    degauss = 0.001 \n')
    
    return text

def write_QE_ELECTRONS():
    text = ('/\n'+
    '&ELECTRONS \n'+
    '    electron_maxstep = 100 \n'+
    '    conv_thr = 1.0d-10 \n'+
    '    mixing_mode = \'plain\' \n'+
    '    mixing_beta = 0.7 \n'+
    '    mixing_ndim = 8 \n'+
    '    diagonalization = \'david\' \n'+
    '    diago_david_ndim = 4 \n'+
    '    diago_full_acc = .true. \n'+
    '/\n')

    return text

def write_QE_ATOMS(QE_ATOMS):
    arq = open(QE_ATOMS)
    text = arq.read()

    return text

def write_QE_CELL(QE_CELL):
    arq = open(QE_CELL)
    text = arq.read()

    return text

def write_QE_KPOINTS(KPOINTS_FILE):  
    arq = open(KPOINTS_FILE)
    text = arq.read()

    return text


def write_QE_input(in_pw_name, nbnds, calc_type, KPOINTS_FILE):

    text = write_QE_CONTROL(calc_type, [])
    text += write_QE_SYSTEM(nbnds)
    text += write_QE_ELECTRONS()
    text += write_QE_CELL(QE_CELL)
    text += write_QE_ATOMS(QE_ATOMS)
    text += write_QE_KPOINTS(KPOINTS_FILE)

    arq = open(in_pw_name, 'w')
    arq.write(text)
    arq.close()
    

def write_pw2bgw_input(pw2bgw_inp, nbnds, qshift, qgrid):
    
    qx, qy, qz = qgrid
    qshift_x, qshift_y, qshift_z = qshift

    text = '&input_pw2bgw \n'
    text += f'prefix = {PREFIX} \n'
    text += 'real_or_complex = 2 \n'
    text += 'wfng_flag = .true. \n'
    text += 'wfng_file = \'wfn.complex\' \n'
    text += 'wfng_kgrid = .true. \n'
    text += f'wfng_nk1 = {qx} \n'
    text += f'wfng_nk2 = {qy} \n'
    text += f'wfng_nk3 = {qz} \n'
    text += f'wfng_dk1 = {qshift_x} '+'\n'
    text += f'wfng_dk2 = {qshift_y} '+'\n'
    text += f'wfng_dk3 = {qshift_z} '+'\n'
    text += 'rhog_flag = .true. \n'
    text += 'rhog_file = \'rho.complex\' \n'
    text += 'vxcg_flag = .true. \n'
    text += 'vxcg_file = \'vxc.complex\' \n'
    text += 'vxc_flag = .true. \n'
    text += 'vxc_file = \'vxc.dat\' \n'
    text += 'vxc_diag_nmin = 1 \n'
    text += f'vxc_diag_nmax = {nbnds} \n'
#    text += 'vxc_offdiag_nmin = 1 \n'
#    text += 'vxc_offdiag_nmax = '+str(nbnds)+' \n'
    text += '/\n'

    arq = open(pw2bgw_inp, 'w')
    arq.write(text)
    arq.close()
    
def write_ph_input(ph_input):
    text = f'''$inputph
    verbosity = 'high'
    prefix = '{PREFIX}'
    outdir = './'
    fildyn = 'dyn'
    fildvscf = 'dvscf'
    electron_phonon='simple'
    trans=.true.
    nogg=.true.
    tr2_ph=1.0d-16
    max_seconds=1d8
    search_sym=.false.
/
0.0 0.0 0.0'''

def get_kpoints(KPOINTS_FILE):
    
    Kpoints_in_file = []
    arq = open(KPOINTS_FILE)
    text = arq.read().split('\n')
    
    for line in text:
        line_split = line.split()
        if len(line_split) == 4:
            kx, ky, kz = float(line_split[0]), float(line_split[1]), float(line_split[2])
            Kpoints_in_file.append([kx, ky, kz])
            
    return Kpoints_in_file
            

def write_epsilon(epsilon_inp):
    
    Kpoints = get_kpoints(KPOINTS_FILE_WFN)

    text = '\n\nepsilon_cutoff 16.0 \n\n'
    text += 'degeneracy_check_override \n\n'
    text += 'begin qpoints\n'
    for kpoint in Kpoints:
        if kpoint == [0.0, 0.0, 0.0]:
            text += ' 0.0 0.0 0.0 1.0 1 \n'
        else:
            kx, ky, kz = kpoint
            text += f' {kx}  {ky}  {kz}   1.0 0 \n'
    text += 'end \n'

    arq = open(epsilon_inp, 'w')
    arq.write(text)
    arq.close()
    
    
def write_sigma(sigma_inp, NminGW, NmaxGW):
    
    Kpoints = get_kpoints(KPOINTS_FILE_WFN)

    text = '\n\degeneracy_check_override 16.0 \n\n'
    text += 'screening_semiconductor \n\n'
    text += 'exact_static_ch 1 \n\n'
    text += f'band_index_min {NminGW}\n'
    text += f'band_index_max {NmaxGW}'+'\n\n'    
    text += 'begin kpoints\n'
    for kpoint in Kpoints:
        kx, ky, kz = kpoint
        text += f' {kx}  {ky}  {kz}   1.0 0 \n'
    text += 'end \n'

    arq = open(sigma_inp, 'w')
    arq.write(text)
    arq.close()
    
def write_kernel(kernel_inp, nvalKernel, ncondKernel):

    text =  f'number_val_bands {nvalKernel} '+'\n'
    text += f'number_cond_bands {ncondKernel} '+'\n\n'
    text += 'use_symmetries_coarse_grid \n\n'
    text += 'screening_semiconductor'

    arq = open(kernel_inp, 'w')
    arq.write(text)
    arq.close()


def write_absorption(abs_inp, nvalBSE, ncondBSE, nvalKernel, ncondKernel, use_momentum):

    text = '\n\nverbosity 3 \n\n'
    text += 'diagonalization \n\n'

#    text += 'degeneracy_check_override \n'

    text += f'number_val_bands_coarse {nvalKernel} '+'\n'
    text += f'number_val_bands_fine {nvalBSE}'+' \n'
    text += f'number_cond_bands_coarse {ncondKernel}'+' \n'
    text += f'number_cond_bands_fine {ncondBSE}'+' \n\n'

    text += 'use_symmetries_coarse_grid \n'
    text += 'no_symmetries_fine_grid \n'
    text += 'no_symmetries_shifted_grid \n\n'

    text += 'screening_semiconductor \n'
    text += 'eqp_co_corrections\n\n'

#    text += 'use_velocity \n'
    if use_momentum == True:
        text += 'use_momentum \n'
        text += 'polarization 0.0 0.0 1.0 \n\n'
    else:
        text += 'use_velocity \n'

    text += 'gaussian_broadening \n'
    text += 'energy_resolution 0.03 \n'
    text += 'delta_frequency 0.001 \n\n'
    
    text += 'write_eigenvectors -1'

    arq = open(abs_inp, 'w')
    arq.write(text)
    arq.close()

