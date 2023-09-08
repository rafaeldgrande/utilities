

from workflow_BGW_config import *

# QE functions

def write_QE_CONTROL(calc_type, prefix, pseudo_dir, extra_lines):

    text = '&CONTROL \n'
    text += f"    prefix = '{prefix}' \n"
    if calc_type == 'scf':
        text += '    calculation = \'scf\' \n'        
    elif calc_type == 'bands':
        text += '    calculation = \'bands\' \n'  
    text += '    wf_collect = .true. \n'  # in the manual -> OBSOLETE - NO LONGER IMPLEMENTED      
    text += '    outdir = \'./\' \n'
    text += '    wfcdir = \'./\' \n'
    text += '    tstress = .true. \n'
    text += '    tprnfor =  .true. \n'
    text += '    pseudo_dir = '+pseudo_dir+'\n'
    for line in extra_lines:
        text += line + '\n'

    return text

def write_QE_SYSTEM(nbnds, dft_Ekin_cutoff, Natoms, Ntypes):

    text = ('/\n'+
    '&SYSTEM \n'+
    '    ibrav = 0 \n'+  # crystal axis provided in input: see card CELL_PARAMETERS
    f'    nat = {Natoms} \n'+
    f'    ntyp = {Ntypes} \n'+
    '    nbnd = '+str(nbnds)+' \n'+
    f'    ecutwfc = {dft_Ekin_cutoff} \n'+
    '    occupations = \'smearing\' \n'+
    '    smearing = \'gaussian\' \n'+
    '    degauss = 0.001 \n'+
    '    nosym = .true. \n')
    
    return text

def write_QE_ELECTRONS():
    text = ('/\n'+
    '&ELECTRONS \n'+
    '    electron_maxstep = 100 \n'+
    '    mixing_mode = \'plain\' \n'+
    '    mixing_beta = 0.7 \n'+
    '    mixing_ndim = 8 \n'+
    '    diagonalization = \'david\' \n'+
    '    diago_david_ndim = 4 \n'+
    '    diago_full_acc = .true. \n'+
    '/\n')

    return text

def write_QE_ATOMS(qe_atoms_file):
    
    # it is a file like this
    #     ATOMIC_SPECIES
    #     C  12.011   C.upf
    #     O  15.999   O.upf
    # ATOMIC_POSITIONS angstrom
    # C             5.0    5.0    5.56
    # O             5.0    5.0    4.43
    
    arq = open(qe_atoms_file)
    text = arq.read()

    return text

def write_QE_CELL(qe_cell_file):
    
    # it is a file like this
    #     CELL_PARAMETERS angstrom
    # 10 0 0
    # 0 10 0
    # 0 0 10
    arq = open(qe_cell_file)
    text = arq.read()

    return text

def write_QE_KPOINTS(KPOINTS_FILE):  
    arq = open(KPOINTS_FILE)
    text = arq.read()

    return text


def write_QE_input(in_pw_name, nbnds, prefix, Natoms, Ntypes, pseudo_dir, calc_type, KPOINTS_FILE, dft_Ekin_cutoff, QE_CELL, QE_ATOMS):

    text = write_QE_CONTROL(calc_type, prefix, pseudo_dir, [])
    text += write_QE_SYSTEM(nbnds, dft_Ekin_cutoff, Natoms, Ntypes)
    text += write_QE_ELECTRONS()
    text += write_QE_CELL(QE_CELL)
    text += write_QE_ATOMS(QE_ATOMS)
    text += write_QE_KPOINTS(KPOINTS_FILE)

    arq = open(in_pw_name, 'w')
    arq.write(text)
    arq.close()
    
    print(f'Finished writing QE input :{in_pw_name}')
    

def write_pw2bgw_input(pw2bgw_inp, nbnds, prefix, qshift, qgrid):
    
    qx, qy, qz = qgrid
    qshift_x, qshift_y, qshift_z = qshift

    text = '&input_pw2bgw \n'
    text += f"prefix = '{prefix}' \n"
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
    text += 'vxc_offdiag_nmin = 1 \n'
    text += 'vxc_offdiag_nmax = '+str(nbnds)+' \n'
    text += '/\n'

    arq = open(pw2bgw_inp, 'w')
    arq.write(text)
    arq.close()
    
    print(f'Finished writing pw2bgw input :{pw2bgw_inp}')
    
def write_ph_input(ph_input, prefix):
    text = f'''phonon_calc
    $inputph
    verbosity = 'high'
    prefix = '{prefix}'
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

    arq = open(ph_input, 'w')
    arq.write(text)
    arq.close()
    
    print(f'Finished writing ph input :{ph_input}')

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

def generate_regular_kgrid(KGRID_FILE, kgrid, qshift):
    
    arq = open(KGRID_FILE, 'w')
    
    Nkx, Nky, Nkz = kgrid
    qx, qy, qz = qshift
    
    Nk_tot = Nkx * Nky * Nkz

    arq.write("K_POINTS crystal\n")
    arq.write(f"{Nk_tot}\n")

    for ik in range(Nkx):
        for jk in range(Nky):
            for kk in range(Nkz):
                arq.write(f"{(ik/Nkx+qx):.8f}  {(jk/Nky+qy):.8f}  {(kk/Nkz+qz):.8f}  1.0\n")            
                
    arq.close()
    print(f'Finished generating kgrid {Nkx}x{Nky}x{Nkz} with qshift=({qx},{qy},{qz}) and wrote it on file {KGRID_FILE}')

def truncation_scheme_input(truncation_scheme):
        # truncation scheme
    if truncation_scheme == 'box':
        return 'cell_box_truncation \n\n'
    elif truncation_scheme == 'wire':
        return 'cell_wire_truncation \n\n'
    elif truncation_scheme == 'slab':
        return 'cell_slab_truncation \n\n'
    else:
        return '\n'


def write_epsilon(epsilon_inp, KPOINTS_FILE, epsilon_cutoff, truncation_scheme, qshift):
    
    Kpoints = get_kpoints(KPOINTS_FILE)

    text = '\n\nverbosity 3 \n\n'
    text += f'epsilon_cutoff {epsilon_cutoff} \n\n'
    text += 'degeneracy_check_override \n\n'
    
    text += truncation_scheme_input(truncation_scheme)
    
    text += 'begin qpoints\n'
    for kpoint in Kpoints:
        if kpoint == [0.0, 0.0, 0.0]:
            qx, qy, qz = qshift
            text += f' {qx}  {qy}  {qz} 1.0 1 \n'
        else:
            kx, ky, kz = kpoint
            text += f' {kx}  {ky}  {kz}   1.0 0 \n'
    text += 'end \n'

    arq = open(epsilon_inp, 'w')
    arq.write(text)
    arq.close()
    
    print(f'Finished writing epsilon input file {epsilon_inp}')
    
    
def write_sigma(sigma_inp, NminGW, NmaxGW, KPOINTS_FILE, truncation_scheme):
    
    Kpoints = get_kpoints(KPOINTS_FILE)

    text = '\n\nverbosity 3 \n\n'
    text += '\n\degeneracy_check_override \n\n'
    text += 'screening_semiconductor \n\n'
    text += 'exact_static_ch 1 \n\n'
    text += truncation_scheme_input(truncation_scheme)
    text += f'band_index_min {NminGW}\n'
    text += f'band_index_max {NmaxGW}'+'\n\n'    
    text += 'begin kpoints\n'
    for kpoint in Kpoints:
        kx, ky, kz = kpoint
        text += f' {kx}  {ky}  {kz}   1.0 \n'
    text += 'end \n'

    arq = open(sigma_inp, 'w')
    arq.write(text)
    arq.close()
    
    print(f'Finished writing sigma input file {sigma_inp}')

    
def write_kernel(kernel_inp, nvalKernel, ncondKernel, truncation_scheme):

    text = '\n\nverbosity 3 \n\n'
    text +=  f'number_val_bands {nvalKernel} '+'\n'
    text += f'number_cond_bands {ncondKernel} '+'\n\n'
    text += truncation_scheme_input(truncation_scheme)
    text += 'use_symmetries_coarse_grid \n\n'
    text += 'screening_semiconductor'

    arq = open(kernel_inp, 'w')
    arq.write(text)
    arq.close()
    
    print(f'Finished writing kernel input file {kernel_inp}')



def write_absorption(abs_inp, nvalBSE, ncondBSE, nvalKernel, ncondKernel, truncation_scheme, use_momentum):

    text = '\n\nverbosity 3 \n\n'
    text += 'diagonalization \n\n'

    text += f'number_val_bands_coarse {nvalKernel} '+'\n'
    text += f'number_val_bands_fine {nvalBSE}'+' \n'
    text += f'number_cond_bands_coarse {ncondKernel}'+' \n'
    text += f'number_cond_bands_fine {ncondBSE}'+' \n\n'

    text += 'use_symmetries_coarse_grid \n'
    text += 'no_symmetries_fine_grid \n'
    text += 'no_symmetries_shifted_grid \n\n'
    
    text += truncation_scheme_input(truncation_scheme)    

    text += 'screening_semiconductor \n'
    text += 'eqp_co_corrections\n\n'

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
    
    print(f'Finished writing absorpiton file {abs_inp}')

