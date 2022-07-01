# QE functions
def get_scf_energy(file_out):

    arq = open(file_out)

    for line in arq:
        linha = line.split()

        if len(linha) > 0:
            if linha[0] == '!':
                E = np.float(linha[4])

    return E

#####################################
###### Writing functions ############
#####################################

def write_QE_CONTROL(PREFIX, pseudo_dir, calc_type):

    text = '&CONTROL \n'
    text += '    prefix = \''+PREFIX+'\' \n'

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

    return text


def write_QE_SYSTEM(ecutDFT, nbnds, nat, ntyp):

    text = ('/\n'+
    '&SYSTEM \n'+
    '    ibrav = 0 \n'
    f'    nat = {nat} \n'+
    f'    ntyp = {ntyp} \n'+
    '    nbnd = '+str(nbnds)+' \n'+
    f'    ecutwfc = {ecutDFT}'+'\n'+
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

def write_CELL_PARAMETERS(A, B, C):
    text = 'CELL_PARAMETERS (angstrom) \n'
    text += f'{A}   0   0 \n'
    text += f'0   {B}   0 \n'
    text += f'0   0   {C} \n'

    return text

# ATOMIC_SPECIES = ['C', 'H', ...]
# ATOMIC_MASSES = [12.011, 1.008, ...]
# PSEUDOS_FILES = ['C.upf', 'H.upf', ...]

def write_QE_ATOMS(ATOMIC_SPECIES, ATOMIC_MASSES, PSEUDOS_FILES, ATOMIC_POSITIONS):
    text = 'ATOMIC_SPECIES \n'
    for i_atomic in range(len(ATOMIC_SPECIES)):
        text += f'   {ATOMIC_SPECIES[i_atomic]}   {ATOMIC_MASSES[i_atomic]}   {PSEUDOS_FILES[i_atomic]}\n'
    text += 'ATOMIC_POSITIONS angstrom \n'
    for i_atomic in range(len(ATOMIC_POSITIONS)):
        text += ATOMIC_POSITIONS[i_atomic]
    return text

def write_QE_KPOINTS(qshift, Nkgrid):

    """Nkgrid = (Nkx, Nky, Nkz)  -> crystal units (dumb way of doing it yet)
    qshift = (qx, qy, qz)"""

    Nkx, Nky, Nkz = Nkgrid
    qx, qy, qz = qshift

    Nktotal = Nkx*Nky*Nkz
    
    text = 'K_POINTS crystal\n'
    text += f'{Nktotal} \n'
    for ikx in range(Nkx):
        for iky in range(Nky):
            for ikz in range(Nkz):
                kx = round(ikx/Nkx + qx, 9)
                ky = round(iky/Nky + qy, 9)
                kz = round(ikz/Nkz + qz, 9)
                text += f'{kx} {ky} {kz} 1 \n'

    return text

def write_pw_input(in_pw_name, PREFIX, A, B, C, ecutDFT, nbnds, pseudo_dir, calc_type, ATOMIC_SPECIES, ATOMIC_MASSES, PSEUDOS_FILES, ATOMIC_POSITIONS, qshift, Nkgrid):

    nat, ntyp = len(ATOMIC_POSITIONS), len(ATOMIC_SPECIES)
    text = write_QE_CONTROL(PREFIX, pseudo_dir, calc_type)
    text += write_QE_SYSTEM(ecutDFT, nbnds, nat, ntyp)
    text += write_QE_ELECTRONS()
    text += write_CELL_PARAMETERS(A, B, C)
    text += write_QE_ATOMS(ATOMIC_SPECIES, ATOMIC_MASSES, PSEUDOS_FILES, ATOMIC_POSITIONS)
    text += write_QE_KPOINTS(qshift, Nkgrid)

    arq = open(in_pw_name, 'w')
    arq.write(text)
    arq.close()

def write_ph_input(PREFIX, ph_input, ph_tr_exp):

    text = 'phonon_calc \n'
    text += '$inputph \n'
    text += '    verbosity = \'high\'\n'
    text += '    prefix = \''+PREFIX+'\'\n'
    text += '    outdir = \'./\'\n'
    text += '    fildyn = \'dyn\'\n'
    text += '    fildvscf = \'dvscf\'\n'
    text += '    electron_phonon=\'simple\'\n'
    text += '    trans=.true.\n'
#    text += '    asr=.true.\n'
    text += '    nogg=.true.\n'
#    text += '    el_ph_nsigma = 15 \n'
#    text += '    el_ph_sigma = 0.002 \n' 
    text += f'    tr2_ph=1.0d-{ph_tr_exp}'+'\n'
    text += '/\n'
    text += '0.0 0.0 0.0'

    arq = open(ph_input, 'w')
    arq.write(text)
    arq.close()

def write_pw2bgw_input(pw2bgw_inp_file, PREFIX, nbnds, qshift, Nkgrid):

    Nkx, Nky, Nkz = Nkgrid
    qx, qy, qz = qshift

    text = '&input_pw2bgw \n'
    text += 'prefix = \'' + f'{PREFIX}'+ '\' \n'
    text += 'real_or_complex = 2 \n'
    text += 'wfng_flag = .true. \n'
    text += 'wfng_file = \'wfn.complex\' \n'
    text += 'wfng_kgrid = .true. \n'
    text += f'wfng_nk1 = {Nkx} '+'\n'
    text += f'wfng_nk2 = {Nky} '+'\n'
    text += f'wfng_nk3 = {Nkz} '+'\n'
    text += f'wfng_dk1 = {qx} '+'\n'
    text += f'wfng_dk2 = {qy} '+'\n'
    text += f'wfng_dk3 = {qz} '+'\n'
    text += 'rhog_flag = .true. \n'
    text += 'rhog_file = \'rho.complex\' \n'
    text += 'vxcg_flag = .false. \n'
    text += 'vxcg_file = \'vxc.complex\' \n'
    text += 'vxc_flag = .true. \n'
    text += 'vxc_file = \'vxc.dat\' \n'
    text += 'vxc_diag_nmin = 1 \n'
    text += 'vxc_diag_nmax = '+str(nbnds)+' \n'
    text += 'vxc_offdiag_nmin = 1 \n'
    text += 'vxc_offdiag_nmax = '+str(nbnds)+' \n'
    text += '/\n'

    arq = open(pw2bgw_inp_file, 'w')
    arq.write(text)
    arq.close()

def read_atomic_pos(arq_atomic_pos):
    ATOMIC_POS = []
    arq = open(arq_atomic_pos)
    for line in arq:
        ATOMIC_POS.append(line)
    return ATOMIC_POS
