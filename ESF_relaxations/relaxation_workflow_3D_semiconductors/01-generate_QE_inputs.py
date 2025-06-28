

from ase.io import read, write
from ase import Atoms
import subprocess

'''
Generate QE inputs for GWBSE calculations
1-scf
2-wfn_co 
3-wfnq_co
'''

pseudodir = '/home1/08948/rdgrande/pseudos/nc-sr-04_pw_standard/'
E_cutoff = 60 # Ry
nbnds = 160
cell_file = '../CELL'
atoms_file = 'ATOMS'
dirs_to_create = ['1-scf', '2-wfn_co', '3-wfnq_co']
calc_types = ['scf', 'bands', 'bands']
prefix = 'si'
kgrids = [(2, 2, 2, 0, 0, 0),    # 1-scf
          (1, 1, 1, 0, 0, 0),    # 2-wfn_co
          (1, 1, 1, 1e-3, 0, 0)]  # 3-wfnq_co

# Define pseudopotentials
pseudopotentials = {
   "Si": "Si.upf"
}

def write_ph_input(filename, prefix):

    text = f"""phonon_calc
$inputph
    verbosity = 'high'
    prefix = '{prefix}'
    outdir = './'
    fildyn = 'dyn'
    fildvscf = 'dvscf'
    electron_phonon='simple'
    trans=.true.
    nogg=.true.
    search_sym=.false.
/
0.0 0.0 0.0"""

    arq = open(filename, 'w')
    arq.write(text)
    arq.close()
    print(f"Finished writing {filename}")

def write_dynmat_input(filename):
    text = """$input
fildyn='dyn'
asr='crystal'
fileig='eigvecs'
filxsf='dynmat.axsf'
/"""
    arq = open(filename, 'w')
    arq.write(text)
    arq.close()
    print(f"Finished writing {filename}")

def read_cell_parameters(filename):
    """
    Reads lattice vectors from a Quantum ESPRESSO-style CELL_PARAMETERS block.

    Returns
    -------
    cell : list of lists
        A 3×3 list representing the lattice vectors.
    """
    cell = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.strip().startswith("CELL_PARAMETERS"):
                cell_lines = lines[i+1:i+4]
                for vec in cell_lines:
                    parts = vec.strip().split()
                    if len(parts) == 3:
                        cell.append([float(x) for x in parts])
                break
    return cell

def read_atomic_positions(filename):
    """
    Reads atomic positions from a Quantum ESPRESSO-style ATOMIC_POSITIONS block.

    Returns
    -------
    symbols : list of str
        Chemical symbols of atoms.
    positions : list of lists
        Atomic positions as floats in Å.
    """
    symbols = []
    positions = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.strip().startswith("ATOMIC_POSITIONS"):
                for l in lines[i+1:]:
                    parts = l.strip().split()
                    if len(parts) >= 4:
                        symbols.append(parts[0])
                        positions.append([float(x) for x in parts[1:4]])
                    elif len(parts) == 0:
                        break
                break
    return symbols, positions

def create_atoms(atoms_file, cell_file):
    """
    Creates an ASE Atoms object from a Quantum ESPRESSO input file.

    Returns
    -------
    atoms : ase.Atoms
    """
    cell = read_cell_parameters(cell_file)
    symbols, positions = read_atomic_positions(atoms_file)
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    return atoms

def dict_qe(PREFIX, calc_type, nbnds, pseudodir, E_cutoff, conv_thr=1e-10):
    
    dict_qe = {
            "control": {
                "calculation": calc_type,
                "prefix": PREFIX,
                "pseudo_dir": pseudodir,
                "verbosity": "high",
                "tstress": True,
                "tprnfor": True,
                "outdir": "./"
            },
            "system": {
                "ecutwfc": E_cutoff,
                "nbnd": nbnds,
                "occupations": "smearing",
                "smearing": "gaussian",
                "degauss": 0.001,
                "nosym": True
            },
            "electrons": {
                "conv_thr": conv_thr,
                'diagonalization': 'david',
                'diago_david_ndim': 4,
                'diago_full_acc': True
            }
        }
    return dict_qe  

def generate_kpoints(kgrid):
    Nkx, Nky, Nkz, qx, qy, qz = kgrid
    
    Nk_tot = Nkx * Nky * Nkz

    text = "K_POINTS crystal\n"
    text += f"{Nk_tot}\n"

    for ik in range(Nkx):
        for jk in range(Nky):
            for kk in range(Nkz):
                text += f"{(ik/Nkx+qx):.8f}  {(jk/Nky+qy):.8f}  {(kk/Nkz+qz):.8f}  1.0\n"
    return text

def write_qe_input_file(input_file, atoms, pseudopotentials, dict_conf, kpoints):
    """
    Writes a Quantum ESPRESSO input file.

    Parameters
    ----------
    input_file : str
        Path to the output file.
    atoms : ase.Atoms
        Atoms object containing atomic positions and cell parameters.
    pseudopotentials : dict
        Dictionary of pseudopotentials.
    dict_conf : dict
        Configuration dictionary for Quantum ESPRESSO.
    """
    write(input_file, atoms, format="espresso-in", 
          pseudopotentials=pseudopotentials, input_data=dict_conf)
    
    # Removing K_POINTS gamma written by ASE
    # sed command to delete the line starting with 'K_POINTS' and the next line
    cmd = f"sed '/^K_POINTS/,+1d' {input_file} > tmp && mv tmp {input_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Write KPOINTS
    with open(input_file, 'a') as f:
        f.write(kpoints)

    print(f"Finished writing {input_file}")
    
def write_pw2bgw_input_file(pw2bgw_input, prefix, kgrid, extra_lines=None):
    kx, ky, kz, qx, qy, qz = kgrid
    
    text = f'''&input_pw2bgw
prefix = '{prefix}'
real_or_complex = 2
wfng_flag = .true.
wfng_file = 'wfn.complex'
wfng_kgrid = .true.
wfng_nk1 = {kx}
wfng_nk2 = {ky}
wfng_nk3 = {kz}
wfng_dk1 = {qx}
wfng_dk2 = {qy}
wfng_dk3 = {qz}'''
    if extra_lines:
        for line in extra_lines:
            text += f'\n{line}'
    text += '\n/'
    
    arq = open(pw2bgw_input, 'w')
    arq.write(text)
    arq.close()
    print(f"Finished writing {pw2bgw_input}")




# running code

# Read cell and atomic positions from files ATOMS and CELL
atoms = create_atoms(atoms_file, cell_file)

for i, dir in enumerate(dirs_to_create):
    input_file = dir + '/qe.in'
    calc_type = calc_types[i]
    dict_conf = dict_qe(PREFIX=prefix, calc_type=calc_type, nbnds=nbnds, 
                        pseudodir=pseudodir, E_cutoff=E_cutoff)
    kpoints = generate_kpoints(kgrids[i])
    write_qe_input_file(input_file, atoms, 
                        pseudopotentials, 
                        dict_conf, kpoints)
    
# Write pw2bgw input file
write_pw2bgw_input_file('1-scf/pw2bgw.in', prefix, kgrids[0])
extra_lines =  ["rhog_flag = .true.",
                "rhog_file = 'RHO'",
                "vxcg_flag = .true.",
                "vxcg_file = 'VXC'",
                "vxc_flag = .false.",
                "vscg_flag = .true.",
                "vscg_file = 'VSC'",
                "vkbg_flag = .true.",
                "vkbg_file = 'VKB'"] 
write_pw2bgw_input_file('2-wfn_co/pw2bgw.in', prefix, kgrids[1], extra_lines)
write_pw2bgw_input_file('3-wfnq_co/pw2bgw.in', prefix, kgrids[2])

write_ph_input('1-scf/ph.in', prefix)
write_dynmat_input('1-scf/dynmat.in')
