
import numpy as np
import subprocess

first_step = 0
last_step = 3
Natoms = 54

ry2ev = 13.605703976

STEPS = []
MAX_ES_FORCE, MEAN_ES_FORCE = [], []
MAX_DFT_FORCE, MEAN_DFT_FORCE = [], []
MAX_TOT_FORCE, MEAN_TOT_FORCE = [], []
EXC_ENERGY, SCF_ENERGY, TOT_ENERGY = [], [], []
MAX_DISPLACEMENT, MEAN_DISPLACEMENT = [], []

def read_dft_forces_qe(file, Natoms):
    
    # dft forces array - Ry/bohr
    dft_forces = np.zeros((3*Natoms))
    
    # getting data from QE output file
    grep_command = ["grep", "force", file]
    
    # result = subprocess.run(grep_command, stdout=subprocess.PIPE)
    # print(result.stdout)
    try:
        grep_output = subprocess.check_output(grep_command, stderr=subprocess.PIPE, text=True)
        print('Grep worked sucessfully!')
    except subprocess.CalledProcessError as e:
        print("Error executing grep:", e)
        print("Did not find the DFT forces!")
        print("DFT forces are set to 0.\n")
        return dft_forces
        
    # filtering - gtting the first Natoms lines
    temp_text = grep_output.split('\n')[:Natoms]
    
    # parsing data
    for iatom in range(Natoms):
        # '     atom    1 type  1   force =     0.00000000    0.00000000    0.02862271'
        line_split = temp_text[iatom].split()
        fx, fy, fz = float(line_split[-3]), float(line_split[-2]), float(line_split[-1])
        
        dft_forces[iatom * 3] = fx
        dft_forces[iatom * 3 + 1] = fy
        dft_forces[iatom * 3 + 2] = fz

    print("")
        
    return dft_forces

for istep in range(first_step, last_step+1):
    STEPS.append(float(istep))

    forces_file = f"step_{istep}/9-excited-state-forces/forces_cart.out_1_1"
    ESforces = np.loadtxt(forces_file, usecols=3, dtype=complex) # shape (3N)
    
    Natoms = len(ESforces) // 3
    mod_forces = []
    for iatom in range(Natoms):
        mod_forces.append(np.linalg.norm(ESforces[3*iatom:3*(iatom+1)]))
    MAX_ES_FORCE.append(max(mod_forces))
    MEAN_ES_FORCE.append(np.mean(mod_forces))
    
    dft_forces = read_dft_forces_qe(f"step_{istep}/1-scf/qe.out", Natoms)
    dft_mod_forces = []
    for iatom in range(Natoms):
        dft_mod_forces.append(np.linalg.norm(dft_forces[3*iatom:3*(iatom+1)]))
    MAX_DFT_FORCE.append(max(dft_mod_forces))
    MEAN_DFT_FORCE.append(np.mean(dft_mod_forces))
    
    tot_force = dft_forces + ESforces
    tot_mod_forces = []
    for iatom in range(Natoms):
        tot_mod_forces.append(np.linalg.norm(tot_force[3*iatom:3*(iatom+1)]))
    MAX_TOT_FORCE.append(max(tot_mod_forces))
    MEAN_TOT_FORCE.append(np.mean(tot_mod_forces))
    
    eigvals_file = f"step_{istep}/8-absorption/eigenvalues_b1.dat"
    eigvals = np.loadtxt(eigvals_file)
    EXC_ENERGY.append(eigvals[0, 0])
    
    qe_out_file = f"step_{istep}/1-scf/qe.out"
    arq = open(qe_out_file)
    for line in arq:
        if "total energy" in line:
            # print(line.split())
            scf_energy = float(line.split()[3]) * ry2ev
            SCF_ENERGY.append(scf_energy)
            break
        
    displacements_file = f"step_{istep}/10-displacements/displacements_Newton_method.dat"
    displacements = np.loadtxt(displacements_file, usecols=(1,2,3))
    mod_displacements = []
    for iatom in range(Natoms):
        mod_displacements.append(np.linalg.norm(displacements[3*iatom:3*(iatom+1)]))
    MAX_DISPLACEMENT.append(max(mod_displacements))
    MEAN_DISPLACEMENT.append(np.mean(mod_displacements))
    
    tot_energy = scf_energy + EXC_ENERGY[-1] * ry2ev
    TOT_ENERGY.append(tot_energy)
    
STEPS = np.array(STEPS)
MAX_ES_FORCE = np.array(MAX_ES_FORCE)
MEAN_ES_FORCE = np.array(MEAN_ES_FORCE)
EXC_ENERGY = np.array(EXC_ENERGY)
SCF_ENERGY = np.array(SCF_ENERGY) - SCF_ENERGY[0]
TOT_ENERGY = SCF_ENERGY + EXC_ENERGY
MAX_DISPLACEMENT = np.array(MAX_DISPLACEMENT)
MEAN_DISPLACEMENT = np.array(MEAN_DISPLACEMENT)
MAX_DFT_FORCE = np.array(MAX_DFT_FORCE)
MEAN_DFT_FORCE = np.array(MEAN_DFT_FORCE)
MAX_TOT_FORCE = np.array(MAX_TOT_FORCE)
MEAN_TOT_FORCE = np.array(MEAN_TOT_FORCE)

print('steps.shape = ', STEPS.shape)
print('MAX_ES_FORCE.shape = ', MAX_ES_FORCE.shape)
print('MEAN_ES_FORCE.shape = ', MEAN_ES_FORCE.shape)
print('MAX_DFT_FORCE.shape = ', MAX_DFT_FORCE.shape)
print('MEAN_DFT_FORCE.shape = ', MEAN_DFT_FORCE.shape)    
print('MAX_TOT_FORCE.shape = ', MAX_TOT_FORCE.shape)
print('MEAN_TOT_FORCE.shape = ', MEAN_TOT_FORCE.shape)

print('exc_energy.shape = ', EXC_ENERGY.shape)
print('scf_energy.shape = ', SCF_ENERGY.shape)
print('tot_energy.shape = ', TOT_ENERGY.shape)
print('max_displacement.shape = ', MAX_DISPLACEMENT.shape)
print('mean_displacement.shape = ', MEAN_DISPLACEMENT.shape)

# print(f"{'Step':>4} {'Max_ES_force':>8} {'Mean_ES_force':>8} {'MaxDFT_force':>8} {'MeanDFT_force':>8'} {'MaxTot_force':>8} {'MeanTot_force':>8} {'ESCF':>10} {'E_exc':>8} {'TotE':>10} {'MaxDisp':>10} {'MeanDisp':>10}")
print(f"{'Step':>4} {'Max_ES_force':>14} {'Mean_ES_force':>14} {'MaxDFT_force':>14} {'MeanDFT_force':>14} "
      f"{'MaxTot_force':>14} {'MeanTot_force':>14} {'ESCF':>10} {'E_exc':>8} {'TotE':>10} {'MaxDisp':>10} {'MeanDisp':>10}")

print("-" * 70)

for i in range(len(STEPS)):
    print(f"{STEPS[i]} ",
          f"{MAX_ES_FORCE[i]:8.4f} ",
          f"{MEAN_ES_FORCE[i]:8.4f} ",
          f"{MAX_DFT_FORCE[i]:8.4f} ",
          f"{MEAN_DFT_FORCE[i]:8.4f} ",
          f"{MAX_TOT_FORCE[i]:8.4f} ",
          f"{MEAN_TOT_FORCE[i]:8.4f} ",
          f"{SCF_ENERGY[i]:10.6f} ",
          f"{EXC_ENERGY[i]:.4f} ",
          f"{TOT_ENERGY[i]:.6f} ",
          f"{MAX_DISPLACEMENT[i]:.4f} ",
          f"{MEAN_DISPLACEMENT[i]:.4f}")