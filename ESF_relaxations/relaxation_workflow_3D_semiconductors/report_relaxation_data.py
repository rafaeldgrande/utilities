
import numpy as np

first_step = 0
last_step = 3

ry2ev = 13.605703976

STEPS = []
MAX_FORCE, MEAN_FORCE = [], []
EXC_ENERGY, SCF_ENERGY, TOT_ENERGY = [], [], []
MAX_DISPLACEMENT, MEAN_DISPLACEMENT = [], []

for istep in range(first_step, last_step+1):
    STEPS.append(float(istep))

    forces_file = f"step_{istep}/9-excited-state-forces/forces_cart.out_1_1"
    forces = np.loadtxt(forces_file, usecols=3, dtype=complex) # shape (3N)
    
    Natoms = len(forces) // 3
    mod_forces = []
    for iatom in range(Natoms):
        mod_forces.append(np.linalg.norm(forces[3*iatom:3*(iatom+1)]))
    MAX_FORCE.append(max(mod_forces))
    MEAN_FORCE.append(np.mean(mod_forces))
    
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
MAX_FORCE = np.array(MAX_FORCE)
MEAN_FORCE = np.array(MEAN_FORCE)
EXC_ENERGY = np.array(EXC_ENERGY)
SCF_ENERGY = np.array(SCF_ENERGY) - SCF_ENERGY[0]
TOT_ENERGY = SCF_ENERGY + EXC_ENERGY
MAX_DISPLACEMENT = np.array(MAX_DISPLACEMENT)
MEAN_DISPLACEMENT = np.array(MEAN_DISPLACEMENT)

print('steps.shape = ', STEPS.shape)
print('max_force.shape = ', MAX_FORCE.shape)
print('mean_force.shape = ', MEAN_FORCE.shape)
print('exc_energy.shape = ', EXC_ENERGY.shape)
print('scf_energy.shape = ', SCF_ENERGY.shape)
print('tot_energy.shape = ', TOT_ENERGY.shape)
print('max_displacement.shape = ', MAX_DISPLACEMENT.shape)
print('mean_displacement.shape = ', MEAN_DISPLACEMENT.shape)

print(f"{'Step':>4} {'MaxF':>8} {'MeanF':>8} {'SCF(dE)':>10} {'ExcE':>8} {'TotE':>10} {'MaxDisp':>10} {'MeanDisp':>10}")
print("-" * 70)

for i in range(len(STEPS)):
    print(f"{STEPS[i]} ",
          f"{MAX_FORCE[i]:8.4f} ",
          f"{MEAN_FORCE[i]:8.4f} ",
          f"{SCF_ENERGY[i]:10.6f} ",
          f"{EXC_ENERGY[i]:.4f} ",
          f"{TOT_ENERGY[i]:.6f} ",
          f"{MAX_DISPLACEMENT[i]:.4f} ",
          f"{MEAN_DISPLACEMENT[i]:.4f}")
    
    
    
    
