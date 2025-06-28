
# 04 parabands, 05-epsilon, 06-sigma, 07-kernel, 08-absorption

nc_coarse, nv_coarse = 36, 19
nc_fine, nv_fine = 36, 19
Nkx, Nky, Nkz = 1, 1, 1
kgrid = Nkx, Nky, Nkz, 0, 0, 0
# for the case of super cell with 54 atoms of Si, there are 4*54 = 216 electrons = 108 bands
Nval = 108 

def generate_kpoints(epsilon_or_sigma, kgrid):
    Nkx, Nky, Nkz, qx, qy, qz = kgrid
    
    if epsilon_or_sigma == 'epsilon':
        text = 'begin qpoints \n'
    else:
        text = 'begin kpoints \n'

    if epsilon_or_sigma == 'epsilon':
        for ik in range(Nkx):
            for jk in range(Nky):
                for kk in range(Nkz):
                    if ik == 0 and jk == 0 and kk == 0:
                        text += f"0.001 0.000 0.000  1.0 1\n"
                    else:
                        text += f"{(ik/Nkx+qx):.8f}  {(jk/Nky+qy):.8f}  {(kk/Nkz+qz):.8f}  1.0    0\n"
    else: 
        for ik in range(Nkx):
            for jk in range(Nky):
                for kk in range(Nkz):
                    text += f"{(ik/Nkx+qx):.8f}  {(jk/Nky+qy):.8f}  {(kk/Nkz+qz):.8f}  1.0\n"
                    
    text += 'end\n'
    return text

parabands_text = f"""input_wfn_file wfn.complex
output_wfn_file WFN.h5
vsc_file VSC
vkb_file VKB
verbosity 2

use_pseudobands
protected_cond_bands {nc_coarse + 10}
accumulation_window 0.02"""

epsilon_text = f"""

use_wfn_hdf5
verbosity 3

epsilon_cutoff 10.0
degeneracy_check_override

"""
epsilon_text += generate_kpoints('epsilon', kgrid)

sigma_text = f"""
verbosity 3

use_wfn_hdf5

degeneracy_check_override

band_index_min {Nval - nv_coarse - 1}
band_index_max {Nval + nc_coarse + 1}

dont_use_vxcdat
screening_semiconductor

"""
sigma_text += generate_kpoints('sigma', kgrid)

kernel_text = f"""
use_wfn_hdf5

number_val_bands {nv_coarse}
number_cond_bands {nc_coarse}

screening_semiconductor

no_symmetries_coarse_grid
"""

absorption_text = f"""
degeneracy_check_override

number_val_bands_fine {nv_fine}
number_val_bands_coarse {nv_coarse}

number_cond_bands_fine {nc_fine}
number_cond_bands_coarse {nc_coarse}

no_symmetries_fine_grid
no_symmetries_shifted_grid
no_symmetries_coarse_grid

use_wfn_hdf5

eqp_co_corrections

diagonalization

screening_semiconductor

use_momentum

energy_resolution 0.05
gaussian_broadening

write_eigenvectors 100"""

arq = open('4-parabands/parabands.inp', 'w')
arq.write(parabands_text)
arq.close()
print(f"Finished writing parabands.inp")

arq = open('5-epsilon/epsilon.inp', 'w')
arq.write(epsilon_text)
arq.close()
print(f"Finished writing epsilon.inp")

arq = open('6-sigma/sigma.inp', 'w')
arq.write(sigma_text)
arq.close()
print(f"Finished writing sigma.inp")

arq = open('7-kernel/kernel.inp', 'w')
arq.write(kernel_text)
arq.close()
print(f"Finished writing kernel.inp")

arq = open('8-absorption/absorption.inp', 'w')
arq.write(absorption_text)
arq.close()
print(f"Finished writing absorption.inp")