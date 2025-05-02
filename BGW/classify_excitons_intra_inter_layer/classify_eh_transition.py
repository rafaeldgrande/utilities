import numpy as np
import h5py
from scipy.fft import ifftn
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

bohr2Ang = 0.529177

# WFN_file = "/Users/rdelgrande/work/Projects/WSe2/WFN.h5"
WFN_file = "/Users/rdelgrande/work/TEMP_DATA/teste_compute_rho_z_2d_bilayer/WFN.h5"
nc = 13
nv = 5

z1 = -0.850290
z2 = 5.650293
z_div = (z1 + z2) / 2 

def split_coeffs(coeffs, gvecs, ngk):
    print("Splitting coefficients and G-vectors for each k-point...")
    list_coeffs = []
    list_gvecs = []
    start = 0
    for count in ngk:
        sub_array = coeffs[:, 0, start:start+count]
        sub_array = sub_array[:, np.newaxis, :]
        list_coeffs.append(sub_array)
        gvecs_sub_array = gvecs[start:start+count, :]
        list_gvecs.append(gvecs_sub_array)
        start += count
    print("Splitting completed.")
    return list_coeffs, list_gvecs

def wavefunction_for_k_point(coeffs, gvecs, grid_size):
    print("Computing real-space wavefunction for a k-point...")
    n_bands, _, _ = coeffs.shape
    psi_k_point = np.zeros((n_bands, *grid_size), dtype=complex)

    for n in range(n_bands):
        grid = np.zeros(grid_size, dtype=complex)
        for i, G in enumerate(gvecs):
            grid[tuple(G)] = coeffs[n, 0, i]
        psi_k_point[n] = ifftn(grid)

    print("Wavefunction computation for k-point completed.")
    return psi_k_point

def wavefunction_real_space(COEFFS, GVECS, grid_size):
    print("Computing real-space wavefunction for all k-points...")
    n_bands, _, _ = coeffs.shape
    n_kpoints = len(COEFFS)
    psi_real_space = np.zeros((n_kpoints, n_bands, *grid_size), dtype=complex)
    
    for ik in range(n_kpoints):
        print(f"Processing k-point {ik+1}/{n_kpoints}...")
        psi_real_space[ik] = wavefunction_for_k_point(COEFFS[ik], GVECS[ik], grid_size)

    print("Real-space wavefunction computation completed.")
    return psi_real_space

# loading data from WFN file
print("Loading data from WFN file...")
with h5py.File(WFN_file, "r") as f:
    coeffs_temp = f["/wfns/coeffs"][:]
    gvecs = f["/wfns/gvecs"][:]
    kpoints = f["/mf_header/kpoints/rk"][:]
    kpoints = kpoints.T
    FFTgrid = f["/mf_header/gspace/FFTgrid"][:]
    ngk = f["/mf_header/kpoints/ngk"][:]
    ifmax = f["/mf_header/kpoints/ifmax"][:]
    lat_vecs = f["/mf_header/crystal/avec"][:]
    alat = f["/mf_header/crystal/alat"][()]

print("Data loading completed.")

Nval = ifmax[0, 0]
grid_size = tuple(FFTgrid)
min_idx = Nval - nv
max_idx = Nval + nc
coeffs = coeffs_temp[min_idx:max_idx, :, :, 0] + 1j * coeffs_temp[min_idx:max_idx, :, :, 1]

print("Splitting coefficients and G-vectors...")
COEFFS, GVECS = split_coeffs(coeffs, gvecs, ngk)

print("Computing real-space wavefunction...")
psi_real = wavefunction_real_space(COEFFS, GVECS, grid_size)

print("Calculating electronic density...")
rho_z = psi_real[:, :, :, :, :] * np.conj(psi_real[:, :, :, :, :])
rho_z = np.real(rho_z.sum(axis=(2, 3)))


norm_factors = rho_z.sum(axis=2, keepdims=True)
rho_z /= norm_factors

print("rho_z shape:", rho_z.shape)

Lz = alat * bohr2Ang * lat_vecs[2, 2]
print("Lz:", Lz)

Lz_plot = np.linspace(0, Lz, grid_size[2])
Lz_plot2 = np.linspace(-Lz/2, Lz/2, grid_size[2])

idx_Lz_eq_0 = np.argmin(np.abs(Lz_plot2 - 0))

idx_shift = np.argmin(np.abs(Lz_plot - z_div))
print("idx_shift:", idx_shift)

idx_cut = np.argmin(np.abs(Lz_plot - (z_div + Lz/2)))
print("idx_cut:", idx_cut)

nz = Lz_plot.shape[0]
nz1 = idx_shift
nz2 = idx_cut - idx_shift
nz3 = nz - (nz1 + nz2)
print("nz1, nz2, nz3:", nz1, nz2, nz3, nz1+nz2+nz3)

print("Reorganizing rho_z...")
rho_new = np.zeros((rho_z.shape))
rho_new[:, :, :nz3] = rho_z[:, :, idx_cut:]
rho_new[:, :, nz3:nz1+nz3] = rho_z[:, :, :idx_shift]
rho_new[:, :, nz1+nz3:] = rho_z[:, :, idx_shift:idx_cut]

nk, nbnds, _ = rho_new.shape

print("Calculating projections...")
projection = np.zeros((rho_new.shape[0], rho_new.shape[1], 2))

projection[:, :, 0] = np.sum(rho_new[:, :, :idx_Lz_eq_0], axis=2)
projection[:, :, 1] = np.sum(rho_new[:, :, idx_Lz_eq_0:], axis=2)

print("Saving wavefunction projections to PDF...")
with PdfPages("wavefunction_projections.pdf") as pdf:
    for ik in range(nk):
        for ib in range(nbnds):
            plt.figure()
            plt.plot(Lz_plot2, rho_new[ik, ib, :])
            plt.title(f"Band {ib+1}, k-point {ik+1}, nL = {projection[ik, ib, 0]:.3f}, nR = {projection[ik, ib, 1]:.3f}")
            plt.axvline(z1 - z_div, color='k', linestyle='--')
            plt.axvline(z2 - z_div, color='k', linestyle='--')
            pdf.savefig()
            plt.close()

print("Calculating x_kcv...")
x_kcv = np.zeros((nk, nc, nv))

for ik in range(nk):
    for ic in range(nc):
        for iv in range(nv):
            cLvL = projection[ik, nv+ic, 0] * projection[ik, nv-iv, 0]
            cRvR = projection[ik, nv+ic, 1] * projection[ik, nv-iv, 1]
            cLvR = projection[ik, nv+ic, 0] * projection[ik, nv-iv, 1]
            cRvL = projection[ik, nv+ic, 1] * projection[ik, nv-iv, 0]
            
            x_kcv[ik, ic, iv] = (cLvL + cRvR - cLvR - cRvL) / (cLvL + cRvR + cLvR + cRvL)

print("Saving x_kcv to HDF5 file...")
with h5py.File("x_kcv.h5", "w") as f:
    f.create_dataset("x_kcv", data=x_kcv)

print('Finished writing x_kcv.h5')
