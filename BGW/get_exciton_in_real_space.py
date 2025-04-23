
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.fft import ifftn

WFN_file = 'WFN.h5'
eigenvectors_file = 'eigenvectors.h5'
q_vec = np.array([0, 0, 0])
Nval = 5
iexc_plot = 0

def split_coeffs(coeffs, gvecs, ngk):
    """
    Split the flattened coefficients array into a list of arrays, one per k-point.
    
    Parameters:
      coeffs (np.ndarray): Coefficients array with shape (n_bands, 1, total_G).
      ngk (array-like): 1D array of length nk, where each element is the number of G-vectors for that k-point.
      
    Returns:
      list: A list of length nk, where each element is an array of shape (n_bands, 1, ng),
            with ng being the number of G-vectors for that k-point.
    """
    list_coeffs = []
    list_gvecs = []
    start = 0
    for count in ngk:
        # Extract a slice for this k-point.
        # coeffs[:, 0, start:start+count] has shape (n_bands, count).
        sub_array = coeffs[:, 0, start:start+count]
        # Add back the singleton dimension in the second axis to get (n_bands, 1, count)
        sub_array = sub_array[:, np.newaxis, :]
        list_coeffs.append(sub_array)
        gvecs_sub_array = gvecs[start:start+count, :]
        list_gvecs.append(gvecs_sub_array)
        start += count
    return list_coeffs, list_gvecs

def wavefunction_for_k_point(coeffs, gvecs, grid_size):
    """Compute real-space wavefunction via inverse FFT."""
    n_bands, _, _ = coeffs.shape
    psi_k_point = np.zeros((n_bands, *grid_size), dtype=complex)

    for n in range(n_bands):
        # Place coefficients in a 3D reciprocal-space grid
        grid = np.zeros(grid_size, dtype=complex)
        for i, G in enumerate(gvecs):
            grid[tuple(G)] = coeffs[n, 0, i]

        # Inverse FFT to get real-space wavefunction
        psi_k_point[n] = ifftn(grid)

    return psi_k_point

def wavefunction_real_space(COEFFS, GVECS, grid_size):
    n_bands, _, _ = coeffs.shape
    n_kpoints = len(COEFFS)
    psi_real_space = np.zeros((n_kpoints, n_bands, *grid_size), dtype=complex)
    
    for ik in range(n_kpoints):
        psi_real_space[ik] = wavefunction_for_k_point(COEFFS[ik], GVECS[ik], grid_size)

    return psi_real_space

import numpy as np
import h5py

def exciton_wavefunction(A_vck, psi_real, kpoints, q_vec, grid_size, Nval, output_filename="exciton_wavefunction.h5"):
    """
    Compute and store the 6D exciton wavefunction Psi(x_e, y_e, z_e, x_h, y_h, z_h) in an HDF5 file.
    
    Parameters:
      A_vck: Exciton coefficients (nk, nc, nv)
      psi_real: Real-space wavefunctions (nk, Nbnds, Nx, Ny, Nz)
      kpoints: K-point coordinates (nk, 3)
      q_vec: Momentum transfer vector (3,)
      grid_size: Tuple (Nx, Ny, Nz)
      Nval: Number of valence bands (highest valence band index is Nval-1)
      output_filename: Name of the output HDF5 file (default: "exciton_wavefunction.h5")
      
    Returns:
      None (writes Psi to HDF5 file)
    """
    nk, nc, nv = A_vck.shape
    Nx, Ny, Nz = grid_size

    # Create real-space grid (assumed to span [0,1) in each dimension)
    x = np.linspace(0, 1, Nx, endpoint=False)
    y = np.linspace(0, 1, Ny, endpoint=False)
    z = np.linspace(0, 1, Nz, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    r_grid = np.stack((X, Y, Z), axis=-1)  # shape: (Nx, Ny, Nz, 3)
    r_prime_grid = r_grid  # assume same grid for the hole

    # Create an HDF5 file to store the exciton wavefunction
    with h5py.File(output_filename, "w") as h5f:
        dset = h5f.create_dataset("Psi", (Nx, Ny, Nz, Nx, Ny, Nz), dtype=np.complex128)

        # Progress tracking
        total_iterations = nk * nc * nv
        update_interval = max(total_iterations // 100, 1)
        counter = 0

        for k in range(nk):
            k_vec = np.squeeze(kpoints[k])
            phase_e = np.exp(1j * (r_grid[...,0]*k_vec[0] +
                                   r_grid[...,1]*k_vec[1] +
                                   r_grid[...,2]*k_vec[2]))
            phase_h = np.exp(-1j * (r_prime_grid[...,0]*(k_vec[0]+q_vec[0]) +
                                    r_prime_grid[...,1]*(k_vec[1]+q_vec[1]) +
                                    r_prime_grid[...,2]*(k_vec[2]+q_vec[2])))
            for c in range(nc):
                idx_c = Nval + c  # conduction band index in psi_real
                phi_ck = psi_real[k, idx_c, :, :, :]  # electron wavefunction
                for v in range(nv):
                    idx_v = Nval - 1 - v  # valence band index in psi_real
                    phi_vk = psi_real[k, idx_v, :, :, :]  # hole wavefunction
                    electron_part = phi_ck * phase_e
                    hole_part = phi_vk * phase_h
                    contrib = A_vck[k, c, v] * electron_part[:, :, :, np.newaxis, np.newaxis, np.newaxis] * \
                              np.conj(hole_part)[np.newaxis, np.newaxis, np.newaxis, :, :, :]
                    
                    # Update Psi directly in HDF5 to save RAM
                    dset[:] += contrib
                    
                    # Progress update
                    counter += 1
                    if counter % update_interval == 0:
                        percent_done = (counter / total_iterations) * 100
                        print(f"Progress: {percent_done:.2f}% done", end='\r')

        print("\nFinished computing and saving Psi to HDF5 file.")



# considering nspin = 1

# Open the HDF5 file
with h5py.File(eigenvectors_file, "r") as f:
    A_vck = f['exciton_data/eigenvectors'][()]   # Adjust the key name based on actual structure
A_vck = A_vck[0, :, :, :, :, 0, 0] + 1j * A_vck[0, :, :, :, :, 0, 1]

# loading data from WFN file
with h5py.File(WFN_file, "r") as f:
    coeffs_temp = f["/wfns/coeffs"][:]  # Check actual dataset name
    gvecs = f["/wfns/gvecs"][:]  # Reciprocal lattice vectors
    kpoints = f["/mf_header/kpoints/rk"][:]  # shape (3, nrk)
    kpoints = kpoints.T  # shape (nrk, 3)
    FFTgrid = f["/mf_header/gspace/FFTgrid"][:]
    ngk = f["/mf_header/kpoints/ngk"][:]
    
grid_size = tuple(FFTgrid)
coeffs = coeffs_temp[:, :, :, 0] + 1j * coeffs_temp[:, :, :, 1]

# Now split:
COEFFS, GVECS = split_coeffs(coeffs, gvecs, ngk)

# running
psi_real = wavefunction_real_space(COEFFS, GVECS, grid_size)

Exciton_wavefunc = exciton_wavefunction(A_vck[iexc_plot], psi_real, kpoints, q_vec, grid_size, Nval)
