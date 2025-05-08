

''' This code calculates rho(ze, zh) for a given set of exciton coefficients and wavefunctions.
    where rho(ze, zh) = integral dx_e dy_e dx_h dy_h |Psi(x_e, y_e, z_e; x_h, y_h, z_h)|^2
    
    it reads exciton coefficients from eigenvectors.h5 and wavefunctions from WFN.h5
    and saves the result in RHO_ZE_ZH_data.h5'''

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.fft import ifftn
import time
from multiprocessing import Pool
import argparse
import psutil

# pip install numpy matplotlib h5py scipy psutil

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()

parser.add_argument('--run_parallel', type=str2bool, default=True, help='Run in parallel or serial mode')
parser.add_argument('--num_processes', type=int, default=1, help='Number of parallel processes to use')

parser.add_argument('--verbose', type=int, default=0, help='Enable verbose output')

parser.add_argument('--WFN_file', type=str, default='WFN.h5', help='Path to WFN file')
parser.add_argument('--eigenvectors_file', type=str, default='eigenvectors.h5', help='Path to eigenvectors file')

# consider Psi(ze, zh) = 0 for zmin < ze < zmax or zmin < zh < zmax
# useful for 2D materials with large vacuum distance
parser.add_argument('--zmin_set_zero', type=float, default=8.5, help='Minimum z value to set to zero')
parser.add_argument('--zmax_set_zero', type=float, default=30 - 1.5, help='Maximum z value to set to zero')
parser.add_argument('--Lz', type=float, default=30, help='Length in z direction')

parser.add_argument('--limitGvecs', type=str2bool, default=False, help='Limit G-vectors for testing')
parser.add_argument('--maxGs', type=int, default=100, help='Maximum number of G-vectors for testing purposes')

parser.add_argument('--limitGrid', type=str2bool, default=False, help='Limit grid size for testing')
parser.add_argument('--Ngrid_teste', type=int, default=4, help='Grid size for testing purposes')

parser.add_argument('--i_exc', type=int, default=0, help='Index of exciton to compute')
parser.add_argument('--limit_BSE_sum_up_to_value', type=float, default=1.0, help='Limit for BSE sum up to value')
# parser.add_argument('--q_vec', type=float, nargs=3, default=[0, 0, 0], help='Momentum transfer vector (q_x, q_y, q_z)')

# usage: python exciton_density_as_function_of_ze_zh.py --run_parallel True --num_processes 16 --verbose 1 --WFN_file WFN.h5 --eigenvectors_file eigenvectors.h5 --zmin_set_zero 8.5 --zmax_set_zero 30 --Lz 30 --limitGvecs True --maxGs 100 --limitGrid False --Ngrid_teste 4 --i_exc 0

args = parser.parse_args()

num_processes = args.num_processes
verbose = args.verbose
run_parallel = args.run_parallel
WFN_file = args.WFN_file
eigenvectors_file = args.eigenvectors_file
zmin_set_zero = args.zmin_set_zero
zmax_set_zero = args.zmax_set_zero
Lz = args.Lz
i_exc = args.i_exc
# q_vec = np.array(args.q_vec)
limitGvecs = args.limitGvecs
maxGs = args.maxGs
limitGrid = args.limitGrid
Ngrid_teste = args.Ngrid_teste
limit_BSE_sum_up_to_value = args.limit_BSE_sum_up_to_value
# Start the timer
start_time = time.time()

q_vec = np.array([0, 0, 0])  # q shift between electron and hole

precision_complex = np.complex64 # or use np.complex128 for more precision


def set_z_zero(z):
    if zmin_set_zero < z < zmax_set_zero:
        return False
    else:
        return True

if run_parallel:
    print(f"Running in parallel with {num_processes} processes.")
else:
    print("Running in serial mode.")

if limitGvecs:
    print("Limiting G-vectors to", maxGs)
else:
    print("Not limiting G-vectors")

if limitGrid:
    print("Limiting grid size to", Ngrid_teste)
else:
    print("Not limiting grid size")
    
def top_n_indexes_all(array, limit_BSE_sum_up_to_value):
    """
    Identify the top N indexes in a 3D array based on their values, and return
    the indexes up to the point where the cumulative sum of squared values exceeds
    a specified limit.
    Parameters:
    -----------
    array : numpy.ndarray
        A 3D numpy array containing the values to be analyzed.
    limit_BSE_sum_up_to_value : float
        The threshold for the cumulative sum of squared values. The function
        will stop collecting indexes once this threshold is exceeded.
    Returns:
    --------
    list of tuple
        A list of tuples representing the indexes of the top values in the array,
        ordered by descending value. The list is truncated to include only the
        indexes needed to reach the specified cumulative sum threshold.
    """
    
    # Flatten the array
    flat_array = array.flatten()
    
    # array size
    N = len(flat_array)
    
    # Get the indexes of the top N values in the flattened array
    flat_indexes = np.argpartition(flat_array, -N)[-N:]
    
    # Sort these indexes by the values they point to, in descending order
    sorted_indexes = flat_indexes[np.argsort(-flat_array[flat_indexes])]
    
    # Convert the 1D indexes back to 3D indexes
    top_indexes = np.unravel_index(sorted_indexes, array.shape)
    
    # Combine the indexes into a list of tuples
    top_indexes = list(zip(*top_indexes))
    
    # now checking how many values we need to store
    counter_indexes = 0
    sum_abs_Akcv2 = 0
    for index in top_indexes:
        counter_indexes += 1 
        sum_abs_Akcv2 += array[index[0], index[1], index[2]]**2
        if sum_abs_Akcv2 > limit_BSE_sum_up_to_value:
            break

    return top_indexes[:counter_indexes]
    
def check_available_memory():
    # Get the total system memory (RAM) in bytes
    total_memory = psutil.virtual_memory().total
    
    # Get the available system memory (RAM) in bytes
    available_memory = psutil.virtual_memory().available

    # Convert bytes to GB for easier reading
    total_memory_gb = total_memory / (1024 ** 3)  # Convert bytes to GB
    available_memory_gb = available_memory / (1024 ** 3)  # Convert bytes to GB

    print(f"Total memory: {total_memory_gb:.2f} GB")
    print(f"Available memory: {available_memory_gb:.2f} GB")

    
def calculate_memory_Psi(Nx, Ny, num_processes, precision_complex):
    # Calculate the size of one Psi in bytes
    size_one_psi = Nx * Ny * Nx * Ny * np.dtype(precision_complex).itemsize  # Item size in bytes

    # Total memory required for all processes
    total_memory = size_one_psi * num_processes  # Total memory in bytes

    # Convert to MB and GB for easier reading
    total_memory_mb = total_memory / (1024 ** 2)  # Convert bytes to MB
    total_memory_gb = total_memory / (1024 ** 3)  # Convert bytes to GB

    # Print memory usage
    print(f"Each exciton wavefunction calculation requires {size_one_psi / (1024 ** 2):.2f} MB")
    print(f"Total memory required for {num_processes} processes: {total_memory_mb:.2f} MB ({total_memory_gb:.2f} GB)")

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
        sub_array = sub_array[:, np.newaxis, :] # shape (n_bands, 1, count)
        gvecs_sub_array = gvecs[start:start+count, :] # shape (count, 3)
        
        if limitGvecs:
            # Limit the number of G-vectors for testing purposes
            sub_array = sub_array[:, :, :maxGs]
            gvecs_sub_array = gvecs_sub_array[:maxGs, :]
            
        list_coeffs.append(sub_array)
        list_gvecs.append(gvecs_sub_array)
        
        start += count
    return list_coeffs, list_gvecs

def wavefunction_for_k_point(coeffs, gvecs, grid_size):
    """Compute real-space wavefunction via inverse FFT."""
    n_bands, _, _ = coeffs.shape
    psi_k_point = np.zeros((n_bands, *grid_size), dtype=precision_complex)

    for n in range(n_bands):
        # Place coefficients in a 3D reciprocal-space grid
        grid = np.zeros(grid_size, dtype=precision_complex)
        for i, G in enumerate(gvecs):
            grid[tuple(G)] = coeffs[n, 0, i]

        # Inverse FFT to get real-space wavefunction
        psi_k_point[n] = ifftn(grid)

    return psi_k_point

def wavefunction_real_space(COEFFS, GVECS, grid_size):
    """
    Computes the wavefunction in real space for a given set of coefficients and G-vectors.
    Parameters:
        COEFFS (list): A list of length nk, where each element is an array of shape (n_bands, 1, ng),
            with ng being the number of G-vectors for that k-point.
        GVECS (numpy.ndarray): A 2D array of shape (n_kpoints, n_coeffs, 3) containing the 
                               G-vectors corresponding to the coefficients for each k-point.
        grid_size (tuple): A tuple of three integers specifying the size of the real-space grid 
                           (nx, ny, nz).
    Returns:
        numpy.ndarray: A 4D array of shape (n_kpoints, n_bands, nx, ny, nz) containing the 
                       wavefunctions in real space, represented as complex numbers.
    """
    n_bands, _, _ = COEFFS[0].shape
    psi_real_space = np.zeros((Nk, n_bands, *grid_size), dtype=precision_complex)
    
    for ik in range(Nk):
        psi_real_space[ik] = wavefunction_for_k_point(COEFFS[ik], GVECS[ik], grid_size)

    return psi_real_space

def exciton_wavefunction_fixed_ze_zh(A_vck, psi_real, kpoints, q_vec, grid_size, Nval, ze_idx, zh_idx):
    """
    Compute Psi(x_e, y_e, x_h, y_h) for fixed z_e and z_h.

    Parameters:
      A_vck: Exciton coefficients with shape (nk, nc, nv)
      psi_real: Real-space wavefunctions with shape (nk, Nbnds, Nx, Ny, Nz)
      kpoints: Array of k-point coordinates with shape (nk, 3)
      q_vec:   Momentum transfer vector (3,)
      grid_size: Tuple (Nx, Ny, Nz) for the real-space grid
      Nval:    Number of valence bands (highest valence band index is Nval-1)
      ze_idx:  Index of fixed z_e in the grid (0 ≤ ze_idx < Nz)
      zh_idx:  Index of fixed z_h in the grid (0 ≤ zh_idx < Nz)

    Returns:
      Psi: 4D exciton wavefunction with shape (Nx, Ny, Nx, Ny)
    """

    # Initialize Psi(x_e, y_e, x_h, y_h)
    
    if verbose > 1:
        print(f"Computing Psi for ze_idx={ze_idx}, zh_idx={zh_idx}...")
        time_start_compute_psi = time.time()
        
        delta_t_phase_calcs = 0.0
        delta_t_kcv_loop = 0.0
    
    Psi = np.zeros((Nx, Ny, Nx, Ny), dtype=precision_complex)
    

    # Iterate over k-points, conduction and valence bands
    for k in range(Nk):
        
        time_phase_start = time.time()
        phase_e = np.exp(1j * (r_grid_xy[...,0] * kpoints[0, k] + r_grid_xy[...,1] * kpoints[1, k]))
        phase_h = np.exp(-1j * (r_grid_xy[...,0] * (kpoints[0, k] + q_vec[0]) +
                                r_grid_xy[...,1] * (kpoints[1, k] + q_vec[1])))
        time_phase_end = time.time()
        if verbose > 1:
            delta_t_phase_calcs += time_phase_end - time_phase_start
        
        if verbose > 1:    
            time_kcv_start = time.time()
        
        psi_k = psi_real[k] # Shape: (Nbands, Nx, Ny, Nz)

        # for c in range(Nc):
        #     idx_c = Nv + c  # Conduction band index
        #     phi_ck = psi_k[idx_c, :, :, ze_idx]  # Electron wavefunction at fixed z_e
            
        #     for v in range(Nv):
        #         idx_v = Nv - 1 - v  # Valence band index
        #         phi_vk = psi_k[idx_v, :, :, zh_idx]  # Hole wavefunction at fixed z_h

        #         # Compute the wavefunction contributions
        #         electron_part = phi_ck * phase_e
        #         hole_part = phi_vk * phase_h

        #         # Outer product: electron part on (x_e, y_e), hole part on (x_h, y_h)
        #         # contrib = A_vck[k, c, v] * np.einsum("ij,kl->ijkl", electron_part, np.conj(hole_part))
        #         contrib = A_vck[k, c, v] * (electron_part[:, :, None, None] * np.conj(hole_part)[None, None, :, :])

        #         Psi += contrib
        
        phi_ck_all = psi_k[Nv:Nv+Nc, :, :, ze_idx] * phase_e[None, :, :]  # shape: (Nc, Nx, Ny)
        phi_vk_all = psi_k[Nv-1::-1, :, :, zh_idx] * phase_h[None, :, :]
        
        if abs(limit_BSE_sum_up_to_value - 1.0) < 1e-10: # here limit_BSE_sum_up_to_value = 1.0
            for c in range(Nc):
                for v in range(Nv):
                    contrib = A_vck[k, c, v] * (phi_ck_all[c][:, :, None, None] * np.conj(phi_vk_all[v])[None, None, :, :])
                    Psi += contrib
        else:
            for c in range(Nc):
                for v in range(Nv):
                    if (k, c, v) not in top_kcv_set:
                        continue
                    contrib = A_vck[k, c, v] * (phi_ck_all[c][:, :, None, None] * np.conj(phi_vk_all[v])[None, None, :, :])
                    Psi += contrib
         
        if verbose > 1:       
            time_kcv_end = time.time()
            delta_t_kcv_loop += time_kcv_end - time_kcv_start

    if verbose > 1:
        time_end_compute_psi = time.time()
        print(f"Time taken to compute Psi for ze_idx={ze_idx}, zh_idx={zh_idx}: {time_end_compute_psi - time_start_compute_psi:.2f} seconds")
        print(f"Time taken for phase calculations: {delta_t_phase_calcs:.2f} seconds")
        print(f"Time taken for kcv loop: {delta_t_kcv_loop:.2f} seconds")

    return Psi


def compute_exciton_wavefunction(params):
    """
    Compute the exciton wavefunction for given parameters.

    This function processes a specific combination of indices and computes
    the exciton wavefunction for fixed electron and hole positions on a grid.
    It returns the indices and the computed sum of the absolute values of the
    wavefunction.

    Args:
        params (tuple): A tuple containing the following elements:
            i (int): Index corresponding to the first dimension.
            j (int): Index corresponding to the second dimension.
            ze_idx (int): Index for the fixed electron position.
            zh_idx (int): Index for the fixed hole position.

    Returns:
        tuple: A tuple containing:
            i (int): The first index from the input.
            j (int): The second index from the input.
            float: The sum of the absolute values of the computed wavefunction.
    """
    
    time_start_here = time.time()
    i, j, ze_idx, zh_idx = params
    # print(f"Processing (ze_idx={ze_idx}, zh_idx={zh_idx}) at index ({i}, {j})...")
    Psi_2d = exciton_wavefunction_fixed_ze_zh(A_vck[i_exc], phi, kpoints, q_vec, grid_size, Nval, ze_idx, zh_idx)
    Psi_2d_sum_sq = np.sum(np.abs(Psi_2d)**2)
    time_end_here = time.time()
    if verbose > 1:
        print(f"Time taken for (ze_idx={ze_idx}, zh_idx={zh_idx}) at index ({i}, {j}): {time_end_here - time_start_here:.2f} seconds")
    return i, j, Psi_2d_sum_sq  # Return index and computed value

def calc_parallel(num_processes, RHO_ZE_ZH, verbose):
    chunck_size = len(params_list) // num_processes + 1 
    params_chunks = np.array_split(params_list, chunck_size)
    total_processed = 0
    
    print(f"Dividing {len(params_list)} ze,zh pairs into {len(params_chunks)} chunks for parallel processing.")
    print(f"Each chunk will handle {len(params_chunks)} ze,zh pairs at most")

    with Pool(processes=num_processes) as pool:
        for param_chunk in params_chunks:
            # Process the chunk in parallel
            partial_results = pool.map(compute_exciton_wavefunction, param_chunk)

            # Update the output array
            for i, j, abs_psi2 in partial_results:
                RHO_ZE_ZH[i, j] = abs_psi2

            # Update and report number of terms calculated
            total_processed += len(param_chunk)
            print(f"Computed {total_processed}/{len(params_list)} terms")
            
    return RHO_ZE_ZH

def calc_serial(RHO_ZE_ZH, verbose=False):
    
    total_iterations = len(params_list)
    iteration_count = 0
    if verbose > 0:
        print(f"Total iterations to be calculated: {total_iterations}")
    
    for params in params_list:
        i, j, abs_psi2 = compute_exciton_wavefunction(params)
        RHO_ZE_ZH[i, j] = abs_psi2
        
        iteration_count += 1
        
        # Print progress every 100 iterations
        if verbose > 0:
            if iteration_count % 100 == 0:
                print(f"Completed {iteration_count} of {total_iterations} = {iteration_count*100/total_iterations:.2f}%  iterations.")
        
    return RHO_ZE_ZH

# Open the HDF5 file
print("Opening eigenvectors file to read exciton coefficients...")
with h5py.File(eigenvectors_file, "r") as f:
    A_vck = f['exciton_data/eigenvectors'][()]   # Adjust the key name based on actual structure
A_vck = A_vck[0, :, :, :, :, 0, 0] + 1j * A_vck[0, :, :, :, :, 0, 1]
Nexc, Nk, Nc, Nv = A_vck.shape
print("Exciton coefficients loaded.")
print("A_vck shape:", A_vck.shape)
print("Nk, Nc, Nv:", Nk, Nc, Nv)
print("Nexc:", Nexc)

# getting important transitions in Akcv
top_kcv = top_n_indexes_all(np.abs(A_vck[i_exc]), limit_BSE_sum_up_to_value)
top_kcv_set = set(top_kcv)  # convert to set for O(1) lookups

if limit_BSE_sum_up_to_value < 1.0:
    print(f"Limiting exciton coefficients to {len(top_kcv)} (instead of {Nk*Nc*Nv}) transitions based limit_BSE_sum_up_to_value = {limit_BSE_sum_up_to_value}.")

# Reading WFN file
print("Opening WFN file to read wavefunction data...")
with h5py.File(WFN_file, "r") as f:
    # coeffs_temp = f["/wfns/coeffs"][:]  # Check actual dataset name
    gvecs = f["/wfns/gvecs"][:]  # Reciprocal lattice vectors
    kpoints = f["/mf_header/kpoints/rk"][:]  # shape (3, nrk)
    kpoints = kpoints.T  # shape (nrk, 3)
    FFTgrid = f["/mf_header/gspace/FFTgrid"][:]
    ngk = f["/mf_header/kpoints/ngk"][:] # shape (nrk,)
    ifmax = f["/mf_header/kpoints/ifmax"][:]
print("Wavefunction data loaded.")

# Extract grid size and process coefficients
Nval = ifmax[0, 0]
grid_size = tuple(FFTgrid)
Nx, Ny, Nz = grid_size
min_idx = Nval - Nv
max_idx = Nval + Nc

# Generate real-space grids for electron and hole
x = np.linspace(0, 1, Nx, endpoint=False)
y = np.linspace(0, 1, Ny, endpoint=False)
X, Y = np.meshgrid(x, y, indexing='ij')
r_grid_xy = np.stack((X, Y), axis=-1)  # Shape: (Nx, Ny, 2)

with h5py.File(WFN_file, "r") as f:
    coeffs = f["/wfns/coeffs"][min_idx:max_idx, :, :, 0] + 1j * f["/wfns/coeffs"][min_idx:max_idx, :, :, 1]
    # shape of coeffs is (Nval+Nc-Nv, nk, 1, ng)
    
Nbands = coeffs.shape[0]

# coeffs = coeffs_temp[min_idx:max_idx, :, :, 0] + 1j * coeffs_temp[min_idx:max_idx, :, :, 1]
print(f"Grid size: {grid_size}")

dz = Lz / Nz
z_real_value = np.arange(0, Lz, dz)

# Create lists of ZE_idx and ZH_idx
print("Generating ZE_idx and ZH_idx arrays...")

if limitGrid:
    Nelectron = Ngrid_teste
    Nhole = Ngrid_teste
else:
    Nelectron = Nz
    Nhole = Nz

ZE_idx = np.linspace(0, Nz-1, Nelectron, dtype=int)
ZH_idx = np.linspace(0, Nz-1, Nhole, dtype=int)
print(f"ZE_idx: {ZE_idx}")
print(f"ZH_idx: {ZH_idx}")

RHO_ZE_ZH = np.zeros((len(ZE_idx), len(ZH_idx)))  # Preallocate array

# Prepare arguments for parallel processing
# print("Preparing parameters for parallel computation...")
# params_list = [(i, j, ze_idx, zh_idx) for i, ze_idx in enumerate(ZE_idx) for j, zh_idx in enumerate(ZH_idx)]

params_list = []
for i, ze_idx in enumerate(ZE_idx):
    if not set_z_zero(z_real_value[ze_idx]):
        continue
    for j, zh_idx in enumerate(ZH_idx):
        if not set_z_zero(z_real_value[zh_idx]):
            continue
        params_list.append((i, j, ze_idx, zh_idx))
        
print(f"ZE_idx shape: {ZE_idx.shape}")
print(f"Total number of abs(Psi(ze, zh))**2 we will compute: {len(params_list)}")

print("Splitting coefficients by k-points...")
COEFFS, GVECS = split_coeffs(coeffs, gvecs, ngk)
del coeffs, gvecs  # Free up memory

# COEFSS shape is (nk, n_bands, 1, ng)
# print max value of ng in COEFFS
print(f"Max number of G-vectors: {max(ngk)}")

if verbose > 0:
    print("len COEFFS:", len(COEFFS))
    print("COEFFS[0].shape:", COEFFS[0].shape)
    print("COEFFS[-1].shape:", COEFFS[-1].shape)

print('################################################')
calculate_memory_Psi(Nx, Ny, num_processes, precision_complex)
check_available_memory()
print('################################################')

print("Coefficients split completed.")

# Compute wavefunctions in real space
print("Computing wavefunctions in real space...")
phi = wavefunction_real_space(COEFFS, GVECS, grid_size)
print("Real-space wavefunctions computed.")

if __name__ == "__main__":

    if run_parallel:
        print("Running parallel computation...")
        RHO_ZE_ZH = calc_parallel(num_processes, RHO_ZE_ZH, verbose)
    else:
        print("Running serial computation...")
        RHO_ZE_ZH = calc_serial(RHO_ZE_ZH, verbose)


# Save the results in an HDF5 file
print("Saving results to RHO_ZE_ZH_data.h5...")
with h5py.File("RHO_ZE_ZH_data.h5", "w") as h5f:
    h5f.create_dataset("RHO_ZE_ZH", data=RHO_ZE_ZH)  # Main 3D array
    h5f.create_dataset("ZE_idx", data=ZE_idx)        # Save ZE_idx for reference
    h5f.create_dataset("ZH_idx", data=ZH_idx)        # Save ZH_idx for reference
print("Data saved to RHO_ZE_ZH_data.h5.")
print("RHO_ZE_ZH shape:", RHO_ZE_ZH.shape)

# End timer
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total execution time: {elapsed_time:.2f} seconds")
