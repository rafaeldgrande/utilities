import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.fft import ifftn
import time
from multiprocessing import Pool
import argparse
import psutil

from mod_exciton_density_as_function_of_ze_zh import *
from itertools import product

# Parallelized IntPhi calculation functions
def compute_intphi_conduction_ze(args):
    """Worker function for parallel conduction band IntPhi calculation."""
    ze, phi_ze, Nc, Nk, volume_element = args
    
    # Initialize result array for this z-level
    IntPhi_c_ze = np.zeros((Nc, Nc, Nk, Nk), dtype=np.complex64)
    
    # Extract phi slices for all k-points and conduction bands at this z-level
    phi_c = phi_ze  # Already sliced: shape (nk, Nc, Nx, Ny)
    
    for ic1 in range(Nc):
        for ic2 in range(Nc):
            # Extract wavefunctions for bands ic1 and ic2
            psi1 = phi_c[:, ic1, :, :]  # shape: (nk, Nx, Ny)
            psi2 = phi_c[:, ic2, :, :]  # shape: (nk, Nx, Ny)
            
            # Compute integral using broadcasting: sum over spatial dimensions (Nx, Ny)
            integrals = np.sum(psi1[:, None, :, :] * np.conj(psi2[None, :, :, :]), axis=(2, 3))
            IntPhi_c_ze[ic1, ic2, :, :] = integrals
    
    # Apply volume element normalization
    IntPhi_c_ze = IntPhi_c_ze * volume_element
    
    return ze, IntPhi_c_ze

def compute_intphi_valence_zh(args):
    """Worker function for parallel valence band IntPhi calculation."""
    zh, phi_zh, Nv, Nk, volume_element = args
    
    # Initialize result array for this z-level
    IntPhi_v_zh = np.zeros((Nv, Nv, Nk, Nk), dtype=np.complex64)
    
    # Extract phi slices for all k-points and valence bands at this z-level
    phi_v = phi_zh  # Already sliced: shape (nk, Nv, Nx, Ny)
    
    for iv1 in range(Nv):
        for iv2 in range(Nv):
            # Map valence band indices (reverse order)
            band_v1 = Nv - 1 - iv1
            band_v2 = Nv - 1 - iv2
            
            # Extract wavefunctions for valence bands
            psi1 = phi_v[:, band_v1, :, :]  # shape: (nk, Nx, Ny)
            psi2 = phi_v[:, band_v2, :, :]  # shape: (nk, Nx, Ny)
            
            # Compute integral using broadcasting
            integrals = np.sum(psi1[:, None, :, :] * np.conj(psi2[None, :, :, :]), axis=(2, 3))
            IntPhi_v_zh[iv1, iv2, :, :] = integrals
    
    # Apply volume element normalization
    IntPhi_v_zh = IntPhi_v_zh * volume_element
    
    return zh, IntPhi_v_zh

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

parser.add_argument('--Ngrid_teste', type=int, default=4, help='Grid size for testing purposes')

parser.add_argument('--i_exc', type=int, default=0, help='Index of exciton to compute')
parser.add_argument('--limit_BSE_sum_up_to_value', type=float, default=1.0, help='Limit for BSE sum up to value')
parser.add_argument('--reduce_kpts_for_test', type=str2bool, default=False, help='Reduce k-points to 2 for testing purposes')

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
limitGvecs = args.limitGvecs
maxGs = args.maxGs
Ngrid_teste = args.Ngrid_teste
limit_BSE_sum_up_to_value = args.limit_BSE_sum_up_to_value
reduce_kpts_for_test = args.reduce_kpts_for_test

# Start the timer
start_time = time.time()


    
if run_parallel:
    print(f"Running in parallel with {num_processes} processes.")
else:
    print("Running in serial mode.")

if limitGvecs:
    print("Limiting G-vectors to", maxGs)
else:
    print("Not limiting G-vectors")

if reduce_kpts_for_test:
    print("WARNING: Reducing k-points to 2 for testing purposes!")
else:
    print("Using all k-points")
    

start_time_loading_data = time.time()
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

# print('TEST FOR top_n_indexes_all FUNCTION')
# for teste_bse_limit in np.arange(0.2, limit_BSE_sum_up_to_value, 0.1):
#     top_kcv = top_n_indexes_all(np.abs(A_vck[i_exc]), teste_bse_limit)
#     print(f"len(top_kcv) {len(top_kcv)} (instead of Nk*Nv*Nc = {Nk*Nc*Nv})  based teste_bse_limit = {teste_bse_limit}.") 

# Reading WFN file
print("Opening WFN file to read wavefunction data...")
with h5py.File(WFN_file, "r") as f:
    # coeffs_temp = f["/wfns/coeffs"][:]  # Check actual dataset name
    gvecs = f["/wfns/gvecs"][:]  # Reciprocal lattice vectors
    kpoints = f["/mf_header/kpoints/rk"][:]  # shape (3, nrk)
    print(f"Original kpoints shape from file: {kpoints.shape}")
    kpoints = kpoints.T  # shape (nrk, 3)
    print(f"After transpose, kpoints shape: {kpoints.shape}")
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

print("Reading wavefunction coefficients from WFN file")
with h5py.File(WFN_file, "r") as f:
    coeffs = f["/wfns/coeffs"][min_idx:max_idx, :, :, 0] + 1j * f["/wfns/coeffs"][min_idx:max_idx, :, :, 1]
    # shape of coeffs is (Nval+Nc-Nv, nk, ng)
Nbands = coeffs.shape[0]
print(f"coeffs shape: {coeffs.shape}")
elapsed_time_loading_data = time.time() - start_time_loading_data
print(f"Time spent loading data: {elapsed_time_loading_data:.2f} seconds")

# volume element for normalization
volume_element = 1 / (Nx * Ny * Nz)

print('Generating real-space grids')
# Generate real-space grids for electron and hole
x = np.linspace(0, 1, Nx, endpoint=False)
y = np.linspace(0, 1, Ny, endpoint=False)
z = np.linspace(0, 1, Nz, endpoint=False)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')  # shape: (Nx, Ny, Nz)
r_grid_xyz = np.stack((X, Y, Z), axis=-1)  # Shape: (Nx, Ny, Nz, 3)
print(f"Real-space grid shape: {r_grid_xyz.shape} (Nx={Nx}, Ny={Ny}, Nz={Nz})")

print("Reading wavefunction coefficients from WFN file")
with h5py.File(WFN_file, "r") as f:
    coeffs = f["/wfns/coeffs"][min_idx:max_idx, :, :, 0] + 1j * f["/wfns/coeffs"][min_idx:max_idx, :, :, 1]
    # shape of coeffs is (Nval+Nc-Nv, nk, ng)
Nbands = coeffs.shape[0]
print(f"coeffs shape: {coeffs.shape}")


# coeffs = coeffs_temp[min_idx:max_idx, :, :, 0] + 1j * coeffs_temp[min_idx:max_idx, :, :, 1]
print(f"Grid size: {grid_size}")

dz = Lz / Nz
z_real_value = np.arange(0, Lz, dz)

# Create lists of ZE_idx and ZH_idx
print("Generating ZE_idx and ZH_idx arrays...")

Nelectron = Nz
Nhole = Nz

ZE_idx = np.linspace(0, Nz-1, Nelectron, dtype=int)
ZH_idx = np.linspace(0, Nz-1, Nhole, dtype=int)
# print(f"ZE_idx: {ZE_idx}")
# print(f"ZH_idx: {ZH_idx}")

RHO_ZE_ZH = np.zeros((len(ZE_idx), len(ZH_idx)))  # Preallocate array

# Prepare arguments for parallel processing
# print("Preparing parameters for parallel computation...")
# params_list = [(i, j, ze_idx, zh_idx) for i, ze_idx in enumerate(ZE_idx) for j, zh_idx in enumerate(ZH_idx)]

params_list = []
for i, ze_idx in enumerate(ZE_idx):
    if not set_z_zero(z_real_value[ze_idx], zmin_set_zero, zmax_set_zero):
        continue
    for j, zh_idx in enumerate(ZH_idx):
        if not set_z_zero(z_real_value[zh_idx], zmin_set_zero, zmax_set_zero):
            continue
        params_list.append((i, j, ze_idx, zh_idx))
        
print(f"ZE_idx shape: {ZE_idx.shape}")
print(f"Total number of abs(Psi(ze, zh))**2 we will compute: {len(params_list)}")

print("Splitting WFN coefficients by k-points...")
COEFFS, GVECS = split_coeffs(coeffs, gvecs, ngk, maxGs, limitGvecs)

print(f"Original coeffs shape: {coeffs.shape}")
print(f"len(COEFFS) after split: {len(COEFFS)}")
print(f"COEFFS[0].shape: {COEFFS[0].shape}")
if len(COEFFS) > 1:
    print(f"COEFFS[1].shape: {COEFFS[1].shape}")
print(f"len(GVECS) after split: {len(GVECS)}")
print(f"kpoints.shape: {kpoints.shape}")
del coeffs, gvecs  # Free up memory

if len(COEFFS) != Nk:
    print(f"ERROR: COEFFS length {len(COEFFS)} does not match Nk {Nk}.")
    print(f"Check if WFN file ({WFN_file}) is consistent with eignvectors file ({eigenvectors_file}).")
    raise ValueError("Mismatch in number of k-points after splitting coefficients.")

# COEFFS shape is (nk, n_bands, 1, ng)
# print max value of ng in COEFFS
print(f"Highest number of G-vectors over k points: {max(ngk)}")
if limitGvecs:
    print("Limiting number of G-vectors for each k-point to be", maxGs)

if verbose > 0:
    print("len COEFFS:", len(COEFFS))
    print("COEFFS[0].shape:", COEFFS[0].shape)
    print("COEFFS[-1].shape:", COEFFS[-1].shape)
    
print("Coefficients split completed.")

if reduce_kpts_for_test:
    print("WARNING: Reducing k-points to 4 for testing purposes!")
    Nk = 5  # Reduce to Nk k-points for testing
    kpoints = kpoints[:, :Nk]  # Keep only first Nk k-points
    COEFFS = [COEFFS[ik] for ik in range(Nk)]
    GVECS = [GVECS[ik] for ik in range(Nk)]
    # phi = phi[:Nk, :, :, :, :]  # Keep only first Nk k-points
    A_vck = A_vck[:, :Nk, :, :]  # Keep only first Nk k-points
else:
    print("Using all k-points")

wavefunction_start_time = time.time()
# Compute wavefunctions in real space
print("Computing wavefunctions in real space...")
phi = wavefunction_real_space(COEFFS, GVECS, grid_size, kpoints, r_grid_xyz)
print("Real-space wavefunctions computed.")
print(f"phi shape: {phi.shape}")
print(f"Expected: (Nk={Nk}, Nbands={Nbands}, Nx={Nx}, Ny={Ny}, Nz={Nz})")
print(f"Nbands = Nv + Nc = {Nv} + {Nc} = {Nv + Nc}")
# shape phi is (nk, n_bands, Nx, Ny, Nz)
elapsed_time_wavefunction = time.time() - wavefunction_start_time
print(f"Time spent to compute wavefunctions in real space : {elapsed_time_wavefunction:.2f} seconds")

# getting important transitions in Akcv

top_kcv = top_n_indexes_all(np.abs(A_vck[i_exc]), limit_BSE_sum_up_to_value)
top_kcv_set = set(top_kcv)  # convert to set for O(1) lookup
# top_kcv[i] = ik, ic, iv to be used later

if len(top_kcv) < Nk * Nc * Nv:
    print(f"Limiting exciton coefficients to {len(top_kcv)} (instead of Nk*Nv*Nc = {Nk*Nc*Nv}) transitions based limit_BSE_sum_up_to_value = {limit_BSE_sum_up_to_value}.")

print('len(top_kcv):', len(top_kcv))
print('top_kcv[:2]:', top_kcv[:2])

# Create masks for valid z-coordinates
ze_valid = np.array([set_z_zero(z_real_value[ze_idx], zmin_set_zero, zmax_set_zero) 
                     for ze_idx in ZE_idx])
zh_valid = np.array([set_z_zero(z_real_value[zh_idx], zmin_set_zero, zmax_set_zero) 
                     for zh_idx in ZH_idx])

# Find valid indices
valid_ze_indices = np.where(ze_valid)[0]
valid_zh_indices = np.where(zh_valid)[0]

print(f"Total valid z-coordinates (ZE): {len(valid_ze_indices)}")
print(f"Total valid z-coordinates (ZH): {len(valid_zh_indices)}")

# # for testing purposes, limit the number of valid indices!
# valid_ze_indices = valid_ze_indices[:5]
# valid_zh_indices = valid_zh_indices[:5]


# Compute IntPhi(ic1, ic2, ik1, ik2, ze) = integral phi[ik1, ic1, :, :, ze] * conj(phi[ik2, ic2, :, :, ze]) dx dy

def compute_IntPhi(phi, ik1, ik2, band1, band2, ze):
    if ik1 >= phi.shape[0] or ik2 >= phi.shape[0]:
        print(f"ERROR: k-point index out of bounds! ik1={ik1}, ik2={ik2}, phi.shape[0]={phi.shape[0]}")
        raise IndexError(f"k-point index out of bounds")
    if band1 >= phi.shape[1] or band2 >= phi.shape[1]:
        print(f"ERROR: band index out of bounds! band1={band1}, band2={band2}, phi.shape[1]={phi.shape[1]}")
        raise IndexError(f"band index out of bounds")
    return np.sum(phi[ik1, band1, :, :, ze] * np.conj(phi[ik2, band2, :, :, ze]))

IntPhi_c = np.zeros((Nc, Nc, Nk, Nk, Nelectron), dtype=np.complex64)
IntPhi_v = np.zeros((Nv, Nv, Nk, Nk, Nhole), dtype=np.complex64)

print(f"Starting conduction band calculations...")
print(f"Nk from A_vck: {Nk}")
print(f"phi.shape[0] (actual k-points): {phi.shape[0]}")
print(f"Nc: {Nc}, Nv: {Nv}")
print(f"Adjusted IntPhi_c shape: {IntPhi_c.shape}")
print(f"Adjusted IntPhi_v shape: {IntPhi_v.shape}")

# Serial conduction band calculations (original code)
print("Computing conduction band integrals (serial)...")
conduction_start_time = time.time()

def conduction_worker(args):
    ze, ic1, ic2 = args
    # Only process valid ze
    if ze not in valid_ze_indices:
        return None
    # Extract wavefunctions for bands ic1 and ic2
    psi1 = phi[:, ic1, :, :, ze]  # shape: (nk, Nx, Ny)
    psi2 = phi[:, ic2, :, :, ze]  # shape: (nk, Nx, Ny)
    # Compute integral using broadcasting: sum over spatial dimensions (Nx, Ny)
    integrals = np.sum(psi1[:, None, :, :] * np.conj(psi2[None, :, :, :]), axis=(2, 3))
    return (ic1, ic2, ze, integrals)

def process_conduction_chunk(chunk):
    """Process a chunk of conduction band arguments and return results."""
    chunk_results = []
    for args in chunk:
        result = conduction_worker(args)
        if result is not None:
            chunk_results.append(result)
    return chunk_results

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# Prepare all combinations of (ze, ic1, ic2) for valid ze
conduction_args = list(product(valid_ze_indices, range(Nc), range(Nc)))

if run_parallel and num_processes > 1:
    print(f"Parallelizing conduction band calculation over {num_processes} processes...")
    
    # Split conduction_args into chunks for better progress tracking
    total_args = len(conduction_args)
    chunk_size = max(1, total_args // (num_processes * 4))  # Create more chunks than processes for better load balancing
    
    conduction_chunks = list(chunks(conduction_args, chunk_size))
    print(f"Split {total_args} argument sets into {len(conduction_chunks)} chunks of ~{chunk_size} sets each")
    
    # Process chunks in parallel with progress reporting
    processed_args = 0
    with Pool(processes=num_processes) as pool:
        # Submit all chunks to the pool
        chunk_results = pool.imap(process_conduction_chunk, conduction_chunks)
        
        # Process results as they complete
        for chunk_idx, chunk_result in enumerate(chunk_results):
            # Fill IntPhi_c with results
            for result in chunk_result:
                ic1, ic2, ze, integrals = result
                IntPhi_c[ic1, ic2, :, :, ze] = integrals
            
            # Update progress
            processed_args += len(chunk_result)
            progress_percent = 100 * processed_args / total_args
            
            print(f"Completed conduction chunk {chunk_idx + 1}/{len(conduction_chunks)}: "
                  f"{processed_args}/{total_args} argument sets ({progress_percent:.1f}%)")
            
            # Estimate remaining time based on current progress
            if chunk_idx > 0:  # After first chunk completion
                elapsed_time_so_far = time.time() - conduction_start_time
                time_per_arg = elapsed_time_so_far / processed_args
                remaining_args = total_args - processed_args
                estimated_remaining_time = time_per_arg * remaining_args
                print(f"  Estimated remaining time: {estimated_remaining_time:.1f} seconds")
else:
    print("Running conduction band calculation in serial mode...")
    for args in conduction_args:
        res = conduction_worker(args)
        if res is not None:
            ic1, ic2, ze, integrals = res
            IntPhi_c[ic1, ic2, :, :, ze] = integrals
        # Progress report for conduction bands
        if (ze + 1) % 10 == 0 or ze == Nelectron - 1:
            print(f"  Conduction bands: processed {ze + 1}/{Nelectron} z-levels ({100*(ze+1)/Nelectron:.1f}%)")

# Apply volume element normalization
IntPhi_c = IntPhi_c * volume_element

conduction_end_time = time.time()
conduction_elapsed = conduction_end_time - conduction_start_time
print(f"Conduction band calculations completed in {conduction_elapsed:.2f} seconds \n\n")

valence_start_time = time.time()

# Prepare all combinations of (zh, iv1, iv2) for valid zh
valence_args = list(product(valid_zh_indices, range(Nv), range(Nv)))

def valence_worker(args):
    zh, iv1, iv2 = args
    # Only process valid zh
    if zh not in valid_zh_indices:
        return None
    # Map valence band indices (reverse order)
    band_v1 = Nv - 1 - iv1
    band_v2 = Nv - 1 - iv2
    # Extract wavefunctions for valence bands
    psi1 = phi[:, band_v1, :, :, zh]  # shape: (nk, Nx, Ny)
    psi2 = phi[:, band_v2, :, :, zh]  # shape: (nk, Nx, Ny)
    # Compute integral using broadcasting
    integrals = np.sum(psi1[:, None, :, :] * np.conj(psi2[None, :, :, :]), axis=(2, 3))
    return (iv1, iv2, zh, integrals)

def process_valence_chunk(chunk):
    """Process a chunk of valence band arguments and return results."""
    chunk_results = []
    for args in chunk:
        result = valence_worker(args)
        if result is not None:
            chunk_results.append(result)
    return chunk_results

if run_parallel and num_processes > 1:
    print(f"Parallelizing valence band calculation over {num_processes} processes...")
    
    # Split valence_args into chunks for better progress tracking
    total_args = len(valence_args)
    chunk_size = max(1, total_args // (num_processes * 4))  # Create more chunks than processes for better load balancing
    
    valence_chunks = list(chunks(valence_args, chunk_size))
    print(f"Split {total_args} argument sets into {len(valence_chunks)} chunks of ~{chunk_size} sets each")
    
    # Process chunks in parallel with progress reporting
    processed_args = 0
    with Pool(processes=num_processes) as pool:
        # Submit all chunks to the pool
        chunk_results = pool.imap(process_valence_chunk, valence_chunks)
        
        # Process results as they complete
        for chunk_idx, chunk_result in enumerate(chunk_results):
            # Fill IntPhi_v with results
            for result in chunk_result:
                iv1, iv2, zh, integrals = result
                IntPhi_v[iv1, iv2, :, :, zh] = integrals
            
            # Update progress
            processed_args += len(chunk_result)
            progress_percent = 100 * processed_args / total_args
            
            print(f"Completed valence chunk {chunk_idx + 1}/{len(valence_chunks)}: "
                  f"{processed_args}/{total_args} argument sets ({progress_percent:.1f}%)")
            
            # Estimate remaining time based on current progress
            if chunk_idx > 0:  # After first chunk completion
                elapsed_time_so_far = time.time() - valence_start_time
                time_per_arg = elapsed_time_so_far / processed_args
                remaining_args = total_args - processed_args
                estimated_remaining_time = time_per_arg * remaining_args
                print(f"  Estimated remaining time: {estimated_remaining_time:.1f} seconds")
else:
    print("Running valence band calculation in serial mode...")
    for args in valence_args:
        res = valence_worker(args)
        if res is not None:
            iv1, iv2, zh, integrals = res
            IntPhi_v[iv1, iv2, :, :, zh] = integrals
        # Progress report for valence bands
        zh = args[0]
        if (zh + 1) % 20 == 0 or zh == Nhole - 1:
            print(f"  Valence bands: processed {zh + 1}/{Nhole} z-levels ({100*(zh+1)/Nhole:.1f}%)")

# Apply volume element normalization
IntPhi_v = IntPhi_v * volume_element

valence_end_time = time.time()
valence_elapsed = valence_end_time - valence_start_time
print(f"Valence band calculations completed in {valence_elapsed:.2f} seconds")

total_intphi_time = conduction_elapsed + valence_elapsed
print(f"Total IntPhi calculation time: {total_intphi_time:.2f} seconds")
print('Finished computing IntPhi_c and IntPhi_v.')

# Optimized density calculation using vectorized operations
print("Calculating density RHO_ZE_ZH (vectorized)...")
start_time_comp_density = time.time()
RHO_ZE_ZH = np.zeros((len(ZE_idx), len(ZH_idx)))  # Preallocate array


def compute_rho_ze_zh(ze, zh, IntPhi_c, IntPhi_v, A_vck, i_exc, top_kcv):
    """
    Optimized vectorized computation of exciton density at specific ze, zh coordinates.
    
    Original nested loops: O(N²) where N = len(top_kcv)
    Optimized version: O(N) using NumPy vectorization
    """
    if len(top_kcv) == 0:
        return 0.0
    
    # Convert top_kcv to arrays for vectorized operations
    top_kcv_array = np.array(top_kcv)  # Shape: (N, 3) where N = len(top_kcv)
    ik_vals = top_kcv_array[:, 0]      # All ik values
    ic_vals = top_kcv_array[:, 1]      # All ic values  
    iv_vals = top_kcv_array[:, 2]      # All iv values
    
    # Extract A_vck coefficients for all transitions at once
    # A_coeffs shape: (N,) - coefficients for all (ik, ic, iv) combinations
    A_coeffs = A_vck[i_exc, ik_vals, ic_vals, iv_vals]
    
    # Extract IntPhi values using advanced indexing
    # IntPhi_c_vals shape: (N, N) - all pairwise combinations
    IntPhi_c_vals = IntPhi_c[ic_vals[:, None], ic_vals[None, :], 
                             ik_vals[:, None], ik_vals[None, :], ze]
    
    # IntPhi_v_vals shape: (N, N) - all pairwise combinations  
    # Note: indices are swapped (ik2, ik1) for IntPhi_v
    IntPhi_v_vals = IntPhi_v[iv_vals[:, None], iv_vals[None, :],
                             ik_vals[None, :], ik_vals[:, None], zh]
    
    # Compute all pairwise products using broadcasting
    # A_coeffs[:, None] * np.conj(A_coeffs[None, :]) creates (N, N) matrix
    A_product = A_coeffs[:, None] * np.conj(A_coeffs[None, :])
    
    # Final vectorized computation - element-wise multiplication and sum
    rho = np.sum(IntPhi_c_vals * IntPhi_v_vals * A_product)
    
    return np.real(rho)

print(f"Computing density for {len(valid_ze_indices)} x {len(valid_zh_indices)} = {len(valid_ze_indices) * len(valid_zh_indices)} valid z-coordinates...")

# Optimized vectorized computation for all valid combinations

# Precompute all (ze, zh) pairs to process
ze_zh_pairs = [(ze, zh) for ze in valid_ze_indices for zh in valid_zh_indices]

def rho_worker(args):
    ze, zh = args
    return ze, zh, compute_rho_ze_zh(ze, zh, IntPhi_c, IntPhi_v, A_vck, i_exc, top_kcv)

if run_parallel and num_processes > 1:
    print(f"Parallelizing density calculation over {num_processes} processes...")
    
    # Split ze_zh_pairs into chunks for better progress tracking
    total_pairs = len(ze_zh_pairs)
    chunk_size = max(1, total_pairs // (num_processes * 4))  # Create more chunks than processes for better load balancing
    
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]
    
    ze_zh_chunks = list(chunks(ze_zh_pairs, chunk_size))
    print(f"Split {total_pairs} pairs into {len(ze_zh_chunks)} chunks of ~{chunk_size} pairs each")
    
    def process_chunk(chunk):
        """Process a chunk of ze_zh pairs and return results."""
        chunk_results = []
        for ze, zh in chunk:
            rho = compute_rho_ze_zh(ze, zh, IntPhi_c, IntPhi_v, A_vck, i_exc, top_kcv)
            chunk_results.append((ze, zh, rho))
        return chunk_results
    
    # Process chunks in parallel with progress reporting
    processed_pairs = 0
    with Pool(processes=num_processes) as pool:
        # Submit all chunks to the pool
        chunk_results = pool.imap(process_chunk, ze_zh_chunks)
        
        # Process results as they complete
        for chunk_idx, chunk_result in enumerate(chunk_results):
            # Store results in the main array
            for ze, zh, rho in chunk_result:
                RHO_ZE_ZH[ze, zh] = rho
            
            # Update progress
            processed_pairs += len(chunk_result)
            progress_percent = 100 * processed_pairs / total_pairs
            
            print(f"Completed chunk {chunk_idx + 1}/{len(ze_zh_chunks)}: "
                  f"{processed_pairs}/{total_pairs} pairs ({progress_percent:.1f}%)")
            
            # Estimate remaining time based on current progress
            if chunk_idx > 0:  # After first chunk completion
                elapsed_time_so_far = time.time() - start_time_comp_density
                time_per_pair = elapsed_time_so_far / processed_pairs
                remaining_pairs = total_pairs - processed_pairs
                estimated_remaining_time = time_per_pair * remaining_pairs
                print(f"  Estimated remaining time: {estimated_remaining_time:.1f} seconds")
                
else:
    print("Running density calculation in serial mode...")
    counter = 0
    total_pairs = len(ze_zh_pairs)
    for ze, zh in ze_zh_pairs:
        # report progress every 100 iterations
        if counter % 100 == 0 or counter == 5:
            print(f"Processed {counter}/{total_pairs} pairs ({100 * counter / total_pairs:.1f}%)")
        RHO_ZE_ZH[ze, zh] = compute_rho_ze_zh(ze, zh, IntPhi_c, IntPhi_v, A_vck, i_exc, top_kcv)
        counter += 1

elapsed_time_comp_density = time.time() - start_time_comp_density
print(f"Density RHO_ZE_ZH computed in {elapsed_time_comp_density:.2f} seconds")

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

print(f"Total time to load data: {elapsed_time_loading_data:.2f} seconds")
print(f"Total time to compute wavefunctions in real space : {elapsed_time_wavefunction:.2f} seconds")
print(f"Total IntPhi calculation time: {conduction_elapsed + valence_elapsed:.2f} seconds")
print(f"Total time to compute density RHO_ZE_ZH: {elapsed_time_comp_density:.2f} seconds")