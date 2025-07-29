
import time
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.fft import ifftn

# parameters
WFN_file = "/Users/rdelgrande/work/TEMP_DATA/teste_compute_rho_z_2d_bilayer/WFN.h5"
nc, nv = 2, 1

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

# load WFN.h5 file

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

# type(COEFFS) = list. Each element is a 3D array of shape (n_bands, 1, n_gvecs)

# create M_{n,n'}(k,q,G) matrix given by eq. 13 of https://arxiv.org/pdf/1111.4429
# M_{n,n'}(k,q,G) = FFT^-1 conj(psi_{n,k+q})(r) * psi_{n',k}(r)


def compute_M_nnprime_fixed_k_q(COEFFS, GVECS, k_idx, kq_idx, verbose=True):
    """
    Computes M_{nn'}(k, q, G) from plane-wave coefficients.
    
    Parameters:
    - COEFFS: list of arrays [n_kpoints] -> each array (n_bands, 1, n_gk)
    - GVECS: list of arrays [n_kpoints] -> each array (n_gk, 3)
    - k_idx: index of the k-point
    - kq_idx: index of the (k+q)-point

    Returns:
    - M: array (n_bands_kq, n_bands_k, n_gk_k) where n_gk_k is number of G-vectors at k
    """
    

    start_time = time.time()
    
    coeffs_k = COEFFS[k_idx][:,0,:]    # shape (n_bands_k, n_gk_k)
    coeffs_kq = COEFFS[kq_idx][:,0,:]  # shape (n_bands_kq, n_gk_kq)
    
    gvecs_k = GVECS[k_idx]   # shape (n_gk_k, 3)
    gvecs_kq = GVECS[kq_idx] # shape (n_gk_kq, 3)

    n_bands_kq, n_gkq = coeffs_kq.shape
    n_bands_k, n_gk = coeffs_k.shape

    # Build a dictionary for quick lookup: G_kq -> index
    gvecs_kq_dict = {tuple(g): i for i, g in enumerate(gvecs_kq)}

    M = np.zeros((n_bands_kq, n_bands_k, n_gk), dtype=complex)
    
    total_loops = n_bands_kq * n_bands_k
    loop_count = 0

    if verbose:
        print(f"Starting computation: {n_bands_kq} x {n_bands_k} = {total_loops} band pairs...")

    
    # M_{n,n'}(k,q,G) = V * sum_{G'} conj(C_{n,k+q}(G+G')) * C_{n',k}(G')
    # not including the V (volume) factor now
    
    if verbose:
        print(f"Creating lookup table for G+G'...")

    # Build a dictionary: G_kq -> index
    gvecs_kq_dict = {tuple(g): i for i, g in enumerate(gvecs_kq)}

    if verbose:
        print(f"Creating lookup table for G+G'...")

    # Precompute matching indices for each G
    lookup = []
    for ig, G in enumerate(gvecs_k):
        matched_igprime = []
        matched_ig_G_plus_Gprime = []
        for igprime, Gprime in enumerate(gvecs_kq):
            G_plus_Gprime = tuple(G + Gprime)
            if G_plus_Gprime in gvecs_kq_dict:
                matched_igprime.append(igprime)
                matched_ig_G_plus_Gprime.append(gvecs_kq_dict[G_plus_Gprime])
        lookup.append( (np.array(matched_igprime, dtype=int), np.array(matched_ig_G_plus_Gprime, dtype=int)) )

    if verbose:
        print(f"Lookup table created (elapsed {time.time()-start_time:.1f}s)")
        print(f"Starting main loops...")

    # Now perform the actual summations
    M = np.zeros((n_bands_kq, n_bands_k, n_gk), dtype=complex)

    for n in range(n_bands_kq):
        for n_prime in range(n_bands_k):
            for ig in range(n_gk):
                igprime_list, ig_GplusGprime_list = lookup[ig]
                if igprime_list.size > 0:
                    contrib = (np.conj(coeffs_kq[n, ig_GplusGprime_list]) *
                               coeffs_k[n_prime, igprime_list])
                    M[n, n_prime, ig] = np.sum(contrib)

            if verbose and ((n*n_bands_k + n_prime+1) % max(1, (n_bands_kq*n_bands_k)//10) == 0):
                print(f"Progress: {(n*n_bands_k+n_prime+1)}/{n_bands_kq*n_bands_k} band pairs processed (elapsed {time.time()-start_time:.1f}s)")

    # for n in range(n_bands_kq): # n
    #     for n_prime in range(n_bands_k): # n'
            
    #         loop_count += 1
    #         if verbose and loop_count % 10 == 0:
    #             elapsed = time.time() - start_time
    #             print(f"Progress: {loop_count}/{total_loops} band pairs processed (elapsed {elapsed:.1f}s)")
            
    #         for ig, G in enumerate(gvecs_k): # G
    #             for igprime, Gprime in enumerate(gvecs_kq): # G'
    #                 G_plus_Gprime = tuple(G + Gprime)  # Adjusted to use G and Gq
    #                 if G_plus_Gprime in gvecs_kq_dict:
    #                     ig_G_plus_Gprime = gvecs_kq_dict[G_plus_Gprime]
    #                     M[n, n_prime, ig] += np.conj(coeffs_kq[n, ig_G_plus_Gprime]) * coeffs_k[n_prime, igprime]

    end_time = time.time()
    if verbose:
        print(f"Finished in {end_time - start_time:.2f} seconds.")

    return M

    
def wrap_to_bz(kvec):
    """
    Wrap a k-vector back to the first Brillouin Zone.
    Assumes BZ is from -0.5 to 0.5 in each direction.
    """
    return (kvec + 0.5) % 1.0 - 0.5


def compute_M(COEFFS, GVECS, k_vectors, k_indices, q_vectors, verbose=True):
    """
    Compute M_{nn'}(k, q, G) for all given k's and q's.
    
    Parameters:
    - COEFFS: list of arrays for each k-point
    - GVECS: list of arrays for each k-point
    - k_vectors: array (n_kpoints, 3) of k-points
    - k_indices: list of indices of k-points in COEFFS/GVECS
    - q_vectors: array (n_qpoints, 3) of q-vectors

    Returns:
    - M_dict: dictionary with keys (k_idx, q_idx) and values = M matrices
    """

    M_dict = {}

    # Build a fast lookup: k-point tuple -> index
    kvec_dict = {tuple(np.round(kvec, 8)): idx for idx, kvec in zip(k_indices, k_vectors)}

    for k_count, (k_idx, kvec) in enumerate(zip(k_indices, k_vectors)):
        for q_count, qvec in enumerate(q_vectors):
            kq_vec = wrap_to_bz(kvec + qvec)
            kq_vec_rounded = tuple(np.round(kq_vec, 8))

            if kq_vec_rounded in kvec_dict:
                kq_idx = kvec_dict[kq_vec_rounded]
                if verbose:
                    print(f"Computing M for k_idx={k_idx} q_idx={q_count} (k+q matches k_idx={kq_idx})")

                M = compute_M_nnprime_fixed_k_q(COEFFS, GVECS, k_idx, kq_idx)
                M_dict[(k_idx, q_count)] = M

            else:
                if verbose:
                    print(f"Warning: (k+q) = {kq_vec} not found in k-points.")

    return M_dict
    


# truncated coulomb potential

# dielectric function

# Calculate kernel

# direct kernel is given by eq. 36 of https://arxiv.org/pdf/1111.4429
# <kcv|Kd|k'c'v'> = sum_{GG'} conj(M_{cc'}(k,q,G) W_{GG'}(q;0) M_{vv'}(k,q,G')
# in the diagonal approximation G=G'
# <kcv|Kd|k'c'v'> = sum_{G} conj(M_{cc'}(k,q,G) W_{GG}(q;0) M_{vv'}(k,q,G)
# W_{GG}(q;0) = v(q+G) / epsilon_{GG}(q)

# exchange kernel is given by eq. 37 of https://arxiv.org/pdf/1111.4429
# <kcv|Kx|k'c'v'> = sum_{G != 0} conj(M_{vc}(k,q,G) v(q+G) M_{v'c'}(k,q,G)
# which is already diagonal in G

