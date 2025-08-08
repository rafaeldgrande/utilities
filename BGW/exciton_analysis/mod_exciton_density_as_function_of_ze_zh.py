
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.fft import ifftn
import time
from multiprocessing import Pool
import argparse
import psutil

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def set_z_zero(z, zmin_set_zero, zmax_set_zero):
    if zmin_set_zero < z < zmax_set_zero:
        return False
    else:
        return True
    


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
        if np.abs(limit_BSE_sum_up_to_value - sum_abs_Akcv2) < 1e-6:
            break
        elif sum_abs_Akcv2 > limit_BSE_sum_up_to_value:
            # If the cumulative sum exceeds the limit, we stop
            break

    return top_indexes[:counter_indexes]

def split_coeffs(coeffs, gvecs, ngk, maxGs, limitGvecs):
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

def wavefunction_for_k_point(coeffs, gvecs, grid_size, kvec, r_grid_xyz, precision_complex=np.complex64):
    """Compute real-space wavefunction via inverse FFT."""
    n_bands, _, _ = coeffs.shape
    psi_k_point = np.zeros((n_bands, *grid_size), dtype=precision_complex)
    
    # print('!!!!!!!!!!!!!')
    # print("r_grid_xyz shape:", r_grid_xyz.shape)
    # print("kvec shape:", kvec.shape)
    
    # Compute Bloch phase: exp(i k⋅r)
    phase = np.exp(1j * np.tensordot(r_grid_xyz, kvec, axes=([3], [0])))  # shape: (Nx, Ny, Nz)

    for n in range(n_bands):
        # Place coefficients in a 3D reciprocal-space grid
        grid = np.zeros(grid_size, dtype=precision_complex)
        for i, G in enumerate(gvecs):
            grid[tuple(G)] = coeffs[n, 0, i]

        # Inverse FFT to get real-space wavefunction
        u_nk = ifftn(grid) # periodic part u_nk(r). Shape is (Nx, Ny, Nz)
        
        # Multiply by Bloch phase to get full wavefunction
        psi_nk = u_nk * phase
        # norm = np.sqrt(np.sum(np.abs(psi_nk)**2))
        # print("norm:", norm)
        psi_k_point[n] = psi_nk #/ norm
        
    return psi_k_point

def wavefunction_real_space(COEFFS, GVECS, grid_size, kvecs, r_grid_xyz, precision_complex=np.complex64):
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
    Nk = kvecs.shape[1]
    n_bands, _, _ = COEFFS[0].shape
    psi_real_space = np.zeros((Nk, n_bands, *grid_size), dtype=precision_complex)
    
    for ik in range(Nk):
        psi_real_space[ik] = wavefunction_for_k_point(COEFFS[ik], GVECS[ik], grid_size, kvecs[:, ik], r_grid_xyz, precision_complex)
        if ik % 10 == 0:
            print(f"Processed k-point {ik+1}/{Nk}")

    return psi_real_space