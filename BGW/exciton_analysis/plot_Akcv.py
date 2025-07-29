
import numpy as np
import matplotlib.pyplot as plt
import h5py

eigenvectors_file = "eigenvectors.h5"  
iexc = 0

# Open the HDF5 file
print("Opening eigenvectors file to read exciton coefficients...")
with h5py.File(eigenvectors_file, "r") as f:
    A_vck = f['exciton_data/eigenvectors'][()]   # Adjust the key name based on actual structure
A_vck = A_vck[0, :, :, :, :, 0, 0] + 1j * A_vck[0, :, :, :, :, 0, 1] # shape (Nexc, Nk, Nc, Nv)
Nexc, Nk, Nc, Nv = A_vck.shape
A_iexc = A_vck[iexc, :, :, :]  # shape (Nk, Nc, Nv)
print("Exciton coefficients loaded.")
print("A_vck shape:", A_vck.shape)
print("Nk, Nc, Nv:", Nk, Nc, Nv)
print("Nexc:", Nexc)

A_iexc_flat_vc_k = A_iexc.reshape(Nk, Nc * Nv)  # shape (Nk, Nc*Nv)
A_iexc_v_to_c_summed_k = np.sum(np.abs(A_iexc)**2, axis=0)  # shape (nc, nv)

# make a heatmap of the exciton coefficients
plt.figure(figsize=(10, 6))
plt.imshow(np.log(np.abs(A_iexc_flat_vc_k)), aspect='auto', cmap='viridis')
plt.colorbar(label='log(Akcv)')
plt.xlabel('transition v -> c index')
plt.ylabel('k-point index')
plt.savefig('exciton_coefficients.png', dpi=300)

# heatmap for A_iexc_v_to_c_summed_k
plt.figure(figsize=(10, 6))
plt.imshow(A_iexc_v_to_c_summed_k, aspect='auto', cmap='viridis')
plt.colorbar(label='sum_k |Akcv|**2')
plt.xlabel('v band index')
plt.ylabel('c band index')
plt.savefig('exciton_coefficients_summed_over_k.png', dpi=300)