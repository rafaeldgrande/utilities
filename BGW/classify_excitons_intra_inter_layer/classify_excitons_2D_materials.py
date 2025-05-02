
import numpy as np
import h5py

eigenvectors_file = 'eigenvectors.h5'
x_kcv_file = 'x_kcv.h5'


# Open the eigenvectors file
with h5py.File(eigenvectors_file, "r") as f:
    A_vck = f['exciton_data/eigenvectors'][()]   # Adjust the key name based on actual structure
    eigenvalues = f['exciton_data/eigenvalues'][()]
A_vck = A_vck[0, :, :, :, :, 0, 0] + 1j * A_vck[0, :, :, :, :, 0, 1]
mod_A2 = np.real(np.abs(A_vck)**2)

Nexc, Nk, Nc, Nv = A_vck.shape
print("A_vck shape:", A_vck.shape)
print("Nk, Nc, Nv:", Nk, Nc, Nv)
print("Nexc:", Nexc)

# Open the x_kcv file
with h5py.File(x_kcv_file, "r") as f:
    x_kcv = f['x_kcv'][()]
print("x_kcv shape:", x_kcv.shape)

x_exc = np.zeros((Nexc))

# calculate x_exc, where x_exc[iexc] = x_kcv[ik, ic, iv] * mod_A2[iexc, ik, ic, iv]
x_exc = np.sum(mod_A2 * x_kcv[None, :, :, :], axis=(1, 2, 3))

print("x_exc shape:", x_exc.shape)

# save x_exc to a file
with h5py.File("x_exc.h5", "w") as f:
    f.create_dataset("x_exc", data=x_exc)
    f.create_dataset("eigenvalues", data=eigenvalues)
    
print('Finished writing x_exc.h5')