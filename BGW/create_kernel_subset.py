import h5py
import numpy as np
import argparse

# author: RRDG
# date: 2025/2/27

''' This script modifies the number of valence and conduction bands in a BSE kernel file.
    It reads  the original file bsemat.h5 and creates a new file bsemat_mod.h5 with the 
    modified number of bands. Assuming restricted kernel
    
    Usage:
    python create_kernel_subset.py nc_new nv_new
    where nc_new and nv_new are the new number of conduction and valence bands, respectively.
'''

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Modify bsemat.h5 by reducing conduction and valence bands.")
parser.add_argument("mv", type=int, help="New number of valence bands")
parser.add_argument("mc", type=int, help="New number of conduction bands")
args = parser.parse_args()

# Open the original file
with h5py.File("bsemat.h5", "r") as f_orig:
    # Read original number of valence and conduction bands
    nv = f_orig["/bse_header/bands/nvb"][()]
    nc = f_orig["/bse_header/bands/ncb"][()]

    print(f"Original kernel: nv = {nv}, nc = {nc}")
    print(f"New kernel: nv = {args.mv}, nc = {args.mc}")

    # Create a new file
    with h5py.File("bsemat_mod.h5", "w") as f_new:
        # Copy all groups and datasets except the ones we need to modify
        for name in f_orig:
            f_orig.copy(name, f_new)

        # Modify specific datasets
        f_new["/bse_header/bands/nvb"][()] = args.mv
        f_new["/bse_header/bands/ncb"][()] = args.mc
        f_new["/bse_header/bands/n1b"][()] = args.mv  # Assuming restricted kernel
        f_new["/bse_header/bands/n2b"][()] = args.mc  # Assuming restricted kernel

        # Modify the kernel matrix datasets
        for mat_name in ["head", "wing", "body", "exchange", "fxc"]:
            mat_path = f"/mats/{mat_name}"
            if mat_path in f_new:
                data = f_new[mat_path][:]

                # Determine new shape
                new_shape = list(data.shape)
                new_shape[1] = args.mv  # Modify valence dimension
                new_shape[3] = args.mc  # Modify conduction dimension

                # Create a new dataset with the modified shape
                del f_new[mat_path]  # Delete the existing dataset
                f_new.create_dataset(mat_path, shape=tuple(new_shape), dtype=data.dtype)

                # Copy relevant data
                f_new[mat_path][:] = data[: args.mv, : args.mv, : args.mc, : args.mc, :, :]

print("New bsemat_mod.h5 file created successfully.")
