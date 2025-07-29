import h5py
import numpy as np
import argparse

# author: RRDG
# date: 2025/2/27

'''
This script modifies the number of valence and conduction bands in a BSE kernel file.
It reads the original file (default: bsemat.h5) and creates a new file (default: bsemat_mod.h5)
with the modified number of bands. Assuming restricted kernel.

Usage:
    python create_kernel_subset.py --nc <new_conduction_bands> --nv <new_valence_bands> \
                                   --input bsemat.h5 --output bsemat_mod.h5
'''

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Modify bsemat.h5 by reducing conduction and valence bands.")
parser.add_argument("--nv", type=int, required=True, help="New number of valence bands")
parser.add_argument("--nc", type=int, required=True, help="New number of conduction bands")
parser.add_argument("--input", type=str, default="bsemat.h5", help="Input kernel file (default: bsemat.h5)")
parser.add_argument("--output", type=str, default="bsemat_mod.h5", help="Output kernel file (default: bsemat_mod.h5)")
args = parser.parse_args()

nv_new = args.nv
nc_new = args.nc
input_file = args.input
output_file = args.output

# Open the original file
with h5py.File("bsemat.h5", "r") as f_orig:
    # Read original number of valence and conduction bands
    nv_orig = f_orig["/bse_header/bands/nvb"][()]
    nc_orig = f_orig["/bse_header/bands/ncb"][()]
    
    print(f"Original kernel: {input_file} with nv = {nv_orig}, nc = {nc_orig}")
    print(f"New kernel: {output_file} with nv = {nv_new}, nc = {nc_new}")
    
    if nc_orig < nc_new:
        raise ValueError(f"New number of conduction bands ({nc_new}) is greater than the original number ({nc_orig}).")
    if nv_orig < nv_new:
        raise ValueError(f"New number of valence bands ({nv_new}) is greater than the original number ({nv_orig}).")



    # Create a new file
    with h5py.File(output_file, "w") as f_new:
        # Copy all groups and datasets except the ones we need to modify
        for name in f_orig:
            f_orig.copy(name, f_new)

        # Modify specific datasets
        f_new["/bse_header/bands/nvb"][()] = nv_new
        f_new["/bse_header/bands/ncb"][()] = nc_new
        f_new["/bse_header/bands/n1b"][()] = nv_new  # Assuming restricted kernel
        f_new["/bse_header/bands/n2b"][()] = nc_new  # Assuming restricted kernel

        # Modify the kernel matrix datasets
        for mat_name in ["head", "wing", "body", "exchange", "fxc"]:
            mat_path = f"/mats/{mat_name}"
            if mat_path in f_new:
                data = f_new[mat_path][:]
                
                # shape of data
                # nk*ns, nk*ns, nc_orig, nc_orig, nv_orig, nv_orig, flavor
                # nk = number of k-points, ns = spin components
                # flavo = 1 (real) or 2 (complex)

                # Determine new shape
                new_shape = list(data.shape)
                new_shape[2] = nc_new  # Modify valence dimension
                new_shape[3] = nc_new  # Modify valence dimension
                new_shape[4] = nv_new  # Modify conduction dimension
                new_shape[5] = nv_new  # Modify conduction dimension
                
                if mat_name == 'head':
                    print(f'Previous shape {data.shape}')
                    print(f"New shape {new_shape}")

                # Create a new dataset with the modified shape
                del f_new[mat_path]  # Delete the existing dataset
                f_new.create_dataset(mat_path, shape=tuple(new_shape), dtype=data.dtype)

                # Copy relevant data
                f_new[mat_path][:] = data[:, :, :nc_new, :nc_new, :nv_new, :nv_new, :]
                
                print('Finished copying', mat_name)

print(f"New {output_file} file created successfully.")
