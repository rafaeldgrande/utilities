
import h5py
import shutil
import os

# Original files
source_header_file = "../../DFT/6-wfn_fi/WFN.h5"
base_file = "../../DFPT/wfn_fi_dfpt/WFN.h5"
output_file = "WFN_fi_mod.h5"

# Step 1: Copy the entire original file to the new one
shutil.copy(base_file, output_file)
print(f"Copied {base_file} to {output_file}")

# Step 2: Replace /mf_header in the new file
with h5py.File(source_header_file, "r") as src, h5py.File(output_file, "a") as dst:
    # Remove old /mf_header if it exist
    if "/mf_header" in dst:
        del dst["/mf_header"]
        print("Deleted existing /mf_header in output file.")

    # Copy new /mf_header from source
    src.copy("/mf_header", dst)
    print("Copied new /mf_header from", source_header_file)
