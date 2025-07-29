
import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file_name', type=str, default='RHO_ZE_ZH_data.h5', help='HDF5 file produced by exciton_density_as_function_of_ze_zh.py')
parser.add_argument('--fig_name', type=str, default='RHO_ZE_ZH_heatmap.png', help='Figure name for the heatmap output')

args = parser.parse_args()
filename = args.file_name
fig_name = args.fig_name

Lz = 30

z1 = -0.850290 + Lz
z2 = 5.650293 + Lz
z_div = (z1 + z2) / 2 

# Load the data
with h5py.File(filename, "r") as h5f:
    RHO_ZE_ZH = h5f["RHO_ZE_ZH"][:]
    ZE_idx = h5f["ZE_idx"][:]
    ZH_idx = h5f["ZH_idx"][:]
    
Nz = len(ZE_idx)
dz = Lz / Nz
    
ZE_idx_new = [i for i in range(2*len(ZE_idx))]
ZH_idx_new = [i for i in range(2*len(ZH_idx))]

ZE_idx_new = np.array(ZE_idx_new) * dz - z_div
ZH_idx_new = np.array(ZH_idx_new) * dz - z_div

RHO_ZE_ZH_new = np.zeros((len(ZE_idx_new), len(ZH_idx_new)))

for i in range(len(ZE_idx)):
    for j in range(len(ZH_idx)):
        RHO_ZE_ZH_new[i, j] = RHO_ZE_ZH[i, j]
        RHO_ZE_ZH_new[i + len(ZE_idx), j] = RHO_ZE_ZH[i, j]
        RHO_ZE_ZH_new[i, j + len(ZH_idx)] = RHO_ZE_ZH[i, j]
        RHO_ZE_ZH_new[i + len(ZE_idx), j + len(ZH_idx)] = RHO_ZE_ZH[i, j]


# Plot the heatmap
plt.figure(figsize=(8, 6))
plt.imshow(RHO_ZE_ZH_new.T, origin='lower', aspect='auto', cmap='plasma',
           extent=[ZE_idx_new[0], ZE_idx_new[-1], ZH_idx_new[0], ZH_idx_new[-1]],
           interpolation='bicubic')
plt.colorbar(label=r"$\rho(z_e, z_h)$")
plt.xlabel(r"$z_e \ (\rm{\AA})$")
plt.ylabel(r"$z_h \ (\rm{\AA})$")

# add contour lines
ze_grid, zh_grid = np.meshgrid(ZE_idx_new, ZH_idx_new)  # Match shape of original data
CS = plt.contour(ze_grid, zh_grid, RHO_ZE_ZH_new, levels=10, colors='white', linewidths=0.8)
# plt.clabel(CS, inline=True, fontsize=8, fmt="%.2f")  # Optional: add labels to the lines

# plot vertical and horizontal lines at z1 and z2
plt.axvline(x=z1 - z_div, color='k', linestyle='--')
plt.axvline(x=z2 - z_div, color='k', linestyle='--')
plt.axhline(y=z1 - z_div, color='k', linestyle='--')
plt.axhline(y=z2 - z_div, color='k', linestyle='--')

# plot diagonal line x = y
plt.plot([-10, 10], [-10, 10], 'k--')

plt.xlim(-7.5, 7.5)
plt.ylim(-7.5, 7.5)

# Save the figure
plt.savefig(fig_name, dpi=300, bbox_inches='tight')


def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()

ze1_idx = find_nearest_index(ZE_idx_new, z1 - z_div)
ze2_idx = find_nearest_index(ZE_idx_new, z2 - z_div)
zh1_idx = find_nearest_index(ZH_idx_new, z1 - z_div)
zh2_idx = find_nearest_index(ZH_idx_new, z2 - z_div)

# Extract the cuts
rho_ze1 = RHO_ZE_ZH_new[ze1_idx, :]
rho_ze2 = RHO_ZE_ZH_new[ze2_idx, :]
rho_zh1 = RHO_ZE_ZH_new[:, zh1_idx]
rho_zh2 = RHO_ZE_ZH_new[:, zh2_idx]

# Diagonal cut: z_e = z_h
z_diag = np.array([val for val in ZE_idx_new if val in ZH_idx_new])
idx_diag = [find_nearest_index(ZE_idx_new, z) for z in z_diag]
rho_diag = np.array([RHO_ZE_ZH_new[i, i] for i in idx_diag])


# Plot the cuts
plt.figure(figsize=(10, 6))

plt.plot(ZH_idx_new, rho_ze1, 'r-', label=fr"$\rho(z_e = z_1, z_h)$")
plt.plot(ZH_idx_new, rho_ze2, 'r--', label=fr"$\rho(z_e = z_2, z_h)$")
plt.plot(ZE_idx_new, rho_zh1, 'b-', label=fr"$\rho(z_e, z_h = z_1)$")
plt.plot(ZE_idx_new, rho_zh2, 'b--', label=fr"$\rho(z_e, z_h = z_2)$")
plt.plot(z_diag, rho_diag, '-.', color='black', label=r"$\rho(z_e = z_h)$")

# vertical lineas at z1 and z2
plt.axvline(x=z1 - z_div, color='k', linestyle='-', linewidth=0.5)
plt.axvline(x=z2 - z_div, color='k', linestyle='-', linewidth=0.5)


# Add text annotations near the top
y_max = max(max(rho_ze1), max(rho_ze2), max(rho_zh1), max(rho_zh2), max(rho_diag))
plt.text(z1 - z_div + 0.1, y_max * 0.98, r"$z_1$", rotation=0, verticalalignment='top')
plt.text(z2 - z_div + 0.1, y_max * 0.98, r"$z_2$", rotation=0, verticalalignment='top')

plt.ylim([0, y_max * 1.03])

plt.xlabel(r"$z \ (\rm{\AA})$")
plt.ylabel(r"$\rho(z_e, z_h)$")
plt.title("Exciton Density Cuts")
plt.legend()
# plt.grid(True)
plt.xlim(-7.5, 7.5)

# Save the cuts figure
plt.savefig("RHO_ZE_ZH_cuts.png", dpi=300, bbox_inches='tight')



