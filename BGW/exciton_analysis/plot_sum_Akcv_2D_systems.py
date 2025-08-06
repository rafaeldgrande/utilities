
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

'''
Plot exciton eigenvector distributions in k space and 
for transitions v to c.

Usage:
python plot_sum_cv_Akcv_2D_systems.py --filename eigenvectors.h5 --i_exc_max 1'''


# === Argument parser ===
parser = argparse.ArgumentParser(description="Plot exciton eigenvector distributions.")
parser.add_argument("--filename", type=str, default="eigenvectors.h5", help="Path to eigenvectors.h5 file")
parser.add_argument("--i_exc_max", type=int, default=1, help="Number of excitons to plot (starting from 0)")
parser.add_argument("--finite_momentum_exciton", type=bool, default=False, help="Is this a finite momentum exciton?")
parser.add_argument("--save_pdf", type=bool, default=True, help="Save plots to PDF")

args = parser.parse_args()

filename = args.filename
i_exc_max = args.i_exc_max
fin_mom_exciton = args.finite_momentum_exciton
save_pdf = args.save_pdf

# === Load data once ===
with h5py.File(filename, "r") as f:
    flavor = f["/exciton_header/flavor"][()]
    assert flavor == 2, "Expected complex eigenvectors (flavor=2)"
    ns = f["/exciton_header/params/ns"][()]
    nc = f["/exciton_header/params/nc"][()]
    nv = f["/exciton_header/params/nv"][()]
    nk = f["/exciton_header/kpoints/nk"][()]
    nevecs = f["/exciton_header/params/nevecs"][()]
    
    Akcv_all = f["/exciton_data/eigenvectors"][:]
    Akcv = Akcv_all[0, :, :, :, :, 0, 0] + 1j * Akcv_all[0, :, :, :, :, 0, 1]

    kpts = f["/exciton_header/kpoints/kpts"][:]  # shape: (3, nk)
    bvec = f["/mf_header/crystal/bvec"][:]       # shape: (3, 3)
    
    if fin_mom_exciton:
        try:
            Qshifts_temp = f["/exciton_header/kpoints/exciton_Q_shifts"][:] #[:, 0]
            Qshift = Qshifts_temp[0, 0]*bvec[0] + Qshifts_temp[0, 1]*bvec[1] + Qshifts_temp[0, 2]*bvec[2]
            print('Qshift found')
            print(f'Qshift = {Qshift}, shape = {Qshift.shape}')
        except:
            print('Qshift not found. Setting fin_mom_exciton to False')
            fin_mom_exciton = False
            
# === Convert k-points to Cartesian coordinates ===
k_cart = [kpts[ik, 0] * bvec[0] + kpts[ik, 1] * bvec[1] + kpts[ik, 2] * bvec[2] for ik in range(len(kpts))]
k_cart = np.array(k_cart)
kx_orig = k_cart[:, 0]
ky_orig = k_cart[:, 1]

# Reciprocal vectors in 2D
b1 = bvec[0, :2]
b2 = bvec[1, :2]

# First BZ hexagon corners
kpts_hex = np.array([
    (1/3) * ( b1 + 2 * b2),
    (1/3) * (2 * b1 + b2),
    (1/3) * ( b1 - b2),
    (1/3) * (-b1 - 2 * b2),
    (1/3) * (-2 * b1 - b2),
    (1/3) * (-b1 + b2)
])
angles = np.arctan2(kpts_hex[:,1], kpts_hex[:,0])
kpts_hex_sorted = kpts_hex[np.argsort(angles)]
kpts_hex_closed = np.vstack([kpts_hex_sorted, kpts_hex_sorted[0]])

# === Plot all excitons and save to a PDF ===
if save_pdf:
    pdf = PdfPages("exciton_Akcv_heatmaps.pdf")

for i_exc in range(i_exc_max):
    Akcv_exc = Akcv[i_exc]  # shape: (nk, nc, nv)
    sum_cv_Akcv = np.sum(np.abs(Akcv_exc)**2, axis=(1, 2))
    sum_k_Akcv = np.sum(np.abs(Akcv_exc)**2, axis=(0))

    # Repeat for shifted k-points
    kx = np.copy(kx_orig)
    ky = np.copy(ky_orig)
    values = np.copy(sum_cv_Akcv)

    for v in [b1, b2, -b1, -b2, b1+b2, -b1-b2]:
        delta_kx, delta_ky = v
        kx = np.append(kx, kx_orig + delta_kx)
        ky = np.append(ky, ky_orig - delta_ky)
        values = np.append(values, sum_cv_Akcv)

    # Plot
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.tripcolor(kx, ky, values, shading='flat', cmap='plasma')
    plt.colorbar(sc, ax=ax, label=r'$\sum_{cv} |A_{kcv}|^2$')

    ax.plot(kpts_hex_closed[:, 0], kpts_hex_closed[:, 1], 'k-', linewidth=1)
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.8, 0.8)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$k_x$')
    ax.set_ylabel(r'$k_y$')
    ax.set_title(f'Exciton {i_exc + 1}')
    
    if fin_mom_exciton:
        ax.arrow(0, 0, Qshift[0], Qshift[1], color='white', width=0.01, head_width=0.04, length_includes_head=True)

    plt.tight_layout()
    if save_pdf:
        pdf.savefig(fig)
    else:
        plt.savefig(f'exciton_{i_exc + 1}_mod_Akcv_summed_cv.png', dpi=300)
    plt.close(fig)
    
    # Plot 2: (c, v) band index amplitude
    fig = plt.figure()
    im = plt.imshow(sum_k_Akcv.T, origin='lower', aspect='auto', cmap='plasma')
    plt.colorbar(im, label=r'$\sum_k |A_{kcv}|^2$')
    plt.xlabel('Conduction')
    plt.ylabel('Valence')
    plt.xticks(np.arange(nc), np.arange(nc)+1)
    plt.yticks(np.arange(nv), np.arange(nv)+1)
    plt.title(f'Exciton {i_exc + 1}')
    plt.tight_layout()

    if save_pdf:
        pdf.savefig(fig)
    else:
        plt.savefig(f'exciton_{i_exc + 1}_mod_Akcv_summed_k.png', dpi=300)
    plt.close(fig)
    
    print('Finished plotting exciton', i_exc + 1, 'of', i_exc_max)

print("Finished")