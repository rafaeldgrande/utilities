
import numpy as np
import matplotlib.pyplot as plt

import argparse
parser = argparse.ArgumentParser()

Nval = 13

def dists_symm(Path):
    
    Dists_symm = [0.0]
    for i in range(1, len(Path)):
        dr = np.linalg.norm( np.array(Path[i]) - np.array(Path[i-1]) )
        Dists_symm.append(Dists_symm[-1] + dr)
        
    return Dists_symm


def read_eqp_dat_file(eqp_file):
    
    bands_dft, bands_qp = [], []
    
    # loading file
    data = np.loadtxt(eqp_file)

    # first getting number of bands in this file. first line is:   0.000000000  0.000000000  0.000000000      13
    Nbnds = int(data[0, 3])
    
    # getting list of band indexes
    band_indexes = data[1:Nbnds+1, 1]
    
    # now we get the k points in this file
    Kpoints = data[0::Nbnds+1] # get lines 0, Nbnds+1, 2*(Nbnds+2), ...
    Kpoints = Kpoints[:, :3] # remove last collumn with 
    
    Nk = len(Kpoints)
    print(f'Number of kpoints {Nk}')

    for ibnd in range(Nbnds):
        temp = data[ibnd+1::Nbnds+1]
        bands_dft.append(temp[:, 2])
        bands_qp.append(temp[:, 3])
                
    return np.array(bands_dft), np.array(bands_qp), Kpoints, Nk, band_indexes

def K_path_from_k_pts_list(k_pts_list):

    K = [0.0]

    for i in range(len(k_pts_list) - 1):
        dK = np.linalg.norm( np.array(k_pts_list[i + 1]) - np.array(k_pts_list[i]) )
        K.append(K[-1] + dK)
        
    return np.array(K)

if __name__ == "__main__":

    # Gettting parameters
    
    # -eqp_file eqp.dat file from BGW 
    # -symm_points File with list of high symmetry points
    parser.add_argument("-eqp", "--eqp_file", help="eqp.dat file from BGW ")
    parser.add_argument("-symm", "--symm_points_file", help="File with list of high symmetry points")

    args = parser.parse_args()

    # print( "eqp.dat file  {} \nsymm points file {}".format(
    #         args.eqp_file,
    #         args.symm_points_file
    #         ))

    if args.eqp_file == None:
        print("Using default name for eqp file : eqp.dat")
        eqp_file = "eqp.dat"
    else:
        eqp_file = args.eqp_file
        print(f"eqp file to be read {eqp_file}")
        
    if args.symm_points_file == None:
        print("Using default name for high symmetry points file: symm_points_file")
        symm_points_file = "symm_points_file"    
    else:
        symm_points_file = args.symm_points_file
        print(f"high symmetry points file to be read: {symm_points_file}")    


    bands_dft, bands_qp, Kpoints, Nk, band_indexes = read_eqp_dat_file(eqp_file)
    K = K_path_from_k_pts_list(Kpoints)

    Nval_index = np.where(band_indexes == Nval)[0][0]
    E0_dft = np.max(bands_dft[Nval_index, :])
    E0_qp = np.max(bands_qp[Nval_index, :])
    bands_dft = bands_dft - E0_dft
    bands_qp = bands_qp - E0_qp

# Gamma = [0, 0, 0]
# X = [1/2, 0, 0]
# Y = [0, 1/2, 0]
# Z = [0, 0, 1/2]
# M = [1/2, 1/2, 0]
# N = [0, 1/2, 1/2]
# O = [1/2, 0, 1/2]
# R = [1/2, 1/2, 1/2]

# # Caminho R-Gamma-X-M-Gamma
# Path_tetra = [Gamma, X, M, R, O, Z, Gamma]
# Names_tetra = [r'$\mathrm{\Gamma}$', 
#          r'$\mathrm{X}$',
#          r'$\mathrm{M}$',
#          r'$\mathrm{R}$',
#          r'$\mathrm{O}$',
#          r'$\mathrm{Z}$',
#          r'$\mathrm{\Gamma}$']                    

    for iband in range(np.shape(bands_dft)[0]):
        plt.plot(K, bands_dft[iband, :], color='gray', alpha = 0.5)
        plt.plot(K, bands_qp[iband, :], color='red', alpha = 1.0)
    plt.show()
  
# plt.xticks(K_path_from_k_pts_list(Path_ortho), Names_ortho)
# plt.xlim([0, K_path_from_k_pts_list(Path_ortho)[-1]])