
import numpy as np
import matplotlib.pyplot as plt

import argparse
parser = argparse.ArgumentParser()

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
    
    print('''Usage : python plot_eqp_file.py -eqp eqp.dat -symm symm_points_file -Nval Nval
          where eqp.dat is the file with the band structure from BGW
          symm_points_file is the file with the high symmetry points
          Nval is the index of last valence band (starting from 1)''')
    
    # -eqp_file eqp.dat file from BGW 
    # -symm_points File with list of high symmetry points
    # -Nval Number of valence band
    parser.add_argument("-eqp", "--eqp_file", help="eqp.dat file from BGW ")
    parser.add_argument("-symm", "--symm_points_file", help="File with list of high symmetry points")
    parser.add_argument("-Nval", "--Nval", help="Number of valence band")

    args = parser.parse_args()

    if args.eqp_file == None:
        print("Using default name for eqp file : eqp.dat")
        eqp_file = "eqp.dat"
    else:
        eqp_file = args.eqp_file
        print(f"eqp file to be read {eqp_file}")
        
    if args.symm_points_file == None:
        print("symm points file not given. Not using it")
        symm_points_file = None    
    else:
        symm_points_file = args.symm_points_file
        print(f"high symmetry points file to be read: {symm_points_file}")    
    
    if args.Nval == None:
        print("Nval not given. Not shifting bands")
        shift_bands = False
    else:
        Nval = int(args.Nval)
        print(f"Setting maximum valence band {Nval} to zero")
        shift_bands = True


    bands_dft, bands_qp, Kpoints, Nk, band_indexes = read_eqp_dat_file(eqp_file)
    K = K_path_from_k_pts_list(Kpoints)

    if shift_bands:
        Nval_index = np.where(band_indexes == Nval)[0][0]
        E0_dft = np.max(bands_dft[Nval_index, :])
        E0_qp = np.max(bands_qp[Nval_index, :])
        bands_dft = bands_dft - E0_dft
        bands_qp = bands_qp - E0_qp
        
    
    if symm_points_file != None:
        
        print('Reading high symmetry K points file')
        with open(symm_points_file) as f:
            lines = f.readlines()
        f.close()
        
        K_symm = []
        K_symm_labels = []
        for line in lines:
            line_split = line.split()
            if len(line_split) == 4:
                kx, ky, kz = float(line_split[1]), float(line_split[2]), float(line_split[3])
                K_symm.append([kx, ky, kz])
                K_symm_labels.append(line_split[0])
                
        K_symm = np.array(K_symm)
        print(f'Number of high symmetry points {len(K_symm)}')
        
        K_symm_dists = K_path_from_k_pts_list(K_symm)                   

    for iband in range(np.shape(bands_dft)[0]):
        plt.plot(K, bands_dft[iband, :], color='gray', alpha = 0.5)
        plt.plot(K, bands_qp[iband, :], color='red', alpha = 1.0)
    
    if symm_points_file != None:
        plt.xticks(K_symm_dists, K_symm_labels)
        plt.xlim([0, K_symm_dists[-1]])
    plt.grid()
    
    plt.show()
  
# plt.xticks(K_path_from_k_pts_list(Path_ortho), Names_ortho)
# plt.xlim([0, K_path_from_k_pts_list(Path_ortho)[-1]])
