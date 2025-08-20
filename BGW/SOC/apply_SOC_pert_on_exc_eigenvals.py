
import numpy as np
import argparse
import h5py

def get_SOC_corrections(corrections_SOC):

    deltaE_SOC = []
    arq = open(corrections_SOC)
    for line in arq:
        line_split = line.split()
        if line_split[0] == 'kpoint':
            deltaE_SOC.append([[]])
        else:
            if len(deltaE_SOC[-1][-1]) < 2:
                deltaE_SOC[-1][-1].append(float(line_split[1]))
            else:
                deltaE_SOC[-1].append([float(line_split[1])])
    arq.close()
    deltaE_SOC = np.array(deltaE_SOC)
    print('SOC corrections loaded:', deltaE_SOC.shape) # shape (Nk, Nbands, Nspin)
    return deltaE_SOC

def get_exciton_info(exciton_file):

    """    
    Return the exciton energy and the eigenvec coefficients Acvk

    Assuming calculations with TD approximation
    Info about file at: http://manual.berkeleygw.org/3.0/eigenvectors_h5_spec/
    Also, just working for excitons with Q = 0 and one spin
    
    Parameters:
    exciton_file = exciton file name (string). ex: eigenvecs.h5
    iexc = Exciton index to be read
    
    Returns:
    Acvk = Exciton wavefunc coefficients. array Akcv[ik, ic, iv] with complex values
    Omega = Exciton energy (BSE eigenvalue) in eV (float)
    """

    print('Reading exciton info from file', exciton_file)

    f_hdf5 = h5py.File(exciton_file, 'r')
    
    flavor_calc = f_hdf5['/exciton_header/flavor'][()]
    eigenvecs   = f_hdf5['exciton_data/eigenvectors'][()]          # (nQ, Nevecs, nk, nc, nv, ns, real or imag part)
    eigenvals   = f_hdf5['exciton_data/eigenvalues'][()] 
    ifmax       = f_hdf5['/mf_header/kpoints/ifmax'][()]  
    
    # print('ifmax:', ifmax)
    # print('ifmax shape:', ifmax.shape)
    nval_index = np.max(ifmax)
    
    if flavor_calc == 2:
        print('Flavor in BGW: complex')
        Akcv = eigenvecs[0,:,:,:,:,0,0] + 1.0j*eigenvecs[0,:,:,:,:,0,1]
    else:
        print('Flavor in BGW: real')
        Akcv = eigenvecs[0,:,:,:,:,0,0]
    # shape of Akcv is (nexc, nk, nc, nv)

    # K points in the fine grid
    Kpoints_bse = f_hdf5['/exciton_header/kpoints/kpts'][()] 
    Nkpoints = f_hdf5['/exciton_header/kpoints/nk'][()]

    return Akcv, Kpoints_bse, Nkpoints, nval_index

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Apply SOC corrections to excitons eigenvalues.')
    parser.add_argument('--corrections_SOC', default='Corrections_SOC.dat', help='SOC corrections file (default: Corrections_SOC.dat)')
    parser.add_argument('--eigenvals_file', default='eigenvalues.dat', help='eigenvalues file output from absorption step of BerkeleyGW (default: eigenvalues.dat)')
    parser.add_argument('--eigenvals_with_SOC', default='eigenvalues_with_SOC.dat', help='Output eigenvalues file with SOC (default: eigenvalues_with_SOC.dat)')
    parser.add_argument('--eigenvectors_file', default='eigenvectors.h5', help='eigenvectors file (default: eigenvectors.h5)')
    parser.add_argument('--plot_soc', type=lambda x: x.lower() == 'true', default=False, help='Plot SOC perturbation on eigenvalues (default: False)')
    args = parser.parse_args()

    corrections_SOC = args.corrections_SOC
    eigenvals_file = args.eigenvals_file
    eigenvals_with_SOC = args.eigenvals_with_SOC
    eigenvectors_file = args.eigenvectors_file
    plot_soc = args.plot_soc
    
    if plot_soc:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
    
    # loading SOC corrections
    deltaE_SOC = get_SOC_corrections(corrections_SOC) # shape (Nk, Nbands, Nspin)

    # loading original eigenvalues data
    eigenvals_data = np.loadtxt(eigenvals_file)

    eigenvals = eigenvals_data[:, 0]  
    dip_moments = eigenvals_data[:, 1:]  # not perturbed by SOC
    
    # loading exciton info
    Akcv, Kpoints_bse, Nkpoints, nval_index = get_exciton_info(eigenvectors_file)
    
    Nexc, Nk, Nc, Nv = Akcv.shape
    print('Akcv shape:', Akcv.shape)
    print('Nk:', Nk, 'Nc:', Nc, 'Nv:', Nv)
    print('nval_index:', nval_index)
    
    bmin = nval_index - Nv 
    bmax = nval_index + Nc
    
    # print('bmin:', bmin, 'bmax:', bmax)
    
    delta_E_SOC_val = deltaE_SOC[:, bmin:nval_index, :][:, ::-1, :] # shape (Nk, Nv, Nspin) - reversed valence bands
    delta_E_SOC_cond = deltaE_SOC[:, nval_index:bmax, :] # shape (Nk, Nc, Nspin)

    delta_E_SOC_kcv = np.zeros((Nk, Nc, Nv, 2))
    
    for ik in range(Nk):
        for ic in range(Nc):
            for iv in range(Nv):
                delta_E_SOC_kcv[ik, ic, iv, :] = delta_E_SOC_cond[ik, ic, :] - delta_E_SOC_val[ik, iv, :]
                
    # print('delta_E_SO_kcv shape:', delta_E_SO_kcv.shape)

    # print('delta_E_SOC_val shape:', delta_E_SOC_val.shape)
    # print('delta_E_SOC_cond shape:', delta_E_SOC_cond.shape)

    # apply SOC corrections to eigenvalues
    
    eigenvales_up, eigenvals_down = np.zeros((Nexc)), np.zeros((Nexc))
    for iexc in range(Nexc):
        pert_up = np.sum(delta_E_SOC_kcv[:, :, :, 0] * abs(Akcv[iexc, :, :, :])**2)
        pert_down = np.sum(delta_E_SOC_kcv[:, :, :, 1] * abs(Akcv[iexc, :, :, :])**2)
        eigenvales_up[iexc] = eigenvals[iexc] + pert_up
        eigenvals_down[iexc] = eigenvals[iexc] + pert_down
        
        if plot_soc:
            plt.plot(eigenvals[iexc], pert_up, 'ro')
            plt.plot(eigenvals[iexc], pert_down, 'bo')
            
    if plot_soc:
        plt.xlabel('Exciton index')
        plt.ylabel('SOC Perturbation (eV)')
        plt.title('SOC Perturbation on Exciton Eigenvalues')
        plt.legend(['Up Spin', 'Down Spin'])
        plt.xlim([1.5, 3.0])
        plt.savefig('SOC_perturbation.png')

    # print('eigenvales_up shape:', eigenvales_up.shape)
    # print('first 5 eigenvales_up:', eigenvales_up[:5])
    # print('eigenvals_down shape:', eigenvals_down.shape)
    # print('first 5 eigenvals_down:', eigenvals_down[:5])
    
    arq = open(eigenvals_with_SOC, 'w')
    for iexc in range(Nexc):
        arq.write(f"{eigenvales_up[iexc]:.9f}   {dip_moments[iexc, 0]:.9f}   {dip_moments[iexc, 1]:.9f}   {dip_moments[iexc, 2]:.9f}\n")
        arq.write(f"{eigenvals_down[iexc]:.9f}   {dip_moments[iexc, 0]:.9f}   {dip_moments[iexc, 1]:.9f}   {dip_moments[iexc, 2]:.9f}\n")
    arq.close()
    print('Finished! Eigenvalues with SOC saved to', eigenvals_with_SOC)
