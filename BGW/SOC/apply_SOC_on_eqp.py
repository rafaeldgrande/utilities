import numpy as np
import argparse

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
    print('SOC corrections loaded. Original shape:', deltaE_SOC.shape) # shape (Nk, Nbands, Nspin)
    return deltaE_SOC

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Apply SOC corrections to EQP data.')
    parser.add_argument('--eqp_original', default='eqp.dat', help='Input EQP file (default: eqp.dat)')
    parser.add_argument('--corrections_SOC', default='Corrections_SOC.dat', help='SOC corrections file (default: Corrections_SOC.dat)')
    parser.add_argument('--eqp_new', default='eqp_with_SOC.dat', help='Output EQP file with SOC (default: eqp_with_SOC.dat)')
    args = parser.parse_args()

    eqp_original = args.eqp_original
    corrections_SOC = args.corrections_SOC
    eqp_new = args.eqp_new
 
    Edft, Eqp, Kpoints, Nk, bands_indexes = read_eqp_dat_file(eqp_original)
    print('Edft shape:', Edft.shape)
    print('Eqp shape:', Eqp.shape)
    print('Kpoints shape:', Kpoints.shape)
    print('Nk:', Nk)
    print('bands_indexes shape:', bands_indexes.shape)
    # print('bands_indexes from function:', bands_indexes)
    Nbnds = len(bands_indexes)
 
    nmin = int(np.min(bands_indexes))
    nmax = int(np.max(bands_indexes))

    # loading SOC corrections
    data_SOC = get_SOC_corrections(corrections_SOC) # shape (Nk, Nbands, Nspin)
    deltaE_SOC = data_SOC[:, nmin-1:nmax, :]
    
    # reshape deltaE_SOC from (Nk, Nbands, Nspin) to (Nbands, Nk, Nspin)
    deltaE_SOC = np.transpose(deltaE_SOC, (1, 0, 2))
    
    print('SOC corrections shape after slicing and transposing:', deltaE_SOC.shape)

    Nbnds_SOC_file, Nk_SOC_file, _ = deltaE_SOC.shape
    
    # give error if Nbnds_SOC doesnt match Nbnds
    if Nbnds_SOC_file != Nbnds:
        raise ValueError(f"Number of bands in SOC file ({Nbnds_SOC_file}) does not match original EQP file ({Nbnds}).")
    if Nk_SOC_file != Nk:
        raise ValueError(f"Number of k-points in SOC file ({Nk_SOC_file}) does not match original EQP file ({Nk}).")
    
    # applying SOC corrections
    Edft_SOC_1 = Edft + deltaE_SOC[:, :, 0]  # first spin component
    Edft_SOC_2 = Edft + deltaE_SOC[:, :, 1]  # second spin component
    Eqp_SOC_1 = Eqp + deltaE_SOC[:, :, 0]  # first spin component
    Eqp_SOC_2 = Eqp + deltaE_SOC[:, :, 1]  # second spin component

    # writing SOC corrections
    arq_with_SOC = open(eqp_new, 'w')

    for ik in range(Nk):
        arq_with_SOC.write(f"    {Kpoints[ik,0]:.8f}   {Kpoints[ik,1]:.8f}    {Kpoints[ik,2]:.8f}     {Nbnds} \n")

        for ibnd in range(Nbnds):
            band_index_1 = int(bands_indexes[ibnd] * 2)
            band_index_2 = band_index_1 + 1
            arq_with_SOC.write(f"1    {band_index_1}   {Edft_SOC_1[ibnd, ik]:.9f}   {Eqp_SOC_1[ibnd, ik]:.9f} \n")
            arq_with_SOC.write(f"1    {band_index_2}   {Edft_SOC_2[ibnd, ik]:.9f}   {Eqp_SOC_2[ibnd, ik]:.9f} \n")

    arq_with_SOC.close()
    print("Finished!")

