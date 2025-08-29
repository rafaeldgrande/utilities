

'''
This code reads two outputs from QE calculations, one using SOC with FR pseudopotentials
and other with spins using SR pseudopotentials (SOC is missing). Here I am calculating the
correction to be used later in GW calculations

Egw_FR = Egw_SR + Edft_FR - Edft_SR

I am assuming the order of k points in the two files is the same. Also assuming it 
is the same order in GW part.

Usage:
python get_SOC_corrections.py --qe_out_SOC_on FR/mos2_bands_soc_fine_grid.dat --qe_out_SOC_off SR/mos2_bands_soc_fine_grid.dat --nval 52

'''

import numpy as np
import argparse


def get_bands_and_spin(bands_file):
    
    print('Reading QE bands output:', bands_file)
    
    with open(bands_file, 'r') as file:
        
        # read first line
        first_line = file.readline()
        # print(first_line.split())
        line_split = first_line.split()
        # print(line_split[2].split(','))
        #  &plot nbnd=  60, nks=    86 /
        nbnds = int(line_split[2].split(',')[0])
        nk = int(line_split[4])
        
        print('Number of bands:', nbnds)
        print('Number of k-points:', nk)
        
        data = np.zeros((nk, nbnds), dtype=float)
        
        ik = -1
        while True:
            line = file.readline()
            line_split = line.split()
            if len(line_split) == 3:
                ik += 1
                ib = 0
            else:
                for val in line_split:
                    data[ik, ib] = float(val)
                    ib += 1
                    # print(ik, ib, bands[ik, ib-1])
            if not line:
                break
    
    return data # shape nk, nb


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Calculate SOC corrections from QE band outputs')
    parser.add_argument('--qe_out_SOC_on', default='FR/mos2_bands_soc_fine_grid.dat', help='Path to QE bands output with FR pseudopotentials (SOC on)')
    parser.add_argument('--qe_out_SOC_off', default='SR/mos2_bands_soc_fine_grid.dat', help='Path to QE bands output with SR pseudopotentials (SOC off)')
    parser.add_argument('--spins_out_SOC_on', default='FR/mos2_bands_soc.dat.3', help='Path to QE spins output with FR pseudopotentials (SOC on)')
    parser.add_argument('--spins_out_SOC_off', default='SR/mos2_bands_soc.dat.3', help='Path to QE spins output with SR pseudopotentials (SOC off)')
    parser.add_argument('--nval', type=int, default=52, help='Index of the valence band. For calculations with spin nval = number of electrons. For calculations without spin nval = number of electrons / 2')
    parser.add_argument('--map_spin', type=bool, default=True, help='Whether to map the spins when calculating the corrections')
    args = parser.parse_args()
    
    qe_out_SOC_on = args.qe_out_SOC_on
    qe_out_SOC_off = args.qe_out_SOC_off
    spins_out_SOC_on = args.spins_out_SOC_on
    spins_out_SOC_off = args.spins_out_SOC_off
    nval = args.nval - 1
    map_spin = args.map_spin
    
    # qe_out_SOC_on = "FR/mos2_bands_soc.dat"
    # qe_out_SOC_off = "SR/mos2_bands_soc.dat"
    # spins_out_SOC_on = "FR/mos2_bands_soc.dat.3"
    # spins_out_SOC_off = "SR/mos2_bands_soc.dat.3"

    bands_SOC_on = get_bands_and_spin(qe_out_SOC_on)
    bands_SOC_off = get_bands_and_spin(qe_out_SOC_off)
    
    # here I will approximate the correction in Gamma at the val band to be zero
    # also assuming the first k point is (0, 0, 0)
    bands_SOC_on -= bands_SOC_on[0, nval]
    bands_SOC_off -= bands_SOC_off[0, nval]

    spins_SOC_on = get_bands_and_spin(spins_out_SOC_on)
    spins_SOC_off = get_bands_and_spin(spins_out_SOC_off)

    # print(spins_SOC_on[0, :10])
    # print(spins_SOC_off[0, :10])
    
    Nk, Nbnds = bands_SOC_on.shape
    
    if map_spin:
        corrections = np.zeros_like(bands_SOC_on) # shape nk, nb

        for ik in range(Nk):
            for ib in range(int(Nbnds/2)):
                ib1 = 2*ib
                ib2 = 2*ib + 1
                
                if spins_out_SOC_on[ik, ib1] <= 0.0:
                    E_SOC_on_spin_down = bands_SOC_on[ik, ib1] 
                    E_SOC_on_spin_up = bands_SOC_on[ik, ib2]
                else:
                    E_SOC_on_spin_down = bands_SOC_on[ik, ib2]
                    E_SOC_on_spin_up = bands_SOC_on[ik, ib1]
                    
                if spins_out_SOC_off[ik, ib1] <= 0.0:
                    E_SOC_off_spin_down = bands_SOC_off[ik, ib1]
                    E_SOC_off_spin_up = bands_SOC_off[ik, ib2]
                else:
                    E_SOC_off_spin_down = bands_SOC_off[ik, ib2]
                    E_SOC_off_spin_up = bands_SOC_off[ik, ib1]

                corrections[ik, ib1] = E_SOC_on_spin_down - E_SOC_off_spin_down
                corrections[ik, ib2] = E_SOC_on_spin_up - E_SOC_off_spin_up

            #     # check if signal of spins_SOC_on[ik, ib1] is the same as spins_SOC_off[ik, ib1]
            #     if np.sign(spins_SOC_on[ik, ib1]) == np.sign(spins_SOC_off[ik, ib1]):
            #         corrections[ik, ib1] = bands_SOC_on[ik, ib1] - bands_SOC_off[ik, ib1]
            #         corrections[ik, ib2] = bands_SOC_on[ik, ib2] - bands_SOC_off[ik, ib2]
            #     # else:
            #     #     corrections[ik, ib1] = bands_SOC_on[ik, ib2] - bands_SOC_off[ik, ib2]
            #     #     corrections[ik, ib2] = bands_SOC_on[ik, ib1] - bands_SOC_off[ik, ib1]
                    
    else:
        corrections = bands_SOC_on - bands_SOC_off # shape nk, nb

    # checking!
    print('Gamma point')
    for ib in range(nval-5, nval+5):
        print(f'Band {ib+1:3d}   SOC on: {bands_SOC_on[0, ib]:8.4f}   SOC off: {bands_SOC_off[0, ib]:8.4f}   Correction: {corrections[0, ib]:8.4f}')
        
    # print('Kpoint ik=120')
    # for ib in range(nval-4, nval+4):
    #     print(f'Band {ib+1:3d}   SOC on: {bands_SOC_on[120, ib]:8.4f}   SOC off: {bands_SOC_off[120, ib]:8.4f}   Correction: {corrections[120, ib]:8.4f}')

    print('Writing SOC corrections to file: Corrections_SOC.dat')
    arq_correction = open('Corrections_SOC.dat', 'w')
    
    for ik in range(Nk):
        arq_correction.write(f'kpoint {ik+1}\n')
        for ibnd in range(Nbnds):
            arq_correction.write(f'{ibnd+1}   {corrections[ik, ibnd]:.6f}\n')

    arq_correction.close()
    
    print('Finished!')


if __name__ == "__main__":
    main()
