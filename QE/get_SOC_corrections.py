

'''
This code reads two outputs from QE calculations, one using SOC with FR pseudopotentials
and other with spins using SR pseudopotentials (SOC is missing). Here I am calculating the
correction to be used later in GW calculations

Egw_FR = Egw_SR + Edft_FR - Edft_SR

I am assuming the order of k points in the two files is the same. Also assuming it 
is the same order in GW part.

Usage:
python get_SOC_corrections.py qe_out_SOC_on qe_out_SOC_off

'''

import numpy as np
import argparse


def get_bands(qe_out):

	bands = []

	with open(qe_out, 'r') as file:
		for line in file:
			line_split = line.split()
			if len(line_split) > 0:
				if line_split[0] == 'k':
					bands.append([])
				elif line_split[0] == 'Writing':
					break	
				else:
					if len(bands) > 0:
						for value in line_split:
							bands[-1].append(float(value)) 

	return np.array(bands)


def main():
	# Set up argument parser
	parser = argparse.ArgumentParser(description='Calculate SOC corrections from QE band outputs')
	parser.add_argument('qe_out_SOC_on', help='Path to QE bands output with FR pseudopotentials (SOC off)')
	parser.add_argument('qe_out_SOC_off', help='Path to QE bands output with SR pseudopotentials (SOC on)')

	args = parser.parse_args()

	qe_out_SOC_on = args.qe_out_SOC_on
	qe_out_SOC_off = args.qe_out_SOC_off

	bands_SOC_on = get_bands(qe_out_SOC_on)
	bands_SOC_off = get_bands(qe_out_SOC_off)

	correction = bands_SOC_on - bands_SOC_off

	arq_correction = open('Corrections_SOC.dat', 'w')

	for ik in range(len(correction)):
		arq_correction.write(f'kpoint {ik+1}\n')
		for ibnd in range(len(correction[0])):
			arq_correction.write(f'{ibnd+1}   {correction[ik, ibnd]}\n')

	arq_correction.close()

	print('Finished!')


if __name__ == "__main__":
	main()
