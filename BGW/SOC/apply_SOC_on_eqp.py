
import numpy as np
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Apply SOC corrections to EQP data.')
	parser.add_argument('--eqp_original', default='eqp.dat', help='Input EQP file (default: eqp.dat)')
	parser.add_argument('--corrections_SOC', default='Corrections_SOC.dat', help='SOC corrections file (default: Corrections_SOC.dat)')
	parser.add_argument('--eqp_new', default='eqp_with_SOC.dat', help='Output EQP file with SOC (default: eqp_with_SOC.dat)')
	args = parser.parse_args()

	eqp_original = args.eqp_original
	corrections_SOC = args.corrections_SOC
	eqp_new = args.eqp_new

	# loading SOC corrections
	deltaE_SOC = []
	arq = open(corrections_SOC)
	for line in arq:
		line_split = line.split()
		if line_split[0] == 'kpoint':
			deltaE_SOC.append([])
		else:
			deltaE_SOC[-1].append(float(line_split[1]))
	arq.close()
	deltaE_SOC = np.array(deltaE_SOC)

	# loading original EQP data
	data = np.loadtxt(eqp_original)

	# figuring out the number of bands
	Nbnds = int(data[0,3])

	# figuring out the number of k points
	Nk = int(len(data)/(Nbnds+1))

	# list of k points
	Kpoints = data[::(Nbnds+1),:3]

	# list of bands indexes
	bands_indexes = data[1:(Nbnds+1):,1]
 
 	# writing SOC corrections
	arq_with_SOC = open(eqp_new, 'w')

	for ik in range(Nk):
		arq_with_SOC.write(f"{Kpoints[ik,0]:.8f}   {Kpoints[ik,1]:.8f}    {Kpoints[ik,2]:.8f}   {Nbnds} \n")

		for ibnd in range(Nbnds):
			Edft = data[ik * (Nbnds+1) + ibnd + 1, 2] + deltaE_SOC[ik, int(bands_indexes[ibnd])-1] 
			Eqp  = data[ik * (Nbnds+1) + ibnd + 1, 3] + deltaE_SOC[ik, int(bands_indexes[ibnd])-1]
			arq_with_SOC.write(f"1    {int(bands_indexes[ibnd])}   {Edft:.9f}   {Eqp:.9f} \n")

	arq_with_SOC.close()
	print("Finished!")
