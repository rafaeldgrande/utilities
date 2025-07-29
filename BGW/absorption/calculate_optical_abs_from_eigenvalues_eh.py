
import numpy as np
import matplotlib.pyplot as plt
import sys

"""Usage ./plot_eigenvalues.py eigenvalues_file gaussian_broadening output_file"""

# Read the eigenvalues from the file
eigenvalues_file = sys.argv[1]
sigma = float(sys.argv[2])
output_file = sys.argv[3]

data = np.loadtxt(eigenvalues_file)
eigvals = data[:,0]
Nexcitons = len(eigvals)
abs_dip2 = data[:,1]

Emax = np.max(eigvals)
Energies = np.arange(0.0, Emax, 0.01)
Absorption = np.zeros(len(Energies))

def gaussian(x, mu, sigma):
    return np.exp(-0.5*((x-mu)/sigma)**2)/(np.sqrt(2*np.pi)*sigma) - np.exp(-0.5*((x+mu)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)

for iexc in range(Nexcitons):
    Absorption += abs_dip2[iexc]*gaussian(Energies, eigvals[iexc], sigma)


np.savetxt(output_file, np.column_stack((Energies, Absorption)))

plt.plot(Energies, Absorption)
plt.savefig('absorption.png')

