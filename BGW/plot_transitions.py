#!/usr/bin/python

'''Usage: 
python plot_transitions.py file.dat Nk Nc Nv Emin Emax
or 
python plot_transitions.py file.dat

Where 
file.dat = Name of eigenvalues file (can be with or without eh interaction)
Nc   = Number of cond bands to be analized 
Nv   = Number of val bands to be analized 
Nk   = Number of k points to be analized 
Emin = Just analyze eigenvalues with energy > Emin. Value in eV
Emax = Just analyze eigenvalues with energy < Emax. Value in eV 


To get all valence (conduction) bands set Nc (Nv) to a negative number.
If 
'''

import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
pdf_report = matplotlib.backends.backend_pdf.PdfPages("output.pdf")

# config_dir='/mnt/d/OneDrive/0-Projects/Codes/utilities/config_files/presentation.mplstyle'
# plt.style.use(config_dir)

Ry2eV = 13.605662285137
bohr2A = 0.529177

plot_each_transition = True  # if True plot abs(k) vs k for each one of v->c transitions

# Default parameters
file_name = 'eigenvalues_noeh.dat'

try:
    file_name = sys.argv[1]
    Nk        = int(sys.argv[2])
    Nc        = int(sys.argv[3])
    Nv        = int(sys.argv[4])
    Emin      = float(sys.argv[5])
    Emax      = float(sys.argv[6])
except:
    print('Usage -> python plot_transitions.py file.dat Nk Nc Nv Emin Emax')
    print('Using default parameters: nc=-1, nv=-1, and nk=-1')
    print('Using default Emin and Emax (1e5 may be enough!) (all values in file)')

    Nc, Nv, Nk = -1, -1, -1
    Emin = 0
    Emax = 1e5 # may be enough!
    file_name = 'eigenvalues_noeh.dat'

# get information in file 
print(f'Reading file {file_name}')
arq = open(file_name)
for line in arq:
    linha = line.split()
    if linha[1] == 'nspin,':  # ex: # nspin, nspinor, nkpt, ncb, nvb =        1       1      96      11      10
        Nk_file = int(linha[-3])
        Nc_file = int(linha[-2])
        Nv_file = int(linha[-1])
        break
arq.close()

# compare file information with input parameters



def compare_input_with_eigenvals_file(N_input, N_file, Name_var):
    N_out = N_input
    if N_input > N_file:
        print(f'Requested {Name_var}({N_input}) greater than {Name_var} in file ({N_file}). Making {Name_var} to be {N_file} now.')
        N_out = N_file
    elif N_input < 0:
        print(f'{Name_var}({N_input}) < 0, then setting {Name_var} to be {N_file}.')
        N_out = N_file
    return N_out

Nk = compare_input_with_eigenvals_file(Nk, Nk_file, 'Nk')
Nc = compare_input_with_eigenvals_file(Nc, Nc_file, 'Nc')
Nv = compare_input_with_eigenvals_file(Nv, Nv_file, 'Nv')

# K indexes
K = [k for k in range(Nk)]

# prepare list with indexes (v, c)
transitions = []
for iv in range(Nv):
    for ic in range(Nc):
        transitions.append((iv + 1, ic + 1))

# abs(dip) for each (v,c) pair
abs_vs_k = []
for i_transition in range(len(transitions)):
    abs_vs_k.append([])

# reading data
data_eigenvals = np.loadtxt(file_name)

for i_line in range(len(data_eigenvals)):

    ik = int(data_eigenvals[i_line, 0])
    ic = int(data_eigenvals[i_line, 1])
    iv = int(data_eigenvals[i_line, 2])

    delta_E = data_eigenvals[i_line, 6]
    if Emin < delta_E < Emax:
        abs_dip = data_eigenvals[i_line, 7]
    else:
        abs_dip = 0.0

    if ic <= Nc and iv <= Nv and ik <= Nk:
        i_transition = transitions.index((iv, ic))
        abs_vs_k[i_transition].append(abs_dip)

data_eigenvals = np.array(data_eigenvals)

# Integrated abs(dip)**2 over k

integrated_absdip = np.zeros((Nc, Nv))

for ic in range(1, Nc + 1):
    for iv in range(1, Nv + 1):
        i_transition = transitions.index((iv, ic))
        integrated_absdip[ic-1, iv-1] = np.sum(abs_vs_k[i_transition])


# Plotting integrated abs**2 over k points for each ic,iv pair in matrix form

plt.figure()
plt.matshow((integrated_absdip))
plt.colorbar()
plt.xlabel('Valence')
plt.ylabel('Conduction')
pdf_report.savefig()


# Plotting abs**2 vs k for each ic,iv pair

if plot_each_transition == True:
    for i_transition in range(len(transitions)):
        plt.figure()
        
        iv, ic = transitions[i_transition]
        plt.title(f'ic={ic} / iv={iv}')
        plt.plot(K, abs_vs_k[i_transition])
        plt.xlabel('k index')
        plt.ylabel('abs(dip)**2')

        pdf_report.savefig()
        plt.close()


pdf_report.close()

print('Finished!')
