

import numpy as np
import matplotlib.pyplot as plt
import sys 

''' Plots data from N files in one plot. Made it to be quick 
    so the plot is not fancy (does not have labels, etc)
    
    Usage:
    
    python input1 input2 ... inputN figure_name.png
    
    Where it will plot the data from inputs file and save it in figure 
    figure_name.png.
    
    Each input must be two collums of values. The data does not be in ascending order
    (the code fix that).'''

config_dir='/mnt/d/OneDrive/0-Projects/Codes/utilities/config_files/presentation.mplstyle'
plt.style.use(config_dir)


Ninputs = len(sys.argv) - 2
out_name = sys.argv[-1]

for i_input in range(Ninputs):

    input_name = sys.argv[i_input + 1]

    data = np.loadtxt(input_name)

    arr_inds = data[:,0].argsort()
    sorted_x = data[:,0][arr_inds]
    sorted_y = data[:,1][arr_inds]

    #plt.plot(data[:, 0], data[:, 1], 'o')
    plt.plot(sorted_x, sorted_y, '-o', label=input_name)

plt.legend(fontsize=10)
plt.savefig(out_name)

