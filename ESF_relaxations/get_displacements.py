
import numpy as np


final_pos = np.loadtxt('final_pos', usecols=(1,2,3)).flatten()
initial_pos = np.loadtxt('initial_pos', usecols=(1,2,3)).flatten()

displacement = final_pos - initial_pos
np.savetxt('displacement', displacement)