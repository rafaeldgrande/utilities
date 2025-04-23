
import numpy as np

# Nkx, Nky, Nkz, qx, qy, qz passed as command line arguments
# Usage: python gera_k_grid.py Nkx Nky Nky qx qy qz

import sys
if len(sys.argv) != 6:
    print("Usage: python gera_k_grid.py Nkx Nky Nkz qx qy qz")
    sys.exit(1)
Nkx = int(sys.argv[1])
Nky = int(sys.argv[2])
Nkz = int(sys.argv[3])
qx = float(sys.argv[4])
qy = float(sys.argv[5])
qz = float(sys.argv[6])

# Nkx, Nky, Nkz = 27, 27, 1
# qx, qy, qz = 0, 0, 0

Nk_tot = Nkx * Nky * Nkz

print("K_POINTS crystal")
print(f"{Nk_tot}")

for ik in range(Nkx):
    for jk in range(Nky):
        for kk in range(Nkz):
            print(f"{(ik/Nkx+qx):.8f}  {(jk/Nky+qy):.8f}  {(kk/Nkz+qz):.8f}  1.0")
