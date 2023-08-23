
import numpy as np

Nkx, Nky, Nkz = 4,4,4
qx, qy, qz = 0, 0, 0

Nk_tot = Nkx * Nky * Nkz

print("K_POINTS crystal")
print(f"{Nk_tot}")

for ik in range(Nkx):
    for jk in range(Nky):
        for kk in range(Nkz):
            print(f"{(ik/Nkx+qx):.8f}  {(jk/Nky+qy):.8f}  {(kk/Nkz+qz):.8f}  1.0")