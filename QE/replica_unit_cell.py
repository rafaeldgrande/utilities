
import numpy as np

arq_old = 'in'
arq_new = 'in_replicated'
N1, N2, N3 = 1, 1, 4
Nreplicas = N1*N2*N3

atomic_pos = []
text = ''

arq_data = open(arq_old)

text_input = arq_data.readlines()

# Modifying number of atom and lattice vectors

for i_line in range(len(text_input)):
    linha = text_input[i_line].split()
    if len(linha) > 0:
        if linha[0] == 'nat':
            Nat = int(linha[2])
            Nat_new = Nat*Nreplicas
            text_input[i_line] = f'    nat = {Nat_new}'+'\n'
        if linha[0] == 'CELL_PARAMETERS':
            linha1 = text_input[i_line + 1].split()
            linha2 = text_input[i_line + 2].split()
            linha3 = text_input[i_line + 3].split()
            a1 = np.array([float(linha1[0]), float(linha1[1]), float(linha1[2])])
            a2 = np.array([float(linha2[0]), float(linha2[1]), float(linha2[2])])   
            a3 = np.array([float(linha3[0]), float(linha3[1]), float(linha3[2])])

            print('Old lattice vectors')
            print('a1', a1)
            print('a2', a2)
            print('a3', a3)

            text_input[i_line + 1] = f' {a1[0]*N1}  {a1[1]*N1}   {a1[2]*N1}'+'\n'
            text_input[i_line + 2] = f' {a2[0]*N2}  {a2[1]*N2}   {a2[2]*N2}'+'\n'
            text_input[i_line + 3] = f' {a3[0]*N3}  {a3[1]*N3}   {a3[2]*N3}'+'\n'

            print('New lattice vectors')
            print('a1', text_input[i_line + 1])
            print('a2', text_input[i_line + 2])
            print('a3', text_input[i_line + 3])


# Reading atomic posistions

Atomic_pos = []

for i_line in range(len(text_input)):
    linha = text_input[i_line].split()
    if len(linha) == 4:
        pos = np.array([ float(linha[1]), float(linha[2]), float(linha[3]) ])
        Atomic_pos.append([linha[0], pos])

arq_data.close()

# Writing new file
arq_saida = open(arq_new, 'w')

for i_line in range(len(text_input)):

    if text_input[i_line].split()[0] == 'K_POINTS': # Now I will write the new atomic positions
        print('hello!')

        for i1 in range(N1):
            for i2 in range(N2):
                for i3 in range(N3):
                    print('here' ,i1, i2, i3)
                    if i1 != 0 or i2 != 0 or i3 != 0:
                        displacement = a1*i1 + a2*i2 + a3*i3
                        print(displacement, 'hi')
                        for i_atom in range(len(Atomic_pos)):
                            pos_new = Atomic_pos[i_atom][1] + displacement
                            text_append = Atomic_pos[i_atom][0] + f' {pos_new[0]}  {pos_new[1]}   {pos_new[2]}'+'\n'
                            arq_saida.write(text_append)



    arq_saida.write(text_input[i_line])
arq_saida.close()

print('Done!')

