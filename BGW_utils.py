
import numpy as np


def read_sigma(arquivo, nband, kpoint, col_2_get):

    E = 0
    arq = open(arquivo)
    flag_read = False
    flag_found = False

    for line in arq:
        linha = line.split()
        if len(linha) > 0:
            if linha[0] == 'k':
                kx, ky, kz = float(linha[2]), float(linha[3]), float(linha[4])
                if kpoint[0] == kx and kpoint[1] == ky and kpoint[2] == kz:
                    flag_read = True
        if len(linha) == 15:
            if linha[0].isnumeric() is True and flag_read is True:
                if int(linha[0]) == nband:
                    flag_found = True
                    E = float(linha[col_2_get])

        if flag_found is True:
            break

    if flag_found is True:
        return E
    else:
        print(arquivo + ' not found')
        return False


def ajuste_nbands(N, N0, a, b):
    return a/(N-N0) + b


def ajuste_nbands2(Entrada, N0, E0, a, b, Econverg):
    N, Ecut = Entrada
    return a/(N-N0) + b/(Ecut - E0) + Econverg


def read_eqp_data(arq, Nval):

    eqp_data = open(arq)

    Kpoints, Ek_n, Ek_n_mf, BandsInd = [], [], [], []

    for line in eqp_data:
        linha = line.split()
        if linha[0] != '1':
            Kpoints.append(
                np.array([float(linha[0]), float(linha[1]), float(linha[2])]))
            TotBands = int(linha[3])
            Ek_n.append([])
            Ek_n_mf.append([])
        else:
            Ek_n[-1].append(float(linha[3]))
            Ek_n_mf[-1].append(float(linha[2]))
            if len(BandsInd) < TotBands:
                BandsInd.append(int(linha[1]))

    eqp_data.close()

    # Convertendo para array
    Ek_n = np.array(Ek_n)
    Ek_n_mf = np.array(Ek_n_mf)

    # Setting maximum of valence band to 0
    if Nval > 0:
        try:
            #print('Nao achei banda de valencia')
            Nval_index = BandsInd.index(Nval)
            Ef = max(Ek_n[:, Nval_index])
            Ef_mf = max(Ek_n_mf[:, Nval_index])
            print('Colocando maximo da banda de val no 0.0')
        except:
            print('Nao colocando maximo da banda de val no 0.0')
            Ef = 0.0
            Ef_mf = 0.0
    else:
        print('Nao colocando maximo da banda de val no 0.0')
        Ef = 0.0
        Ef_mf = 0.0

    Ek_n = Ek_n - Ef
    Ek_n_mf = Ek_n_mf - Ef_mf

    return Kpoints, Ek_n, BandsInd, Ef, Ek_n_mf, Ef_mf


def kpoint_is_between_Kpoints(K0, Kf, kpoint):
    # checa se kpoint esta entre K0 e Kf

    parallel = False

    # parallel?
    dK = Kf - K0
    dk = kpoint - K0

    norm_dK, norm_dk = np.linalg.norm(dK), np.linalg.norm(dk)

    if norm_dK > 0:
        if norm_dk > 0:
            cosseno = np.dot(dk, dK)/(norm_dK*norm_dk)
        else:
            cosseno = 1.0
    else:
        print('Erro: K0 = Kf. Olhar caminho definido!')
        print(Kf, K0)
        cosseno = 0

    if abs(cosseno - 1) <= 1e-3:
        parallel = True

    if parallel is True:
        if (kpoint >= K0).all() and (Kf >= kpoint).all():
            return True
        elif (kpoint <= K0).all() and (Kf <= kpoint).all():
            return True
        else:
            #print('nao esta entre K0 e Kf')
            return False
    else:
        #print('nao eh paralelo')
        return False


def bandstructure(Path, Ek_n, Kpoints, qshift):

    K, E_K = [], []

    for i in range(len(Ek_n[0])):
        E_K.append([])

    Kpoints_plot = []

    for iPath in range(len(Path) - 1):

        K0 = Path[iPath]
        Kf = Path[iPath + 1]

        print(K0, Kf)

        TempKpoints, TempKdists, TempKInds = [], [], []

        for iKpoints in range(len(Kpoints)):
            kpoint = Kpoints[iKpoints] - qshift
            if kpoint_is_between_Kpoints(K0, Kf, kpoint) is True:
                TempKdists.append(np.linalg.norm(kpoint - K0))
                TempKpoints.append(Kpoints[iKpoints])
                TempKInds.append(iKpoints)

        NewList = list(TempKdists)
        NewList.sort()

        for iList in range(len(NewList)):
            newKpoint = TempKpoints[TempKdists.index(NewList[iList])]
            Kpoints_plot.append(newKpoint)
            NewKInd = TempKInds[TempKdists.index(NewList[iList])]
            if len(Kpoints_plot) == 1:
                K.append(0)
            else:
                dk = np.linalg.norm(Kpoints_plot[-1] - Kpoints_plot[-2])
                K.append(K[-1] + dk)

            for iBands in range(len(Ek_n[0])):
                E_K[iBands].append(Ek_n[NewKInd][iBands])

    return K, E_K


def read_Emf_QE(Emf_file, Nval):

    arq = open(Emf_file)

    flag_get_data = False
    Emf = []

    for line in arq:
        linha = line.split()
        if len(linha) > 0:
            if linha[0] == 'k':
                Emf.append([])
                flag_get_data = True
            elif linha[0] == 'Writing':
                flag_get_data = False
            elif flag_get_data is True:
                for i in range(len(linha)):
                    Emf[-1].append(float(linha[i]))

    # Setting maximum of valence band to 0

    Emf = np.array(Emf)
    Ef = max(Emf[:, Nval])
    Emf = Emf - Ef

    return Emf, Ef


def get_espectro_abs(data_file, norm):
    arq = open(data_file)

    omega, eps1, eps2, JDOS = [], [], [], []

    for line in arq:
        linha = line.split()
        if len(linha) == 4:
            if linha[0] != '#':
                omega.append(float(linha[0]))
                eps2.append(float(linha[1]))
                eps1.append(float(linha[2]))
                JDOS.append(float(linha[3]))

    if norm is True:
        maxEps2 = max(eps2)
    else:
        maxEps2 = 1

    omega, eps1, eps2, JDOS = np.array(omega), np.array(
        eps1), np.array(eps2)/maxEps2, np.array(JDOS)

    return([omega, eps2, eps1, JDOS])


def read_BSE_eigenvalues(eigenvalues_file):

    arq = open(eigenvalues_file)
    return np.loadtxt(arq)


# QE functions
def get_scf_energy(file_out):

    arq = open(file_out)

    for line in arq:
        linha = line.split()

        if len(linha) > 0:
            if linha[0] == '!':
                E = np.float(linha[4])

    return E
