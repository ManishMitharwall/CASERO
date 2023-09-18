import numpy as np

def create_anihilate_e(Cg, atom, spin, B, BM, Bm):
    # Create/anihilate electron
    Cgc = np.zeros(len(BM))
    Cga = np.zeros(len(Bm))
    for j in range(len(Cg)):
        S = B[j, :].copy()
        k = S.copy()
        # Create
        if S[atom] == -1*spin:
            S[atom] = 2
            row_indx = np.where((BM == S).all(axis=1))[0]
            if row_indx.size != 0:
                Cgc[row_indx] = Cg[j]

        elif S[atom] == 0:
            S[atom] = spin
            row_indx = np.where((BM == S).all(axis=1))[0]
            if row_indx.size != 0:
                Cgc[row_indx] = Cg[j]

        # Anihilate
        elif S[atom] == spin:
            S[atom] = 0
            row_indx = np.where((Bm == S).all(axis=1))[0]
            if row_indx.size != 0:
                Cga[row_indx] = Cg[j]

        elif S[atom] == 2:
            S[atom] = -1*spin
            row_indx = np.where((Bm == S).all(axis=1))[0]
            if row_indx.size != 0:
                Cga[row_indx] = Cg[j]

    return Cgc, Cga

import numpy as np

def spectral_function(C, E, B, CM, EM, BM, Cm, Em, Bm, sites, eta, w, T, States_DOS_max, Charged_states_max):
    kb = 0.000086173324  # ev/K
    beta = 1 / (kb * 1)
    h = 0.6582  # ev*fs
    DOS = np.zeros_like(w)

    if T > 0:
        # set up Boltzman factors
        p = np.zeros(States_DOS_max)
        Q = 0
        for n in range(States_DOS_max):  # loop over states of the neutral system
            p[n] = np.exp(-beta * E[n, 0])
            Q += p[n]
        p /= Q
        for n in range(States_DOS_max):  # loop over states of the neutral system
            for j in sites:  # loop over sites of active space
                Cc, Ca = create_anihilate_e(C[:, n], abs(j)-1, np.sign(j), B, BM, Bm)
                
                for m in range(Charged_states_max):  # loop over states of the +1 charged system
                    DOS += p[n] * eta * (np.dot(CM[:, m], Cc) ** 2) / ((w - (EM[m, 0] - E[n, 0])) ** 2 + eta * eta)

                    DOS += p[n] * eta * (np.dot(Cm[:, m], Ca) ** 2) / ((w - (E[n, 0] - Em[m, 0])) ** 2 + eta * eta)

    else:  # if T == 0
        # set up Boltzman factors
        indsGS = np.where(np.abs(E[:] - E[0]) < 0.00001)[0]  # indices of Ground State
        for n in range(len(indsGS)):  # loop over states of the neutral system
            for j in sites:  # loop over sites of active space
                Cc, Ca = create_anihilate_e(C[:, n], abs(j)-1, np.sign(j), B, BM, Bm)
                for m in range(Charged_states_max):  # loop over states of the +1 charged system
                    DOS += eta * (np.dot(CM[:, m], Cc) ** 2) / ((w - (EM[m] - E[n])) ** 2 + eta * eta)
                    DOS += eta * (np.dot(Cm[:, m], Ca) ** 2) / ((w - (E[n] - Em[m])) ** 2 + eta * eta)
                    
        DOS /= len(indsGS)

    return DOS



