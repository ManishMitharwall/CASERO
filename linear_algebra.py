import scipy.sparse.linalg as lg
import scipy as sp
import numpy as np

def eigenval_H(H):
    print('Now  Lanczos')
    if H.shape[0] > 10:
        E, C = lg.eigsh(H , k=8, which='SA')
        # print("-----$$$$$$$$$$$$$$$$$$$$$$$")
    else:
        E, C = sp.linalg.eigh(H.toarray())
    return E,C

def eigen_S2_value(S2,E,C):
    eigenS = []
    for i in range(len(E)):
        eigenS += [(C[:,i].T)@(S2@C[:,i])]
    sval =  [0.5*(-1+np.sqrt(1+4*i)) for i in  eigenS]
    return sval