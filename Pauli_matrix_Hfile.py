import numpy as np

def Pauli_matrix_Hfile(site):
    Sx = [[0.5 , site, 0, 0, -site],[0.5, -site, 0, 0, site]]
    Sy = [[-0.5*1j, site, 0, 0, -site], [0.5*1j, -site, 0, 0, site]]
    Sz = [[0.5 , site, 0, 0, site], [-0.5, -site, 0, 0, -site]]
    return(Sx,Sy,Sz)


