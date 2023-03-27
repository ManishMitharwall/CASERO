import numpy as np
from Pauli_matrix_Hfile import Pauli_matrix_Hfile
from product_of_two_1e_Hfile import product_Of_two_1e_Hfile
from clean_Hfile import clean_Hfile

def readHfile(Hfile_file):
    with open(Hfile_file,'r') as f1:
        Hfile = []
        for line in f1:
            if not line.strip():
                break
            else:
                split = line.split()
                Hfile.append([float(split[0]),int(float(split[1])),int(float(split[2])),int(float(split[3])),int(float(split[4]))])
    return np.array(Hfile)

def Hfile_2_one_two_body_Hfile(Hfile):
    Hfile1body = Hfile[Hfile[:,2]==0]
    Hfile2body = Hfile[Hfile[:,2]!=0]
    return Hfile1body, Hfile2body

def one_and_two_body_Sfile(Nst):
    SxT,SyT, SzT = [],[] ,[]
    for site in range(1,(Nst+1)):
        Sx1,Sy1,Sz1 = Pauli_matrix_Hfile(site)
        SxT += Sx1; SyT += Sy1; SzT += Sz1 
    SxT, SyT, SzT = np.array(SxT), np.array(SyT), np.array(SzT)
    S2file= product_Of_two_1e_Hfile(SxT,SxT)
    S2file = np.array(S2file, dtype=np.complex_)
    S2file=np.vstack((S2file,product_Of_two_1e_Hfile(SyT,SyT)))
    S2file=np.vstack((S2file,product_Of_two_1e_Hfile(SzT,SzT)))
    S2file = S2file.real
    S2file = clean_Hfile(S2file)
    S2file1body = S2file[S2file[:,2]==0]
    S2file2body = S2file[S2file[:,2]!=0]
    return S2file1body, S2file2body

def scissoring_opertor(Hfile,Nel,Nst,Sz,c_value):
    occ_orb = int((Nel+Sz)/2)
    for j in range(occ_orb+1,Nst+1):
        add_line = [c_value,j,0,0,j]
        Hfile = np.vstack((Hfile,add_line))
    return(Hfile)


