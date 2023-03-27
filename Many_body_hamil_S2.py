import numpy as np
from construct_basis import construct_basis
from sparce_matrix import *
from create_H_S_array import *
from linear_algebra import *
import time
from print_level import *
from natural_orbitals import natural_Orbitals
from Ras_construct_basis import RAS_construct_basis


def Many_body_hamil_S2(Hfile,Nel,Nst,Sz,max_oc,Processors,S2_spin=False):
    st = time.time()
    Print_input_parameters(Nst,Nel,Sz,max_oc,Processors)

    indsgarb = np.argwhere(np.abs(Hfile[:,0]) < 0.0000001)
    Hfile=np.delete(Hfile,indsgarb,0)
    ifBM_Bm = 0
    Hfile1body, Hfile2body= Hfile_2_one_two_body_Hfile(Hfile)
    S2file1body, S2file2body = one_and_two_body_Sfile(Nst)

    os = 0                     #Default Value
    Hsites = Hfile[:,1:4]
    indoe = np.argwhere(Hsites[:,1]==0)
    Hsitesoe = Hsites[indoe,0]
    # if min(Hsitesoe) < 0:
    #     os = 1
    # else:
    #     os = 0


    if  ifBM_Bm ==1 or ifBM_Bm ==0:
        B,BM, Bm = construct_basis(Nst,Nel,max_oc,Sz,ifBM_Bm)
        #B,BM,Bm=RAS_construct_basis(4,4,10,10,0,excitation=1)
        if os ==0:
            print("\nNow we build the hamiltonian matrix")
            H,HM,Hm = create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst,Processors)
            if S2_spin == True:                
                S2 = create_hamiltonian_sparse_matrix_general_nonreduced(S2file1body,S2file2body,B,BM,Bm,0,Nst,Processors)
        else:
            H,HM,Hm = create_hamiltonian_sparse_matrix_general_openshell(Hfile,B,BM,Bm,ifBM_Bm)

    E,C = eigenval_H(H)
    if S2_spin ==True:
        sval=eigen_S2_value(S2,E,C)
        Print_energies_and_S_value(E,sval)
        et = time.time()
        print("\nWeight of slater determinant of Ground state")
        for i in range(len(B)):
            if np.abs(C[i,0]) > 0.05:
                print(B[i], np.round((C[i,0]*C[i,0]),3), np.round(C[i,0],3))
        print("\nWeight of slater determinant of Excited state")
        for k in range(len(B)):
            if np.abs(C[k,1]) > 0.05:
                print(B[k], np.round((C[k,1]*C[k,1]),3), np.round(C[k,1],3))
        
        Print_time_elepsed(st,et)

    return B,C





