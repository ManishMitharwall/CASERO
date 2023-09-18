import numpy as np
import scipy as sp
import scipy.sparse.linalg as lg
from itertools import combinations
from multiprocess import Pool
import math
import time


def Print_input_parameters(Nst,Nel,Sz,max_oc,Processors):
    if type(Sz) == int:
        Sz = Sz % 2
    print('\n*******************INPUT  PARAMETERS*************************')
    print(f'Number of processor used is                     {Processors}')
    print(f'Maximum occupation of each site is              {max_oc}')
    print(f'Number of sites are                             {Nst}')
    print(f'Total number of electrons are                   {Nel}')
    print(f'Sz value is                                     {Sz}')
    print('*************************************************************')

def Print_energies_and_S_value(E,sval):
    print('\n\tEnergy(ev)\t Spin')
    for i in range(len(E)):
        print('\t%.4f \t %.2f ' % (E[i],sval[i]))
    print('Energy gap of ground state and excited state is     %.4f ' %(E[1]-E[0]))

def Print_time_elepsed(st,et):
    hours, rem = divmod(et-st, 3600)
    minutes, seconds = divmod(rem, 60)
    print("\nTime Taken: {:0>2}H:{:0>2}M:{:05.2f}S\n".format(int(hours),int(minutes),seconds))

def Print_Title():
    title = '''
+---------------------------------------------------------------------------------+
|                              D S P M M                                          |
|                           ================                                      |
|                        A      CASCI      CODE                                   |
+---------------------------------------------------------------------------------+'''
    print(title)

def eigenval_H(H):
    # print('Now  Lanczos')
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

def construct_basis(Nst,Nel,max_oc,Sz):
    def nchoosek(n, k):
        return math.comb(n, k)
    if not type(Sz) is list:
        Sz=[Sz]
        
    Nstps = Nst
    Nst = 2*Nst
    Nbas = nchoosek(Nst,Nel)
    n= [j for j in range(1,(Nst+1))]
    basisel = list()
    basisel += list(combinations(n, Nel))
    B = np.zeros([Nbas,Nstps])
    ind_garb = []
    for j in range(Nbas):
        vbasis = np.zeros([1,Nst])
        for l in basisel[j]:
            vbasis[0,int(l)-1]=1
        vbasis = np.reshape(vbasis, (2,Nstps),order="F")
        spinz= sum(vbasis[0])- sum(vbasis[1])   
        vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1])
        if spinz in Sz and  ( max(vbasis) < (max_oc+1) ):
            B[j] = vbasis
        else:
            ind_garb.append(j)
    B= np.delete(B,ind_garb,0)
    return B.astype(int)

def Hfile_2_one_two_body_Hfile(Hfile):
    Hfile1body = Hfile[Hfile[:,2]==0]
    Hfile2body = Hfile[Hfile[:,2]!=0]
    return Hfile1body, Hfile2body

def Pauli_matrix_Hfile(site):
    Sx = [[0.5 , site, 0, 0, -site],[0.5, -site, 0, 0, site]]
    Sy = [[-0.5*1j, site, 0, 0, -site], [0.5*1j, -site, 0, 0, site]]
    Sz = [[0.5 , site, 0, 0, site], [-0.5, -site, 0, 0, -site]]
    return(Sx,Sy,Sz)

def product_Of_two_1e_Hfile(O1,O2):
    n1 = O1.shape[0]; n2 = O2.shape[0]; O3 = []
    for l1 in range(n1):
        w1 = O1[l1,0]; c1= O1[l1,1]; d1= O1[l1,4]
        for l2 in range(n2):
            w2 = O2[l2,0]; c2= O2[l2,1]; d2= O2[l2,4]
            if d1 != c2 and d1 != d2:
                O3_new =  [w1*w2 , c1 , c2 , d2 , d1]
                O3.append(O3_new)
            if d1 == c2:
                O3_new =  [w1*w2 ,c1, 0 , 0, d2]
                O3.append(O3_new)
                if d1 != d2:
                    O3_new = [-1*(w1*w2), c1, c2, d1 , d2]
                    O3.append(O3_new)
    return O3

def clean_Hfile(Hfile):
    Hfile = np.array(Hfile)
    indsgarb = np.argwhere(np.abs(Hfile[:,0]) < 0.00000001)
    Hfile=np.delete(Hfile,indsgarb,0)
    a=0
    while a < len(Hfile):
        kk = np.where((Hfile[a,1:]==Hfile[:,1:]).all(axis=1))
        kk = np.concatenate(kk)
        if len(kk)>=2:
            coeff = np.sum(Hfile[kk,0],axis=0)
            Hfile[a,0]=coeff
            Hfile= np.delete(Hfile,kk[1:],axis=0)
        a+=1
    indsgarb1 = np.argwhere(np.abs(Hfile[:,0]) < 0.00000001)
    Hfile=np.delete(Hfile,indsgarb1,0)
    return(Hfile)

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

def apply_hamiltonian_1body_nonreduced(Ri, Hfile, Nst):
    R = []
    try:
        Hfile = np.array(Hfile)
        nrow_h, ncol_h = Hfile.shape
    except:
        Hfile= [Hfile]
        Hfile = np.array(Hfile)
        nrow_h, ncol_h = Hfile.shape

    try:
        Ri = np.array(Ri)
        nrow_R, ncol_R = Ri.shape
    except:
        Ri= [Ri]
        Ri = np.array(Ri)
        nrow_R, ncol_R = Ri.shape
    
    for a in range(nrow_h):
        i = int(abs(Hfile[a,1]))
        l = int(abs(Hfile[a,4]))

        si = np.sign(Hfile[a,1])
        sl = np.sign(Hfile[a,4])

        if si == -1:
            i = i + Nst
        if sl == -1:
            l = l + Nst
       
        for el in range(nrow_R):
            b0 = Ri[el,1:]
            b= b0
            V1 = np.array(b)
            V1[V1==-1]=0;V1[V1==2]=1
            V2 = np.array(b)
            V2[V2==1]=0;V2[V2==2]=1
            V2= np.abs(V2)
            V = np.concatenate((V1,V2), axis=0)
            V, signo = apply_operator_1body(i,l,V)
            if len(V) !=0:
                V1,V2 =  V[0:Nst], np.multiply(0.5,V[Nst:]) 
                V= np.add(V1,V2)
                V[V==1.5]=2;V[V==0.5]=-1
                R_new = np.insert(V,0,(Ri[el,0]*Hfile[a,0]*signo))                     
                R.append(R_new)       
    return R

def apply_operator_1body(i,l,V):
    signo = 1
    # destroy l
    signo = signo*((-1)**(sum(V[:l-1])))
    V[l-1] = 0*(V[l-1]==1) + 2*(V[l-1] != 1)
    # create i
    signo = signo*((-1)**(sum(V[:i-1])))
    V[i-1]= 1*(V[i-1]==0) + 2*(V[i-1] !=0)
    if 2 in V or -2 in V:
        V = []
    return V, signo


def apply_operator_2body(i,j,k,l,V):
    signo = 1
    # destroy l
    signo = signo*((-1)**(sum(V[:l-1])))
    V[l-1] = 0*(V[l-1]==1) + 2*(V[l-1] != 1)
    # destro k
    signo = signo*((-1)**(sum(V[:k-1])))
    V[k-1] = 0*(V[k-1]==1) + 2*(V[k-1] != 1)
    # create j
    signo = signo*((-1)**(sum(V[:j-1])))
    V[j-1]= 1*(V[j-1]==0) + 2*(V[j-1] !=0)
    # create i
    signo = signo*((-1)**(sum(V[:i-1])))
    V[i-1]= 1*(V[i-1]==0) + 2*(V[i-1] !=0)
    if 2 in V or -2 in V:
        V = []
    return V, signo

def apply_hamiltonian_1body(Ri, Hfile, Nst):
    R = []
    try:
        Ri = np.array(Ri)
        nrow_R,ncol_R = Ri.shape
    except:
        Ri= [Ri]
        Ri = np.array(Ri)
        nrow_R, ncol_R = Ri.shape      
    try:
        Hfile = np.array(Hfile)
        nrow_h, ncol_h = Hfile.shape
    except:
        Hfile= [Hfile]
        Hfile = np.array(Hfile)
        nrow_h, ncol_h = Hfile.shape
    
    for a in range(nrow_h):
        i = int(np.abs(Hfile[a,1]))
        l = int(np.abs(Hfile[a,4]))

        si = np.sign(Hfile[a,1])
        sl = np.sign(Hfile[a,4])

        if si == -1:
            i = i + Nst
        if sl == -1:
            l = l + Nst

        for el in range(nrow_R):
            b0 = Ri[el,1:]
            b= b0
            V1 =b
            V1=[0 if x==-1 else x for x in V1]
            V1=[1 if x==2 else x for x in V1]

            V2 = b
            V2=[0 if x==1 else x for x in V2]
            V2=[1 if x==2 else x for x in V2]
            V2 = [abs(x) for x in V2]
            V = V1 +V2
            V, signo = apply_operator_1body(i,l,V)

            if len(V) !=0:
                V1,V2 =  V[0:Nst], np.multiply(0.5,V[Nst:]) 
                V= np.add(V1,V2)
                V = [2 if x==1.5 else x for x in V]
                V = [-1 if x==0.5 else x for x in V]
                R_new = [Ri[el,0]*Hfile[a,0]*signo] +V                     
                R.append(R_new)

            if i !=l:
                V1 = b
                V1 = [0 if x ==-1 else x for x in V1]
                V1 = [1 if x == 2 else x for x in V1]
                V2 = b
                V2 = [0 if x ==1 else x for x in V2]
                V2 = [1 if x == 2 else x for x in V2]
                V2 = [abs(x) for x in V2]
                V = V1 + V2
                V, signo = apply_operator_1body(l,i,V)

                if len(V) !=0:
                    V1,V2 =  V[0:Nst], np.multiply(0.5,V[Nst:])
                    V= np.add(V1,V2)
                    V = [2 if x==1.5 else x for x in V]
                    V = [-1 if x==0.5 else x for x in V]
                    R_new = [Ri[el,0]*Hfile[a,0]*signo] + V
                    R.append(R_new)
            i = i + Nst
            l = l + Nst
            V1 =b
            V1=[0 if x==-1 else x for x in V1]
            V1=[1 if x==2 else x for x in V1]

            V2 = b
            V2=[0 if x==1 else x for x in V2]
            V2=[1 if x==2 else x for x in V2]
            V2 = [abs(x) for x in V2]
            V = V1 +V2
            V, signo = apply_operator_1body(i,l,V)

            if len(V) !=0:
                V1,V2 =  V[0:Nst], np.multiply(0.5,V[Nst:]) 
                V= np.add(V1,V2)
                V = [2 if x==1.5 else x for x in V]
                V = [-1 if x==0.5 else x for x in V]
                R_new = [Ri[el,0]*Hfile[a,0]*signo] +V
                R.append(R_new)

            if i !=l:
                V1 = b
                V1 = [0 if x ==-1 else x for x in V1]
                V1 = [1 if x == 2 else x for x in V1]
                
                V2 = b
                V2 = [0 if x ==1 else x for x in V2]
                V2 = [1 if x == 2 else x for x in V2]
                V2 = [abs(x) for x in V2]
                V = V1 + V2
                V, signo = apply_operator_1body(l,i,V)

                if len(V) !=0:
                    V1,V2 =  V[0:Nst], np.multiply(0.5,V[Nst:])
                    V= np.add(V1,V2)
                    V = [2 if x==1.5 else x for x in V]
                    V = [-1 if x==0.5 else x for x in V]
                    R_new = [Ri[el,0]*Hfile[a,0]*signo] + V
                    R.append(R_new)
    return R

def apply_operator_2body(i,j,k,l,V):
    signo = 1
    # destroy l
    signo = signo*((-1)**(sum(V[:l-1])))
    V[l-1] = 2*(V[l-1] != 1)    #  0*(V[l-1]==1) + 2*(V[l-1] != 1)
    # destro k
    signo = signo*((-1)**(sum(V[:k-1])))
    V[k-1] = 2*(V[k-1] != 1)    #  0*(V[k-1]==1) + 2*(V[k-1] != 1)
    # create j
    signo = signo*((-1)**(sum(V[:j-1])))
    V[j-1]= 1*(V[j-1]==0) + 2*(V[j-1] !=0)
    # create i
    signo = signo*((-1)**(sum(V[:i-1])))
    V[i-1]= 1*(V[i-1]==0) + 2*(V[i-1] !=0)

    if 2 in V or -2 in V:
        V = []
    return V, signo


def apply_hamiltonian_2body(Ri, Hfile, Nst):
    R = []  ## state after hamiltonian is applied

    V1_in = [0 if x == -1 else 1 if x == 2 else x for x in Ri[1:]]
    V2_in = [0 if x == 1 else 1 if x == 2 else abs(x) for x in Ri[1:]]
    basis_V = V1_in + V2_in  
    R_new = [] 
    for  row in Hfile:
        i, j, k, l = int(row[1]), int(row[2]), int(row[3]), int(row[4])
        if i < 0: i = abs(i) + Nst
        if j < 0: j = abs(j) + Nst
        if k < 0: k = abs(k) + Nst
        if l < 0: l = abs(l) + Nst
        V_out, signo = apply_operator_2body(i, j, k, l, basis_V.copy())
        if V_out:
            V1_out, V2_out = V_out[:Nst], np.multiply(0.5, V_out[Nst:])
            V_out = [2 if x == 1.5 else -1 if x == 0.5 else x for x in np.add(V1_out, V2_out)]
            R_new = [Ri[0] * row[0] * signo] + V_out
            R.append(R_new)
    return R


def create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,Nst,Processors):
    global call, call_m, call_M
    B = np.array(B)
    nr,nc,nv = [],[],[]   
    n = [x for x in range(len(B))]
    print(len(n))
    def call(j):
        Ri = np.insert(B[j,:],0,1)
        R1= apply_hamiltonian_1body(Ri,Hfile1body,Nst)
        R2 = apply_hamiltonian_2body(Ri,Hfile2body,Nst)
        R = R1 + R2
        R = np.array(R)
        nn = len(R)
        j1,rw1,RR = [],[],[]
        for inn in range(nn):
            R_new = [R[inn,1:]]
            row_indx = np.where((B==R_new).all(axis=1))
            row_indx=np.concatenate(row_indx)
            if len(row_indx)!=0:
                j1.append(j); rw1+=list(row_indx); RR += [R[inn,0]]
        return j1, rw1,RR

    pool = Pool(Processors)
    x=[]
    x.append(pool.map(call,n))
    x = np.array(x,dtype=object)
    nr = x[:,:,0]; nc = x[:,:,1]; nv = x[:,:,2]; 
    nr = np.concatenate(nr);nr = np.concatenate(nr);nr= list(nr)
    nc = np.concatenate(nc);nc = np.concatenate(nc);nc = list(nc)
    nv = np.concatenate(nv);nv = np.concatenate(nv); nv = list(nv)
    nr = np.array(nr).astype(int)
    nc = np.array(nc).astype(int)
    H = sp.sparse.csc_matrix((nv,(nr,nc)))
    H.eliminate_zeros()
    return H


def create_hamiltonian_sparse_matrix_general_nonreduced(Hfile1body,Hfile2body,B,Nst,Processors):
    global call
    B = np.array(B)
    nr,nc,nv = [],[],[]    
    n = [x for x in range(len(B))]
    def call(j):
        Ri = np.insert(B[j,:],0,1)
        R1= apply_hamiltonian_1body_nonreduced(Ri,Hfile1body,Nst)
        R2 = apply_hamiltonian_2body(Ri,Hfile2body,Nst)
        R = R1 + R2
        R = np.array(R)
        nn = len(R)
        j1,rw1,RR = [],[],[]
        for inn in range(nn):
            R_new = [R[inn,1:]]
            row_indx = np.where((B==R_new).all(axis=1))
            row_indx=np.concatenate(row_indx)
            if len(row_indx)!=0:
                j1.append(j); rw1+=list(row_indx); RR += [R[inn,0]]
        return j1, rw1,RR

    x=[]
    pool = Pool(Processors)
    x.append(pool.map(call,n))
    x = np.array(x,dtype=object)
    nr = x[:,:,0]; nc = x[:,:,1]; nv = x[:,:,2]
    nr = np.concatenate(nr);nr = np.concatenate(nr);nr= list(nr)
    nc = np.concatenate(nc);nc = np.concatenate(nc);nc = list(nc)
    nv = np.concatenate(nv);nv = np.concatenate(nv); nv = list(nv)
    H = sp.sparse.csc_matrix((nv,(nr,nc)))
    H.eliminate_zeros()
    return H


def Many_body_hamil_S2(Hfile,Nel,Nst,Sz,max_oc,Processors,S2_spin=False):
    st = time.time()
    Print_input_parameters(Nst,Nel,Sz,max_oc,Processors)
    indsgarb = np.argwhere(np.abs(Hfile[:,0]) < 0.0000001)
    Hfile=np.delete(Hfile,indsgarb,0)
    Hfile1body, Hfile2body= Hfile_2_one_two_body_Hfile(Hfile)
    S2file1body, S2file2body = one_and_two_body_Sfile(Nst)
    B = construct_basis(Nst,Nel,max_oc,Sz)
    #B,BM,Bm=RAS_construct_basis(4,4,10,10,0,excitation=1)

    # print("\nNow we build the hamiltonian matrix")
    H = create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,Nst,Processors)
    if S2_spin == True:                
        S2 = create_hamiltonian_sparse_matrix_general_nonreduced(S2file1body,S2file2body,B,Nst,Processors)

    E,C = eigenval_H(H)
    if S2_spin ==True:
        sval=eigen_S2_value(S2,E,C)
        Print_energies_and_S_value(E,sval)
        et = time.time()
        print("\nWeight of slater determinant of Ground state")
        
        for i in range(len(B)):
            if np.abs(C[i,0]) > 0.05:
                print(f'{str(round( C[i,0]**2,3)):<5}{str(round(C[i,0],3)):>15}{str(B[i]):>30}')
        print("\nWeight of slater determinant of Excited state")
        for k in range(len(B)):
            if np.abs(C[k,1]) > 0.05:
                print(f'{str(round( C[k,0]**2,3)):<5}{str(round(C[k,0],3)):>15}{str(B[k]):>30}')     
        Print_time_elepsed(st,et)
    return B,C,E

def Many_body_hamil_S2_DOS(Hfile,Nel,Nst,Sz,max_oc,Processors,S2_spin=False):
    st = time.time()
    Print_input_parameters(Nst,Nel,Sz,max_oc,Processors)
    indsgarb = np.argwhere(np.abs(Hfile[:,0]) < 0.0000001)
    Hfile=np.delete(Hfile,indsgarb,0)
    Hfile1body, Hfile2body= Hfile_2_one_two_body_Hfile(Hfile)
    S2file1body, S2file2body = one_and_two_body_Sfile(Nst)
    B = construct_basis(Nst,Nel,max_oc,Sz)

    H = create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,Nst,Processors)
    if S2_spin == True:                
        S2 = create_hamiltonian_sparse_matrix_general_nonreduced(S2file1body,S2file2body,B,Nst,Processors)
    if H.shape[0] > 1000:
        E, C = lg.eigsh(H , k=2000, which='SA')
        print("SPARD")
    else:
        E,C = sp.linalg.eigh(H.toarray())

    if S2_spin ==True:
        sval=eigen_S2_value(S2,E,C)
        Print_energies_and_S_value(E,sval)
        et = time.time()
        print("\nWeight of slater determinant of Ground state")
        
        for i in range(len(B)):
            if np.abs(C[i,0]) > 0.05:
                print(f'{str(round( C[i,0]**2,3)):<5}{str(round(C[i,0],3)):>15}{str(B[i]):>30}')
        print("\nWeight of slater determinant of Excited state")
        for k in range(len(B)):
            if np.abs(C[k,1]) > 0.05:
                print(f'{str(round( C[k,0]**2,3)):<5}{str(round(C[k,0],3)):>15}{str(B[k]):>30}')     
        Print_time_elepsed(st,et)
    return B,C,E


# Hfile = np.array([[-3,1,0,0,2],[4,1,-1,-1,1],[4,2,-2,-2,2]])
# Nst = int(np.max(Hfile[:,1:]))
# Nel = Nst
# Sz=0
# max_oc=2
# P=1

# B0,C0,E0=Many_body_hamil_S2(Hfile,Nel,Nst,[-2,1,0,1,2],max_oc,P,S2_spin=1)
