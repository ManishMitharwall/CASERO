import numpy as np
from apply_hamiltonian_1body import apply_hamiltonian_1body
from apply_hamiltonian_1body_nonreduced import apply_hamiltonian_1body_nonreduced
from apply_hamiltonian_general_openshell import apply_hamiltonian_general_openshell
from apply_hamiltonian_2body import apply_hamiltonian_2body
from itertools import chain
from scipy import sparse
from multiprocess import Pool


def create_hamiltonian_sparse_matrix_general(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst,Processors):
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
    H = sparse.csc_matrix((nv,(nr,nc)))
    H.eliminate_zeros()


    if ifBM_Bm ==1:
        Bm = np.array(Bm)
        nr,nc,nv = [],[],[]   
        n = [x for x in range(len(Bm))]
        def call_m(j):
            Ri = np.insert(Bm[j,:],0,1)
            R1= apply_hamiltonian_1body(Ri,Hfile1body,Nst)
            R2 = apply_hamiltonian_2body(Ri,Hfile2body,Nst)
            R = R1 + R2
            R = np.array(R)
            nn = len(R)
            j1,rw1,RR = [],[],[]
            for inn in range(nn):
                R_new = [R[inn,1:]]
                row_indx = np.where((B==R_new).all(axis=1))
                j1.append(j); rw1+= row_indx; RR += [R[inn,0]]
            return j1, rw1,RR

        pool = Pool(Processors)
        x=[]
        x.append(pool.map(call_m,n))
        x = np.array(x,dtype=object)
        nr = x[:,:,0]; nc = x[:,:,1]; nv = x[:,:,2]; 
        nr = np.concatenate(nr);nr = np.concatenate(nr);nr= list(nr)
        nc = np.concatenate(nc);nc = np.concatenate(nc);nc = list(nc)
        nv = np.concatenate(nv);nv = np.concatenate(nv); nv = list(nv)
        nc = list(chain.from_iterable(nc))
        Hm = sparse.csc_matrix((nv,(nr,nc)))
        Hm.eliminate_zeros()


        BM = np.array(BM)
        nr,nc,nv = [],[],[]   
        n = [x for x in range(len(BM))]
        def call_M(j):
            Ri = np.insert(BM[j,:],0,1)
            R1= apply_hamiltonian_1body(Ri,Hfile1body,Nst)
            R2 = apply_hamiltonian_2body(Ri,Hfile2body,Nst)
            R = R1 + R2
            R = np.array(R)
            nn = len(R)
            j1,rw1,RR = [],[],[]
            for inn in range(nn):
                R_new = [R[inn,1:]]
                row_indx = np.where((B==R_new).all(axis=1))
                j1.append(j); rw1+= row_indx; RR += [R[inn,0]]
            return j1, rw1,RR
        pool = Pool(Processors)
        x=[]
        x.append(pool.map(call_M,n))
        x = np.array(x,dtype=object)
        nr = x[:,:,0]; nc = x[:,:,1]; nv = x[:,:,2]; 
        nr = np.concatenate(nr);nr = np.concatenate(nr);nr= list(nr)
        nc = np.concatenate(nc);nc = np.concatenate(nc);nc = list(nc)
        nv = np.concatenate(nv);nv = np.concatenate(nv); nv = list(nv)
        nc = list(chain.from_iterable(nc))
        HM = sparse.csc_matrix((nv,(nr,nc)))
        HM.eliminate_zeros()
    else:
        Hm,HM = np.zeros(2),np.zeros(2)
    return H, HM, Hm          


def create_hamiltonian_sparse_matrix_general_nonreduced(Hfile1body,Hfile2body,B,BM,Bm,ifBM_Bm,Nst,Processors):
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
    H = sparse.csc_matrix((nv,(nr,nc)))
    H.eliminate_zeros()


    if ifBM_Bm ==1:
        try:
            Bm_row, Bm_col = Bm.shape
        except:
            Bm = np.array(Bm)
            Bm_row, Bm_col = Bm.shape

        nr,nc,nv = [],[],[]
        for j in range(Bm_row):
            Ri = np.insert(Bm[j,:],0,1)
            R1= apply_hamiltonian_1body_nonreduced(Ri,Hfile1body,Nst)
            R2 = apply_hamiltonian_2body(Ri,Hfile2body,Nst)
            R = R1 + R2
            R = np.array(R)
            nn = len(R)

            for inn in range(nn):
                R_new = [R[inn,1:]]
                row_indx = np.where((B==R_new).all(axis=1))
                j1 +=[j]; rw1+= row_indx; RR+=[R[inn,0]]
        nc = list(chain.from_iterable(nc))
        Hm = sparse.csc_matrix((nv,(nr,nc)))
        Hm.eliminate_zeros()


        BM = np.array(BM)


        nr,nc,nv = [],[],[]
        for j in range(len(BM)):
            Ri = np.insert(BM[j,:],0,1)
            R1= apply_hamiltonian_1body_nonreduced(Ri,Hfile1body,Nst)
            R2 = apply_hamiltonian_2body(Ri,Hfile2body,Nst)
            R = R1 + R2
            R = np.array(R)
            nn = len(R)
            for inn in range(nn):
                R_new = [R[inn,1:]]
                row_indx = np.where((B==R_new).all(axis=1))
                j1 +=[j]; rw1+= row_indx; RR+=[R[inn,0]]
        nc = list(chain.from_iterable(nc))
        HM = sparse.csc_matrix((nv,(nr,nc)))
        HM.eliminate_zeros()
    else:
        Hm,HM = np.zeros(2),np.zeros(2)
  
    return H 


def create_hamiltonian_sparse_matrix_general_openshell(Hfile,B,BM,Bm,ifBM_Bm):
    try:
        B_row, B_col = B.shape
    except:
        B = np.array(B)
        B_row, B_col = B.shape

    nr,nc,nv = [],[],[]    
    for j in range(B_row):
        Ri = np.insert(B[j,:],0,1)
        R= apply_hamiltonian_general_openshell(Ri,Hfile)
        R = np.array(R)
        nn = len(R)
        for inn in range(nn):
            R_new = [R[inn,1:]]
            row_indx = np.where((B==R_new).all(axis=1))
            j1 +=[j]; rw1+= row_indx; RR+=[R[inn,0]]
    nc = list(chain.from_iterable(nc))
    H = sparse.csc_matrix((nv,(nr,nc)))
    H.eliminate_zeros()


    if ifBM_Bm ==1:
        try:
            Bm_row, Bm_col = Bm.shape
        except:
            Bm = np.array(Bm)

            Bm_row, Bm_col = Bm.shape

        nr,nc,nv = [],[],[]
        for j in range(Bm_row):
            Ri = np.insert(Bm[j,:],0,1)
            R= apply_hamiltonian_general_openshell(Ri,Hfile)
            R = np.array(R)
            nn = len(R)

            for inn in range(nn):
                R_new = [R[inn,1:]]
                row_indx = np.where((B==R_new).all(axis=1))
                j1 +=[j]; rw1+= row_indx; RR+=[R[inn,0]]
        nc = list(chain.from_iterable(nc))
        Hm = sparse.csc_matrix((nv,(nr,nc)))
        Hm.eliminate_zeros()

        try:
            BM_row, BM_col = BM.shape
        except:
            BM = np.array(BM)

            BM_row, BM_col = BM.shape

        nr,nc,nv = [],[],[]
        for j in range(BM_row):
            Ri = np.insert(BM[j,:],0,1)
            R= apply_hamiltonian_general_openshell(Ri,Hfile)
            R = np.array(R)
            nn = len(R)
            for inn in range(nn):
                R_new = [R[inn,1:]]
                row_indx = np.where((B==R_new).all(axis=1))
                j1 +=[j]; rw1+= row_indx; RR+=[R[inn,0]]
        nc = list(chain.from_iterable(nc))
        HM = sparse.csc_matrix((nv,(nr,nc)))
        HM.eliminate_zeros()
    else:
        Hm,HM = np.zeros(2),np.zeros(2)
  
    return H, HM, Hm                  


