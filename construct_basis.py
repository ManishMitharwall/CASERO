import numpy as np
import math
from itertools import combinations


def construct_basis(Nst,Nel,max_oc,Sz,ifBM_Bm):
    def nchoosek(n, k):
        return math.comb(n, k)

    if not type(max_oc) is list:
        max_oc= [max_oc]

    BM,Bm = np.zeros([2,2,2]), np.zeros([2,2,2])
    Nstps = Nst
    Nst = 2*Nst
    Nbas = nchoosek(Nst,Nel)
    n= [j for j in range(1,(Nst+1))]
    basisel = list()
    basisel += list(combinations(n, Nel))
    B = np.zeros([Nbas,Nstps])
    ind_garb = []

    if len(max_oc)==1:
        for j in range(Nbas):
            vbasis = np.zeros([1,Nst])
            for l in basisel[j]:
                vbasis[0,int(l)-1]=1
            vbasis = np.reshape(vbasis, (2,Nstps),order="F")
            spinz= sum(vbasis[0])- sum(vbasis[1])   
            vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1])
            if spinz==Sz and  max(vbasis) < (sum(max_oc)+1):
                B[j] = vbasis
            else:
                ind_garb.append(j)
        B= np.delete(B,ind_garb,0)

    else:
        Nlevel, Nelec = max_oc[0], max_oc[1]
        for j in range(Nbas):
            vbasis = np.zeros([1,Nst])
            for l in basisel[j]:
                vbasis[0,int(l)-1]=1
            vbasis = np.reshape(vbasis, (2,Nstps),order="F")
            spinz= sum(vbasis[0])- sum(vbasis[1])   
            vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1] )
            if spinz==Sz and  sum(abs(vbasis[0:Nlevel]))==Nelec:
                B[j] = vbasis
            else:
                ind_garb.append(j)
        B= np.delete(B,ind_garb,0)
    
    if ifBM_Bm ==1:
        seed = np.zeros([1,Nst])
        Nbas = nchoosek(Nst,Nel+1)
        n= [j for j in range(1,(Nst+1))]
        basisel = list()
        basisel += list(combinations(n, Nel+1))
        BM = np.zeros([Nbas,Nstps])
        ind_garb = []
        
        for j in range(Nbas):
            vbasis = np.zeros([1,Nst])
            for l in basisel[j]:
                vbasis[0,int(l)-1]=1
            vbasis = np.reshape(vbasis, (2,Nstps),order="F")
            spinz= sum(vbasis[0])- sum(vbasis[1])   
            vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1] )
            if spinz==(Sz+1):  # and  sum(abs(vbasis[0:Nlevel]))==Nelec:
                BM[j] = vbasis
            else:
                ind_garb.append(j)
        BM= np.delete(BM,ind_garb,0)

        seed = np.zeros([1,Nst])
        Nbas = nchoosek(Nst,Nel-1)
        n= [j for j in range(1,(Nst+1))]
        basisel = list()
        basisel += list(combinations(n, Nel-1))
        Bm = np.zeros([Nbas,Nstps])
        ind_garb = []

        for j in range(Nbas):
            vbasis = np.zeros([1,Nst])
            for l in basisel[j]:
                vbasis[0,int(l)-1]=1
            vbasis = np.reshape(vbasis, (2,Nstps),order="F")
            spinz= sum(vbasis[0])- sum(vbasis[1])   
            vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1] )
            if spinz==(Sz-1)  and  max(vbasis) < (sum(max_oc)+1):
                Bm[j] = vbasis
            else:
                ind_garb.append(j)
        Bm= np.delete(Bm,ind_garb,0)
    return B, BM, Bm

        

# print(construct_basis(8,8,2,0,0))
