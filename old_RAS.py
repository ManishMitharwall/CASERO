import numpy as np
import math
from itertools import combinations




def RAS_construct_basis(Nst,Nel,block,Sz,ifBM_Bm):
    def nchoosek(n, k):
        return math.comb(n, k)
    Nelmin1= block[0][0];states_deep = block[0][1:]
    state_free = block[1]
    Nelmax3=block[2][0];states_high = block[2][1:]

    BM,Bm = np.zeros([2,2,2]), np.zeros([2,2,2])
    Nstps = Nst
    Nst = 2*Nst
    Nbas = nchoosek(Nst,Nel)
    n= [j for j in range(1,(Nst+1))]
    basisel = list()
    basisel += list(combinations(n, Nel))
    B = []

    for j in range(Nbas):
        vbasis = np.zeros([1,Nst])
        for l in basisel[j]:
            vbasis[0,int(l)-1]=1
        vbasis = np.reshape(vbasis, (2,Nstps),order="F")
        spinz= sum(vbasis[0])- sum(vbasis[1])   
        vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1])
        if spinz==Sz and np.sum(np.abs(vbasis[states_deep]))>= Nelmin1 and np.sum(np.abs((vbasis[states_high])))<= Nelmax3:
            B.append(vbasis)
    B = np.array(B)
    return B, BM, Bm

        

B,BM,Bm=RAS_construct_basis(8,8,[[3,0,1],[2,3,4,5],[1,6,7]],0,0)
print(len(B))


