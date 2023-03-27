import numpy as np
import scipy as sp
from itertools import combinations
from operator import itemgetter
from multiprocess import Pool
from Orbital_tools import write_molden
import time
from print_level import *


def natural_Orbitals(B,C,Nst,Nel,P):
    st = time.time()
    global call
    def call(se):
        Rho = np.zeros([2*Nst,2*Nst])
        seq1 = np.insert(combss[se,:],0,j); seq2 = np.insert(combss[se,:],0,k)
        iseq1,seq1 = zip(*sorted(enumerate(seq1), key=itemgetter(1)))
        iseq2,seq2 = zip(*sorted(enumerate(seq2), key=itemgetter(1)))
        I = np.eye(len(iseq1))
        sg1 = np.linalg.det(I[:,iseq1])
        I = np.eye(len(iseq2))
        sg2 = np.linalg.det(I[:,iseq2])
        b1 = np.zeros(2*Nst)
        b2 = np.zeros(2*Nst)
        for x in seq1:
            b1[x-1]=1
        b1u = b1[:Nst]; b1d= b1[Nst:(2*Nst)+1]
        b1 = b1u - 2*b1d; 
        b1[b1==-1] = 2;b1[b1==-2] = -1
        for x in seq2:
            b2[x-1]=1
        b2u = b2[:Nst]; b2d= b2[Nst:(2*Nst)+1]
        b2 = b2u - 2*b2d; 
        b2[b2==-1] = 2;b2[b2==-2] = -1        
        ib1 = np.where((B==b1).all(axis=1))
        ib2 = np.where((B==b2).all(axis=1))
        ib1 = np.concatenate(ib1)
        ib2 = np.concatenate(ib2)
        if len(ib1)!=0 and len(ib2)!=0:
            Rho[j-1,k-1] =  sg1*sg2*Cg[ib1]*Cg[ib2]  
            Rho[k-1,j-1] = Rho[j-1,k-1]
        return Rho

###########################
    
    Cg = C.copy()
    print('\nNATURAL ORBITALS\n')
    Rho = np.zeros([2*Nst,2*Nst])
    Rho_trash = []
    
    for j in range(1,(2*Nst)+1):
        for k in range(j,(2*Nst+1)):
            vv = [x for x in range(1,(2*Nst)+1)]; vv = [x for x in vv if x!=(j) and x!=(k)]
            combss = list()
            combss += list(combinations(vv, Nel-1))
            combss = np.array(combss)
            n = [x for x in range(len(combss))]
            pool = Pool(P)    
            Rho_trash.append(pool.map(call,n))

    prefac = 1
    Rho_trash = np.array(Rho_trash,dtype=object)
    Rho_trash = np.concatenate(Rho_trash)
    Rho = Rho_trash.sum(axis=0)
    Rho = prefac*Rho
    Enat, Cnat = sp.linalg.eigh(np.add(Rho[:Nst,:Nst],Rho[Nst:,Nst:]))


    idx = np.argsort(Enat)[::-1]

    Enat = Enat[idx];Enat = np.round(Enat,2)
    Cnat = Cnat[:,idx]
    # Enat.sort(); Enat = Enat[::-1]; Enat = np.round(Enat,2)


    print("\nOccupation of Natural Orbitals are:")
    print(Enat)

    write_molden("cas_orbital.molden.input",Cnat,'natural_orbitals.molden.input')
    et = time.time()
    Print_time_elepsed(st,et)

