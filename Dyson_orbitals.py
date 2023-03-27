import numpy as np
from itertools import combinations
from operator import itemgetter
from scipy.linalg import block_diag
from Orbital_tools import getMOLDENcoefs
from multiprocess import Pool, cpu_count

def Dyson_Orbitals(C,Cm,CM,B,Bm,BM,nstat):
    global Fm_process_chunk, FM_process_chunk
    Nst = int(np.size(B[0,:])) 
    Nel = int(np.sum(np.abs(B[0,:])))
    Fm = np.zeros((2*Nst,1))
    FM = np.zeros((2*Nst,1))
    print('\nDoing Dyson orbitals')


    for j in range(1,(2*Nst)+1):
        vv = [x for x in range(1,(2*Nst)+1)]; vv = [x for x in vv if x!=(j)]
        combss = list()
        combss += list(combinations(vv, Nel-1))
        combss = np.array(combss)
        for se in range(len(combss)):
            seq1 = np.insert(combss[se,:],0,j); seq2 = combss[se,:]
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
            ib2 = np.where((Bm==b2).all(axis=1))
            ib1 = np.concatenate(ib1)
            ib2 = np.concatenate(ib2)
            if len(ib1)!=0 and len(ib2)!=0:
                Fm[j-1]=Fm[j-1]+ sg1*sg2*C[ib1]*Cm[ib2]


    for j in range(1,(2*Nst)+1):
        vv = [x for x in range(1,(2*Nst)+1)]; vv = [x for x in vv if x!=(j)]
        combss = list()
        combss += list(combinations(vv, Nel))
        combss = np.array(combss)
        for se in range(len(combss)):
            seq1 = np.insert(combss[se,:],0,j); seq2 = combss[se,:]
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
            ib1 = np.where((BM==b1).all(axis=1))
            ib2 = np.where((B==b2).all(axis=1))
            ib1 = np.concatenate(ib1)
            ib2 = np.concatenate(ib2)
            if len(ib1)!=0 and len(ib2)!=0:
                FM[j-1]=FM[j-1]+ sg1*sg2*CM[ib1]*C[ib2]

#######################################################################

    atomic_basis, Energies = getMOLDENcoefs("cas_orbital.molden.input")
    f_out_file = open(f'dysons_{nstat}_orbitals.molden.input','w')

    a_len=len(atomic_basis)
    FM = block_diag(atomic_basis,atomic_basis)@FM
    Fm = block_diag(atomic_basis,atomic_basis)@Fm
    FM1 = FM[:a_len]; FM2= FM[a_len:]
    Fm1 = Fm[:a_len]; Fm2= Fm[a_len:]
    if np.sum(np.abs(FM1)) > np.sum(np.abs(FM2)):
        FM = FM1
    else:
        FM = FM2

    if np.sum(np.abs(Fm1)) > np.sum(np.abs(Fm2)):
        Fm = Fm1
    else:
        Fm = Fm2
    
    print(f"\nNorm of -1e and +1e For n = {nstat}")
    print('\t',np.round(np.linalg.norm(Fm),3),'\t',np.round(np.linalg.norm(FM),3))

    FM = FM / np.linalg.norm(FM)
    Fm = Fm / np.linalg.norm(Fm)
    FM = FM.T
    Fm = Fm.T
    F = np.concatenate((Fm,FM))

    with open("cas_orbital.molden.input",'r') as mfile:
        for line in mfile:
            f_out_file.write(line)
            if r'MO' in line: break   
    counter= -0.1500
    for p in range(len(F)):
        f_out_file.write(' Sym=      1a\n Ene= '+str(counter)+'\n'+ ' Spin= Alpha\n Occup= 1.000000 \n')
        counter+=0.010
        for t in range(int(a_len)):
            f_out_file.write(f' {t+1} \t {F[p,t]}\n')
    f_out_file.close()



