import numpy as np
import scipy as sp
from Orbital_tools import write_molden



def Natural_Transition_orbitals(C0,B0,C1,B1,s1,s2):
    # AE=  [[ 0.1847,   0.5470,   0.1469,  -0.5584],[-0.3814,   0.4335,   0.4101,   0.4064],[-0.5661,  -0.1135,  -0.5570,   0.1520],  [-0.1847,  -0.5470,   0.1469,  -0.5584],[0.3814,  -0.4335,   0.4101,   0.4064],[0.5661,   0.1135,  -0.5570,   0.1520]]

    Nst = int(np.size(B0[0,:])) 
    Nel0 = int(np.sum(np.abs(B0[0,:])))
    s1 = 0.5*(1-s1); s2 = 0.5*(1-s2)
    
    Nelu0 = len(np.argwhere(B0[0]==1)) + len(np.argwhere(B0[0]==2))
    Neld0 = len(np.argwhere(B0[0]==-1)) + len(np.argwhere(B0[0]==2))
    Nelu1 = len(np.argwhere(B1[0]==1)) + len(np.argwhere(B1[0]==2))
    Neld1 = len(np.argwhere(B1[0]==-1)) + len(np.argwhere(B1[0]==2))
    
    B0u = np.array(B0);  B0d = np.array(B0); B0u[B0u==2]=1; B0u[B0u==-1]=0
    B0d[B0d==1]=0; B0d[B0d==-1]=1; B0d[B0d==2]=1
    B0t = np.concatenate((B0u,B0d),axis=1)
    
    B1u = np.array(B1);  B1d = np.array(B1); B1u[B1u==2]=1; B1u[B1u==-1]=0
    B1d[B1d==1]=0; B1d[B1d==-1]=1; B1d[B1d==2]=1
    B1t = np.concatenate((B1u,B1d),axis=1)
    
    tole = 0.000001
    C0_important = np.argwhere(np.abs(C0)> tole)
    F2 = np.zeros((Nst,Nst))
    F= np.zeros((Nst,1))
 
    
    for A in range(Nst):
        for B in range(Nst):
            for J in range(len(C0_important)):
                J = C0_important[J]
                b = B0t[J,:];b = np.concatenate(b)
                ib = np.argwhere(b==1)

                for l in range(Nel0):
                    if ib[l] == (B+(s2*Nst)):
                        ibp = ib
                        ibp[l] = A + s1*Nst
                        bp = np.zeros((2*Nst))
                        bp[ibp] = 1
                        indbp = np.where(([bp] == B1t).all(axis=1));indbp= np.concatenate(indbp)
                        if len(indbp)>0:
                            F2[A,B] = np.add(F2[A,B] ,C1[indbp]*C0[J])

                            
    x =F2@(F2.T)
    L,U = np.linalg.eigh(x)
    L = np.round(L,3)
    print("\nNATURAL TRANSITION ORBITALS\n",L)
    
    write_molden("cas_orbital.molden.input",U,'NTO_orbitals.molden.input')
    




    
