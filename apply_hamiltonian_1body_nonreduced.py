import numpy as np
from apply_operator import apply_operator_1body


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


# print(apply_hamiltonian_1body([2,-1,1], [-3,1,0,0,2], 2) )
