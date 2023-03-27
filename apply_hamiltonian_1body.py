import numpy as np
from apply_operator import apply_operator_1body


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


# print(apply_hamiltonian_1body([2,-1,1], [-3,1,0,0,2], 2) )