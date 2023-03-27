import numpy as np
from apply_operator import apply_operator_1body
from  apply_operator import apply_operator_2body

def apply_hamiltonian_general_openshell(Ri,Hfile):
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

    V = np.zeros([1,2*Nst])
    Nst = np.max(Hfile[:,1:])
    Nst = np.abs(int(Nst))

    R = []
    for a in range(nrow_h):
  

        
        for el in range(nrow_R):
            b0 = Ri[el,1:]
            b = b0
            i, j, k, l = int(abs(Hfile[a,1])), int(abs(Hfile[a,2])), int(abs(Hfile[a,3])), int(abs(Hfile[a,4]))

            if j !=0:
                si, sj, sk, sl = np.sign(Hfile[a,1]), np.sign(Hfile[a,2]), np.sign(Hfile[a,3]), np.sign(Hfile[a,4])
                if si ==-1:
                    i = i + Nst
                if sj == -1:
                    j = j + Nst
                if sk == -1:
                    k = k + Nst
                if sl == -1:
                    l = l + Nst

                V1 = b
                V1 = [0 if x ==-1 else x for x in V1]
                V1 = [1 if x == 2 else x for x in V1]

                V2 = b
                V2 = [0 if x ==1 else x for x in V2]
                V2 = [1 if x == 2 else x for x in V2]
                V2 = [abs(x) for x in V2]
                V = V1 + V2
                V, signo = apply_operator_2body(i,j,k,l,V)

                if len(V) !=0:
                    V1,V2 =  V[0:Nst], np.multiply(0.5,V[Nst:])
                    V= np.add(V1,V2)
                    V = [2 if x==1.5 else x for x in V]
                    V = [-1 if x==0.5 else x for x in V]
                    R_new = [Ri[el,0]*Hfile[a,0]*signo] + V
                    R.append(R_new)

            else:
                si = np.sign(Hfile[a,1])
                if si == -1:
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


