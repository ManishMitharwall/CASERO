import numpy as np
# from apply_operator import apply_operator_2body


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


# print(apply_hamiltonian_2body([1,2,0], np.array([[0.5,1,-1,-1,1],[0.5,2,-2,-2,2]]),2))
