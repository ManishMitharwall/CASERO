import numpy as np

def apply_operator_1body(i,l,V):
    signo = 1
    # destroy l
    signo = signo*((-1)**(sum(V[:l-1])))
    V[l-1] = 0*(V[l-1]==1) + 2*(V[l-1] != 1)
    # create i
    signo = signo*((-1)**(sum(V[:i-1])))
    V[i-1]= 1*(V[i-1]==0) + 2*(V[i-1] !=0)

    if 2 in V or -2 in V:
        V = []

    return V, signo


def apply_operator_2body(i,j,k,l,V):
    signo = 1
    # destroy l
    signo = signo*((-1)**(sum(V[:l-1])))
    V[l-1] = 0*(V[l-1]==1) + 2*(V[l-1] != 1)
    # destro k
    signo = signo*((-1)**(sum(V[:k-1])))
    V[k-1] = 0*(V[k-1]==1) + 2*(V[k-1] != 1)
    # create j
    signo = signo*((-1)**(sum(V[:j-1])))
    V[j-1]= 1*(V[j-1]==0) + 2*(V[j-1] !=0)
    # create i
    signo = signo*((-1)**(sum(V[:i-1])))
    V[i-1]= 1*(V[i-1]==0) + 2*(V[i-1] !=0)

    if 2 in V or -2 in V:
        V = []
    return V, signo