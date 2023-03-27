from spin_matrix import p_matrix
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def interaction_two_spins(s1,s2,J,D1,D2):
    def eigen(A):
        eigenValues, eigenVectors = sp.linalg.eig(A)
        idx = np.argsort(eigenValues)
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:,idx]
        eigenValues=[np.real(x) for x in eigenValues]
        return eigenValues, eigenVectors

    def solves(H,S):
        eigenValues,eigenVectors=eigen(H)
        solv=np.linalg.solve(eigenVectors,(np.dot(S,eigenVectors)))
        svec=np.diag(solv)
        sval =  [0.5*(-1+np.sqrt(1+4*i)) for i in  svec]
        svec= [np.real(x) for x in svec]
        sval= [np.real(x) for x in sval]
        return  svec, sval
    jlist = J 
    if not type(jlist) is list:
        jlist= [J]
    amatrix=p_matrix(s1)
    Sx1,Sy1,Sz1=amatrix[0],amatrix[1],amatrix[2]
    bmatrix=p_matrix(s2)
    Sx2,Sy2,Sz2=bmatrix[0],bmatrix[1],bmatrix[2]

    n1=len(Sx1)
    n2=len(Sx2)
    Sx1, Sy1, Sz1 = np.kron(Sx1,np.eye(n2)), np.kron(Sy1,np.eye(n2)), np.kron(Sz1,np.eye(n2))
    Sx2, Sy2, Sz2 = np.kron(np.eye(n1),Sx2), np.kron(np.eye(n1),Sy2), np.kron(np.eye(n1),Sz2)

    Sxt,Syt,Szt = Sx1+Sx2, Sy1+Sy2,Sz1+Sz2
    S2, S21, S22 = np.dot(Sxt,Sxt)+np.dot(Syt,Syt)+np.dot(Szt,Szt) ,np.dot(Sx1,Sx1)+np.dot(Sy1,Sy1)+np.dot(Sz1,Sz1),np.dot(Sx2,Sx2)+np.dot(Sy2,Sy2)+np.dot(Sz2,Sz2)
    
    H0 = D1[0]*np.dot(Sx1,Sx1) + D1[1]*np.dot(Sy1,Sy1) + D1[2]*np.dot(Sz1,Sz1) + D2[0]*np.dot(Sx2,Sx2) + D2[1]*np.dot(Sy2,Sy2) + D2[2]*np.dot(Sz2,Sz2) 


    for J in jlist:
        H = H0 + J*(np.dot(Sx1,Sx2) + np.dot(Sy1,Sy2) + np.dot(Sz1,Sz2))
        eigenValues=eigen(H)[0]

        eigS2,sval=solves(H,S2)
        eigS21, sval1=solves(H,S21)
        eigS22, sval2= solves(H,S22)
        eigSz1,svalZ1= solves(H,Sz1)[0],solves(H,Sz1)[0]
        print(f' For value of J {J}')
        for i in range(len(sval)):   
            print('%.6f \t %.6f \t %.6f \t %.6f \t %.6f \n ' % (eigenValues[i], sval[i], sval1[i], sval2[i], svalZ1[i]))
            x= eigenValues[i]-eigenValues[0]
            plt.scatter(J,x, c='red')

    plt.xlabel("J (meV)",size=14)
    plt.ylabel("Energy (meV)",size=14)
    plt.show()

interaction_two_spins(1.0,0.5,[0,1,2],[0,0,5],[1,1,1])
