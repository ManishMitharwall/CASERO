import numpy as np
from itertools import combinations

class RAS_CI:
    def __init__(self,CI=None,spin=None):
        self.CI=CI
        self.spin=spin

    def B_RAS(self,Nst,Nel,n_RAS,excitation):
        Nstps = Nst
        Nst = 2*Nst   
        n= [j for j in range(1,(Nst+1))]
        basisel = list()

        if n_RAS==1:
            for i in range(0,excitation+1):
                basisel += list(combinations(n, Nel+i))
        if n_RAS==2:
            for i in range(-1*excitation,excitation+1):
                basisel += list(combinations(n, Nel+i))
        if n_RAS==3:
            for i in range(-1*excitation,1):
                basisel += list(combinations(n, Nel+i))    
        
        B = [] 
        for j in range(len(basisel)):
            vbasis = np.zeros([1,Nst])
            for l in basisel[j]:
                vbasis[0,int(l)-1]=1
            vbasis = np.reshape(vbasis, (2,Nstps),order="F")
            spinz= sum(vbasis[0])- sum(vbasis[1])   
            vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1])
            B.append(RAS_CI(vbasis,spinz))
        return B
    
    def cons_R_BASIS(self,CAS_el,CAS_orb,RAS1_orb,RAS3_orb,Sz,excitation=1):
        Nel_RAS2,Nst_RAS2=CAS_el,CAS_orb
        Nel_RAS1,Nst_RAS1=(2*RAS1_orb-excitation),RAS1_orb
        Nel_RAS3,Nst_RAS3=(excitation),RAS3_orb
        Nel=Nel_RAS1+Nel_RAS2+Nel_RAS3

        RAS1 = self.B_RAS(Nst_RAS1,Nel_RAS1,n_RAS=1,excitation=excitation)
        RAS2= self.B_RAS(Nst_RAS2,Nel_RAS2,n_RAS=2,excitation=excitation)
        RAS3= self.B_RAS(Nst_RAS3,Nel_RAS3,n_RAS=3,excitation=excitation)

        RAS_AC=[]
        com_RAS=[]
        for i in RAS1:
            for j in RAS3:
                RAS_AC.append(RAS_CI([list(i.CI)]+[list(j.CI)],i.spin+j.spin))
        for k in RAS_AC:
            for l in RAS2:
                if ((k.spin+l.spin)==Sz) and (np.sum(np.abs(k.CI))+np.sum(np.abs(l.CI)))==Nel:
                    com_RAS.append(RAS_CI(k.CI[0]+list(l.CI)+k.CI[1],k.spin+l.spin))
        return com_RAS

    def RAS_construct_basis(self,CAS_el,CAS_orb,RAS1_orb,RAS3_orb,Sz,excitation=1):
        BAS=RAS_CI().cons_R_BASIS(CAS_el,CAS_orb,RAS1_orb,RAS3_orb,Sz,excitation=1)
        B=[basis.CI for basis in BAS]
        B=np.array(B)
        BM,Bm = np.zeros([2,2,2]), np.zeros([2,2,2])
        return B,BM,Bm




a=RAS_CI().RAS_construct_basis(4,4,3,3,0,excitation=1)

print(len(a[0]))



