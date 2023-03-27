import numpy as np
from itertools import combinations
import math



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
        return B



class CAS_CI:
    def __init__(self,CI=None,spin=None):
        self.CI=CI
        self.spin=spin

    def CAS_construct_basis(self,Nst,Nel,max_oc,Sz):
        def nchoosek(n, k):
            return math.comb(n, k)
        if not type(max_oc) is list:
            max_oc= [max_oc]
        Nstps = Nst
        Nst = 2*Nst
        Nbas = nchoosek(Nst,Nel)
        n= [j for j in range(1,(Nst+1))]
        basisel = list()
        basisel += list(combinations(n, Nel))
        B = np.zeros([Nbas,Nstps])
        ind_garb = []

        if len(max_oc)==1:
            for j in range(Nbas):
                vbasis = np.zeros([1,Nst])
                for l in basisel[j]:
                    vbasis[0,int(l)-1]=1
                vbasis = np.reshape(vbasis, (2,Nstps),order="F")
                spinz= sum(vbasis[0])- sum(vbasis[1])   
                vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1])
                if spinz==Sz and  max(vbasis) < (sum(max_oc)+1):
                    B[j] = vbasis
                else:
                    ind_garb.append(j)
            B= np.delete(B,ind_garb,0)

        else:
            Nlevel, Nelec = max_oc[0], max_oc[1]
            for j in range(Nbas):
                vbasis = np.zeros([1,Nst])
                for l in basisel[j]:
                    vbasis[0,int(l)-1]=1
                vbasis = np.reshape(vbasis, (2,Nstps),order="F")
                spinz= sum(vbasis[0])- sum(vbasis[1])   
                vbasis = vbasis[0]-vbasis[1] + np.multiply( 2*vbasis[0],vbasis[1] )
                if spinz==Sz and  sum(abs(vbasis[0:Nlevel]))==Nelec:
                    B[j] = vbasis
                else:
                    ind_garb.append(j)
            B= np.delete(B,ind_garb,0)
        return B

def FCIDUMP_Hfile(FCI_file):
    ''' Convert the give FCIDUMP file into out Hfile format from Eh to eV units.'''
    def FCIDUMP2Hfile(s,u,v,w,z):
        FCIDUMP = [s,u,v,w,z]  
        Hfile_line = list()
        if w==z and z==0:
            a = FCIDUMP
            Hfile_line += [(a[0],a[2],a[3],a[4],a[1])]
            if u ==v and v== w and w==z and u==0:
                Hfile_line = []
        else:
            a = FCIDUMP
            b = [(a[0],a[1],a[3],a[4],a[2]),(a[0],a[2],a[3],a[4],a[1]),(a[0],a[1],a[4],a[3],a[2]),(a[0],a[2],a[4],a[3],a[1])]
            c = sorted(set(b), key=b.index)
            toremove = []
            for k in c:
                allsign = [(k[0],k[1],k[2],k[3],k[4]),(k[0],k[1],-1*k[2],-1*k[3],k[4]),(k[0],-1*k[1],k[2],k[3],-1*k[4]),(k[0],-1*k[1],-1*k[2],-1*k[3],-1*k[4])]            
                for l in allsign:    
                    if l[-1]!=l[-2] and l[-3]!=l[-4] and l not in toremove:
                        Hfile_line += [l]
                    toremove += [l]
                    toremove += [(l[0],l[2],l[1],l[4],l[3])]
        return Hfile_line

    # print(FCIDUMP2Hfile(0.4,2,1,2,1))
    Hfile = list()
    with open(FCI_file) as f1:
        f1.readline();f1.readline();f1.readline();f1.readline()
        for line in f1:
            split = line.split()
            c=format(float(split[0]), '.8f')
            coeff =27.2114*float(c) ;d = int(split[1]);e=int(split[2]);f=int(split[3]);g=int(split[4])
            Hfile_line= FCIDUMP2Hfile(coeff,d,e,f,g)
            Hfile += Hfile_line
    with open('Hfile.txt','w') as f2:
        for j in Hfile:
            f2.write('%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n ' % (j[0], j[1], j[2],j[3],j[4]))    
    return np.array(Hfile)    


def clean_Hfile(Hfile):
    ''' To clean the Hfile. Remove the repeated entries and int less then thresold'''
    Hfile = np.array(Hfile)
    indsgarb = np.argwhere(np.abs(Hfile[:,0]) < 0.00000001)
    Hfile=np.delete(Hfile,indsgarb,0)
    a=0
    while a < len(Hfile):
        kk = np.where((Hfile[a,1:]==Hfile[:,1:]).all(axis=1))
        kk = np.concatenate(kk)
        if len(kk)>=2:
            coeff = np.sum(Hfile[kk,0],axis=0)
            Hfile[a,0]=coeff
            Hfile= np.delete(Hfile,kk[1:],axis=0)
        a+=1
    indsgarb1 = np.argwhere(np.abs(Hfile[:,0]) < 0.00000001)
    Hfile=np.delete(Hfile,indsgarb1,0)
    return(Hfile)

'''__________________ Tools for construction of S^2 matrix for given number of sites__________________'''
def product_Of_two_1e_Hfile(O1,O2):
    n1 = O1.shape[0]; n2 = O2.shape[0]; O3 = []
    for l1 in range(n1):
        w1 = O1[l1,0]; c1= O1[l1,1]; d1= O1[l1,4]
        for l2 in range(n2):
            w2 = O2[l2,0]; c2= O2[l2,1]; d2= O2[l2,4]
            if d1 != c2 and d1 != d2:
                O3_new =  [w1*w2 , c1 , c2 , d2 , d1]
                O3.append(O3_new)
            if d1 == c2:
                O3_new =  [w1*w2 ,c1, 0 , 0, d2]
                O3.append(O3_new)
                if d1 != d2:
                    O3_new = [-1*(w1*w2), c1, c2, d1 , d2]
                    O3.append(O3_new)
    return O3

def Pauli_matrix_Hfile(site):
    Sx = [[0.5 , site, 0, 0, -site],[0.5, -site, 0, 0, site]]
    Sy = [[-0.5*1j, site, 0, 0, -site], [0.5*1j, -site, 0, 0, site]]
    Sz = [[0.5 , site, 0, 0, site], [-0.5, -site, 0, 0, -site]]
    return(Sx,Sy,Sz)

def one_and_two_body_Sfile(Nst):
    SxT,SyT, SzT = [],[] ,[]
    for site in range(1,(Nst+1)):
        Sx1,Sy1,Sz1 = Pauli_matrix_Hfile(site)
        SxT += Sx1; SyT += Sy1; SzT += Sz1 
    SxT, SyT, SzT = np.array(SxT), np.array(SyT), np.array(SzT)
    S2file= product_Of_two_1e_Hfile(SxT,SxT)
    S2file = np.array(S2file, dtype=np.complex_)
    S2file=np.vstack((S2file,product_Of_two_1e_Hfile(SyT,SyT)))
    S2file=np.vstack((S2file,product_Of_two_1e_Hfile(SzT,SzT)))
    S2file = S2file.real
    S2file = clean_Hfile(S2file)
    S2file1body = S2file[S2file[:,2]==0]
    S2file2body = S2file[S2file[:,2]!=0]
    return S2file1body, S2file2body

'''____________________ Done for construction of S^2 file._________________________'''


