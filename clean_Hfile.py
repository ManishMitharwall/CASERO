import numpy as np

def clean_Hfile(Hfile):
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

 









# print(clean_Hfile(Hfile))

