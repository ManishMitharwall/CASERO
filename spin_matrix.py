import numpy as np
import cmath
def p_matrix(s):
    s=float(s)
    n = int(2*s+1)
    sx=np.empty([n,n])
    sy=np.empty([n,n],dtype=complex)
    sz=np.empty([n,n])
    for a in range(1,n+1):
        for b in range(1,n+1):
            sx[a-1,b-1]=( 0.5*((a==b+1) + (a+1==b))*np.sqrt((s+1)*(a+b-1)-a*b))
            sy[a-1,b-1] =  1j*(0.5*((a==b+1) - (a+1==b))*np.sqrt((s+1)*(a+b-1)-a*b))
            sz[a-1,b-1] = (s+1-a)*(a==b)

    return sx,sy,sz
    

