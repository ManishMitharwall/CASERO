#!/auto/vestec1-elixir/home/manishkumar/.conda/envs/kpython310/bin/python3.1
import numpy as np
from CASIRO import Many_body_hamil_S2_DOS
from spect_fun import spectral_function
import matplotlib.pyplot as plt 
import time, os, re
os.environ['OMP_NUM_THREADS'] = "4"

Hfile = np.loadtxt('Hfile.txt')

N = 7
T = 0.5
t = 0.2
e = -2.33
Nel = N + 2
Nst = N + 2


# Hamiltonian of the chain
for j in range(3, N+2):
    Hfile = np.vstack((Hfile, [T, j, 0, 0, j+1]))
    Hfile = np.vstack((Hfile, [e, j, 0, 0, j]))

Hfile = np.vstack((Hfile, [e, Nst, 0, 0, Nst]))

# Hamiltonian of the coupling
Hfile = np.vstack((Hfile, [t, 1, 0, 0, 3]))
Hfile = np.vstack((Hfile, [t, 2, 0, 0, 3])) # for now we have the same coupling with both levels of the impurity


eta = 0.05
# Szpossible = list(np.arange(-Nst,Nst+1))
Szpossible = [-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7]

B0, C0, E0 = Many_body_hamil_S2_DOS(Hfile, Nel, Nst, Szpossible, 2, 10)
BM, CM, EM = Many_body_hamil_S2_DOS(Hfile, Nel+1, Nst, Szpossible, 2, 10)
Bm, Cm, Em = Many_body_hamil_S2_DOS(Hfile, Nel-1, Nst, Szpossible, 2, 10)



w = np.linspace(-8, 4, 1000)
sitesmol=[1,-1,2,-2]
siteschain = [3,-3,4,-4,5,-5,6,-6,7,-7,8,-8,9,-9]
sitesall = [sitesmol, siteschain]
DOSmol = spectral_function(C0,E0,B0,CM,EM,BM,Cm,Em,Bm,sitesmol,eta,w,0,2,2000)
DOSchain = spectral_function(C0,E0,B0,CM,EM,BM,Cm,Em,Bm,siteschain,eta,w,0,2,2000)
np.savetxt('moldos.txt', DOSmol)
np.savetxt('chaindos.txt',DOSchain)


plt.plot(w,DOSmol)
plt.plot(w,DOSchain)
plt.savefig('dosplot.png',dpi=300)




