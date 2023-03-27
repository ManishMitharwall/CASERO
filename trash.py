import numpy as np
import sys


fci = np.zeros((10,10,10,10))

head = ''''''
with open(sys.argv[1]) as f1:
    for i in range(4):
        head  += f1.readline()
    for line in f1:
        split = line.split()
        c=format(float(split[0]), '.8f')
        coeff =27.2114*float(c) 
        d = int(split[1]);e=int(split[2]);f=int(split[3]);g=int(split[4])
        if d==e and f==g and d==0:
            ncore = '{:.12e}' .format(float(split[0])) 
            break
        fci[d-1,e-1,f-1,g-1] =  float(split[0])


M = 6
n = 4
T = 1  # this is the hopping inside the chain
e = (fci[2,2,0,0]+ fci[3,3,0,0] )/2 # this is the energy of the metal levels: it should be "in gap".
T = T / 27.2114

for j in range(M, M+n-1):
    fci[j+1,j , 0, 0] += T
    fci[j,j , 0, 0] += e
fci[M+n-1,M+n-1 , 0, 0] += e


t_max = 0.1/27.2114  # Hopping of the chain with the HOMO/LUMO.
t = np.hstack((np.linspace(0, t_max, int(M/2)), np.flip(np.linspace(0, t_max, int(M/2)))))

for j in range(0,M):
    fci[M,j,0,0] += t[j]



fo = open('chechfci','w')

fo.write(head)
ij = 0
for i in range(len(fci)):
    for j in range(0, i+1):
        kl = 0
        for k in range(0, i+1):
            for l in range(0, k+1):
                if ij >= kl:
                    coef = '{:.12e}' .format(fci[i,j,k,l])
                    fo.write('{:<30} {:>5} {:>5} {:>5} {:>5} \n'.format(coef, i+1, j+1, k+1, l+1 ))
                kl += 1
        ij += 1

for i in range(len(fci)):
    for j in range(0, i+1):
        coef = '{:.12e}' .format(fci[i,j,0,0])
        fo.write('{:<30} {:>5} {:>5} {:>5} {:>5} \n'.format(coef, i+1, j+1, 0, 0))
fo.write('{:<30} {:>5} {:>5} {:>5} {:>5}' .format(ncore, 0, 0, 0, 0))
            
