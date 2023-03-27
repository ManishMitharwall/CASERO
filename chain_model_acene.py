import numpy as np
import sys
from CASero import CASERO_CODE
from sparce_matrix import *
my = CASERO_CODE()

Hfile = my.FCI_2_Hf(sys.argv[1])
M = 6
n = 4
nprocs = 16
# add a chain of n sites

# Hamiltonian of the chain
T = 1  # this is the hopping inside the chain
e = 5.0  # this is the energy of the metal levels: it should be "in gap".

# For that we can do the following: Take the lines in the Hfile of the form
# [e1, a, 0, 0, a] and [e2, b, 0, 0, b] and define e = 0.5*(e1+e2),
# where a and b are the indices of the HOMO and the LUMO, respectively.
# in the case of Huckel_Hubbard this very easy... It's just 0!
# In the ab initio case, I think we ought to do this average.

for j in range(M+1, M+n):
    Hfile = np.vstack((Hfile, [T, j, 0, 0, j+1]))
    Hfile = np.vstack((Hfile, [e, j, 0, 0, j]))

Hfile = np.vstack((Hfile, [e, M+n, 0, 0, M+n]))

# Connection chain-molecule
# let's say that M+1 is the level that connects with the molecule.
# We have to define hoppings with the molecular orbitals!
# This is probably the trickiest part of the model.
# I made up a simple model: the largest hopping is with HOMO / LUMO
# and then the coupling decreases linearly.
# More exactly, we should do this as a function of the orbital energy.
# For some cases we can try to change just by hand.
t_max = 0.1  # Hopping of the chain with the HOMO/LUMO.
# This is the parameter whose impact we should explore.
t = np.hstack((np.linspace(0, t_max, int(M/2)), np.flip(np.linspace(0, t_max, int(M/2)))))

for j in range(1,M+1):
    Hfile = np.vstack((Hfile, [t[j-1], j, 0, 0, M+1]))

with open('newHfile.txt','w') as f2:
    for j in Hfile:
        f2.write('%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n ' % (j[0], j[1], j[2],j[3],j[4]))    

exit()


# With this everything is ready for the calculation...!
# In order to be able to make an analysis of it,
# we need to construct the S2 operator of the molecule alone.
# If t is small, this will still give a good quantum number.
# I do it here by imitating what we already have but truncating to levels 1,...,M:

# Generate Hfile of S^2 operator

B = my.CI_cas(M+n, M+n, 2, 0, 0)
S2file1body, S2file2body = my.Sfile_12(M)
S2mol = create_hamiltonian_sparse_matrix_general_nonreduced(S2file1body,S2file2body,B,[0],[0],0,M+n,Processors=nprocs)



