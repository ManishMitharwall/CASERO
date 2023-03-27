import numpy as np
import re

def getMOLDENcoefs(fname):
    coeff = []
    with open(fname,'r') as mfile:
        while r"MO" not in mfile.readline():continue
        Energies =[];co_line = []
        for line in mfile:
            if re.search(r'Ene',line):
                split = line.split('=');
                if co_line: coeff.append(co_line); co_line = []
                Energies.append(float(split[-1]))
            if '=' not in line:
                    line = ' '.join(line.split()); split = line.split(' ')
                    co_line.append(float(split[-1]))
            if r'[' in line: break
        coeff.append(co_line)
        coeff = np.array(coeff).T
    return coeff, Energies



def write_molden(m_file,T_mat,wfile):
    atomic_orb , energies = getMOLDENcoefs(m_file)
    new_basis = atomic_orb @ T_mat
    f_out_file = open(wfile,'w')
    with open(m_file,'r') as mfile:
        for line in mfile:
            f_out_file.write(line)
            if r'MO' in line: break  

    counter= -0.1500
    for p in range(len(T_mat)):
        f_out_file.write(' Sym=      1a\n Ene= '+str(counter)+'\n'+ ' Spin= Alpha\n Occup= 1.000000 \n')
        counter+=0.010
        for t in range(len(new_basis)):
            f_out_file.write(f' {t+1}     {new_basis[t,p]}\n')
    f_out_file.close()
     

