import numpy as np
import pyscf
from pyscf.tools import fcidump
from pyscf import fci
from NTO_orbitals import Natural_Transition_orbitals
from Dyson_orbitals import Dyson_Orbitals
from natural_orbitals import natural_Orbitals
from construct_basis import construct_basis
from clean_Hfile import clean_Hfile
from natural_orbitals import natural_Orbitals
from Orbital_tools import write_molden
from Orbital_tools import getMOLDENcoefs
from print_level import *
import time

class DSPMM_CODE:
    def __init__(self) -> None:
        Print_Title()
    def Run_Dysons(self,C,Cm,CM,B,Bm,BM,nstat):
        Dyson_Orbitals(C,Cm,CM,B,Bm,BM,nstat)
    def Run_NTOs(self,C0,B0,C1,B1,s1,s2):
        Natural_Transition_orbitals(C0,B0,C1,B1,s1,s2)
    def construct_basis(self,Nst,Nel,max_oc,Sz,ifBM_Bm):
        return construct_basis(Nst,Nel,max_oc,Sz,ifBM_Bm)
    def clean_Hfile(self,Hfile):
        return clean_Hfile(Hfile)
    def Run_Natural_orb(self,B,C,Nst,Nel,P):
        natural_Orbitals(B,C,Nst,Nel,P)

    '-----------------------Now for pyscf use-----------------------------------------------------------'
    def Int_z_FCIDUMP(self, FCIDUMP):
        result = fcidump.read(FCIDUMP)
        h1e = result['H1']
        h2e = result['H2']
        norb = result['NORB']
        nelec = result['NELEC']
        ecore = result['ECORE']
        Sz = result['MS2']
        return h1e, h2e, norb, nelec, ecore, Sz

    def Run_FCI(self,h1e, h2e, norb, nelec, nroots=8, ecore=0 ,max_space=30, max_cycle=100, verbos = True):
        st = time.time()
        Print_input_parameters(norb,nelec,nelec,2,'YES')
        Energies, fcivec = fci.direct_spin1.kernel(h1e, h2e, norb, nelec, nroots=nroots, ecore=ecore ,max_space=max_space, max_cycle=max_cycle)
        if verbos:
            print('\nSTATE \t\t ENERGY(eV) \t\t Sz \t\t S^2')
            for i, c in enumerate(fcivec):
                S2 = pyscf.fci.spin_op.spin_square(c, norb, nelec)[0]
                Szz = 0.5*(-1+np.sqrt(1+4*S2))
                print('Ψ = %d \t\t  %.4f \t\t  %.2f \t\t %.2f' % (i, Energies[i]*27.2114, Szz, S2))
            print(f'Energy gap of ground state and excited state is  {np.round((Energies[1]-Energies[0])*27.2114,5)}')
            self.CI_large(fcivec[0], norb, nelec, weight=0.08, state=0)
            self.CI_large(fcivec[1], norb, nelec, weight=0.08, state=1)
        et = time.time();Print_time_elepsed(st,et)
        return Energies, fcivec

    def Run_FCI_S0(self,h1e, h2e, norb, nelec, nroots=8, ecore=0 ,max_space=30, max_cycle=100, verbos = True):
        st = time.time()
        Print_input_parameters(norb,nelec,nelec,2,'YES')
        Energies, fcivec = fci.direct_spin0.kernel(h1e, h2e, norb, nelec, nroots=nroots, ecore=ecore ,max_space=max_space, max_cycle=max_cycle)
        if verbos:
            print('\nSTATE \t\t ENERGY(eV) \t\t Sz  \t\t S^2')
            for i, c in enumerate(fcivec):
                S2 = pyscf.fci.spin_op.spin_square(c, norb, nelec)[0]
                Szz = 0.5*(-1+np.sqrt(1+4*S2))
                print('Ψ = %d \t\t  %.4f \t\t  %.2f \t\t %.2f' % (i, Energies[i]*27.2114, Szz, S2))
            self.CI_large(fcivec[0], norb, nelec, weight=0.08, state=0)
            self.CI_large(fcivec[1], norb, nelec, weight=0.08, state=1)
            #self.CI_large(fcivec[2], norb, nelec, weight=0.08, state=2)
            #self.CI_large(fcivec[3], norb, nelec, weight=0.08, state=3)
            #self.CI_large(fcivec[4], norb, nelec, weight=0.08, state=4)
        et = time.time();Print_time_elepsed(st,et)
        return Energies, fcivec        

    def py_basis(self,C0, norb, nelec):
        fcivec = fci.addons.large_ci(C0,norb, nelec, tol= getattr(pyscf.__config__, 'fci_addons_large_ci_tol', 0.0) ,return_strs=False)
        coff_arr, alpha, beta = [], [], []
        for i in range(len(fcivec)):
            coff_arr.append(fcivec[i][0])
            alpha.append(fcivec[i][1])
            beta.append(fcivec[i][2])
        abasis = np.zeros((len(alpha),norb))
        bbasis = np.zeros((len(alpha),norb))
        for i in range(len(alpha)):
            abasis[i][alpha[i]] = 3
            bbasis[i][beta[i]] = -1
        basis = abasis + bbasis
        basis[basis ==3 ] = 1
        coff_arr = np.array(coff_arr)
        return basis.astype(int)
    
    def CI_large(self,C0, norb, nelec, weight=0.08, state =0):
        fcivec = fci.addons.large_ci(C0,norb, nelec, tol= getattr(pyscf.__config__, 'fci_addons_large_ci_tol', weight) ,return_strs=False)
        coff_arr, alpha, beta = [], [], []
        for i in range(len(fcivec)):
            coff_arr.append(fcivec[i][0])
            alpha.append(fcivec[i][1])
            beta.append(fcivec[i][2])
        abasis = np.zeros((len(alpha),norb))
        bbasis = np.zeros((len(alpha),norb))
        for i in range(len(alpha)):
            abasis[i][alpha[i]] = 3
            bbasis[i][beta[i]] = -1
        basis = abasis + bbasis
        basis[basis ==3 ] = 1
        coff_arr = np.array(coff_arr)
        idx = np.argsort(np.abs(coff_arr))[::-1]
        basis = basis[idx];coff_arr = coff_arr[idx]; basis = basis.astype(int)
        print(f'\nWeight of slater determinant of Ψ = {state} \nWEIGHT \t\t COEFF \t\t\t CI')
        for i in range(len(basis)):
            print(f'{str(round( coff_arr[i]**2,3)):<5}{str(round(coff_arr[i],3)):>15}{str(basis[i]):>30}')
            # print(f'{round( coff_arr[i]**2,3)} \t {round(coff_arr[i],3)} \t {basis[i]} ')
    
    def pyscf_no(self,C0, norb, nelec ):
        p1 = fci.direct_spin1.make_rdm1(C0,norb, nelec)
        n_nat,c_nat = np.linalg.eigh(p1)
        idx = np.argsort(n_nat)[::-1]
        n_nat = n_nat[idx];n_nat = np.round(n_nat,2)
        c_nat = c_nat[:,idx]
        print("\nOccupation of Natural Orbitals are:"); print(n_nat)
        write_molden("cas_orbital.molden.input",c_nat,'natural_orbitals.molden.input')

    def py_trdm1(self, C0, n, m, norb, nelec):
        trdm = fci.direct_spin1.trans_rdm1(C0[n],C0[m], norb, nelec)
        atomic_basis, Energies = getMOLDENcoefs("cas_orbital.molden.input")
        f_out_file = open(f'tdm_{n}_{m}.molden.input','w')
        a_len=len(atomic_basis)
        n_rho = (atomic_basis@trdm@atomic_basis.T)
        n_rho = np.diagonal(n_rho);print(f'Norm Of tdm {n} to {m} is {np.linalg.norm(n_rho)}')
        with open("cas_orbital.molden.input",'r') as mfile:
            for line in mfile:
                f_out_file.write(line)
                if r'MO' in line: break
        counter= -0.1500
        for p in range(1):
            f_out_file.write(' Sym=      1a\n Ene= '+str(counter)+'\n'+ ' Spin= Alpha\n Occup= 1.000000 \n')
            counter+=0.010
            for t in range(int(a_len)):
                f_out_file.write(f' {t+1} \t {n_rho[t]}\n')
        f_out_file.close()



    def _unpack_NTO_el(self, nel):
        if (nel % 2) == 0:
            a_el = int((nel // 2)+1)
            b_el = int((nel // 2)-1) 
        else:
            a_el = int((nel // 2)+2)
            b_el = int((nel // 2)-1)
        return (a_el, b_el ) 

    def py_NTO(self, C0, h1e, h2e, norb, nelec, ecore = 0):
        B0 = self.py_basis( C0[0], norb, nelec)
        nto_elec = self._unpack_NTO_el(nelec)
        E1,C1 = self.Run_FCI( h1e, h2e, norb, nto_elec, ecore = ecore , verbos= False)
        B1 = self.py_basis( C1[0], norb, nto_elec)
        C0_gd = np.concatenate(C0[0])
        C1_gd = np.concatenate(C1[0])
        self.Run_NTOs(C0_gd,B0,C1_gd,B1,1,-1)

    def py_Dysons(self, C0, h1e, h2e, norb, nelec, ecore = 0, nstat = 1 ):
        Em,Cm = self.Run_FCI( h1e, h2e, norb, int(nelec-1), ecore = ecore )
        EM,CM = self.Run_FCI( h1e, h2e, norb, int(nelec+1), ecore = ecore )
        B0 = self.py_basis( C0[0], norb, nelec)
        BM = self.py_basis( CM[0], norb, int(nelec+1))
        Bm = self.py_basis( Cm[0], norb, int(nelec-1))
        for i in range(nstat):
            C0_gd = np.concatenate(C0[0])
            Cm_gd = np.concatenate(Cm[i])
            CM_gd = np.concatenate(CM[i])
            self.Run_Dysons(C0_gd,Cm_gd,CM_gd,B0,Bm,BM,i)
        



            
    
