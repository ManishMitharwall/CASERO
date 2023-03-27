#!/auto/vestec1-elixir/home/manishkumar/.conda/envs/kpython310/bin/python3.1
import sys
sys.path.append('/storage/praha1/home/manishkumar/mybin/module/DSPMM')

import argparse
import numpy as np
from create_H_S_array import *
from Many_body_hamil_S2 import Many_body_hamil_S2
from NTO_orbitals import Natural_Transition_orbitals
from Dyson_orbitals import Dyson_Orbitals
from natural_orbitals import natural_Orbitals
from print_level import *

def run():
    Hfile_file=sys.argv[1]
    Hfile=readHfile(Hfile_file)
    # Hfile = np.array([[-3,1,0,0,2],[4,1,-1,-1,1],[4,2,-2,-2,2]])
    parser = argparse.ArgumentParser(description='Perform Many body calculations')
    parser.add_argument('FCIDUMP_FILE',metavar='FCIDUMP_FILE',help='FCIDUMP_FILE_PLEASE')
    parser.add_argument('-p','-P',metavar='NProcs', default=1,help='Number of processor used')
    parser.add_argument('-Nst','-NST','-nst',metavar='No_of_sites', default= int(np.max(Hfile[:,1:])), help='Number of sites')
    parser.add_argument('-Nel','-nel','-NEL', metavar='No_of_electron', default=int(np.max(Hfile[:,1:])), help='Number of electron')
    parser.add_argument('-Sz', '-sz','-SZ', metavar='SZ',default=0, help='Spin of system (Sz)')
    parser.add_argument('-mocc', metavar='mac_occ', default=2, help='max occupation at site')
    args = parser.parse_args()
    Nst=int(args.Nst)
    Nel=int(args.Nel)
    P = int(args.p)
    max_oc=int(args.mocc)
    Sz=int(args.Sz)
    if Sz== 0:
        SzM=1
        Szm=1
    if Sz==1:
        SzM=0
        Szm=0

    Print_Title()
    cut_lumo = False
    if cut_lumo:
        c_value = -1.0        #Give  -Ve value
        print(f"Changing the Unocc level energies by {c_value}\n")
        Hfile=scissoring_opertor(Hfile,Nel,Nst,Sz,c_value)
        
    

    # Hfile = np.array([[-3,1,0,0,2],[4,1,-1,-1,1],[4,2,-2,-2,2]])
    B0,C0=Many_body_hamil_S2(Hfile,Nel,Nst,Sz,max_oc,P,S2_spin=1)
    natural_Orbitals(B0,C0[:,0],Nst,Nel,P)
    B1,C1=Many_body_hamil_S2(Hfile,Nel,Nst,Sz+2,max_oc,P,S2_spin=0)
    Natural_Transition_orbitals(C0[:,0],B0,C1[:,0],B1,1,-1)
    print("\nDYSONS ORBITALS")
    Bm,Cm = Many_body_hamil_S2(Hfile,Nel-1,Nst,Szm,max_oc,P,S2_spin=1)
    BM,CM = Many_body_hamil_S2(Hfile,Nel+1,Nst,SzM,max_oc,P,S2_spin=1)
    for n in range(4):
        C0_gd = C0[:,0]; Cm_gd = Cm[:,n]; CM_gd = CM[:,n]
        Dyson_Orbitals(C0_gd,Cm_gd,CM_gd,B0,Bm,BM,nstat=n)
    print("\n\t\t\tHAPPY LANDING!")
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")


if __name__== '__main__':
    run()
