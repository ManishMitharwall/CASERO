#!/auto/vestec1-elixir/home/manishkumar/.conda/envs/kpython310/bin/python3.1
import sys
sys.path.append('/storage/praha1/home/manishkumar/mybin/module/DSPMM')
from DSPMM import DSPMM_CODE
my = DSPMM_CODE()


h1e, h2e, norb, nelec, ecore, Sz = my.Int_z_FCIDUMP( sys.argv[1])    # DOING  Many body calculation for neutral system
e0,C0 = my.Run_FCI( h1e, h2e, norb, nelec, ecore = ecore)          # DOING  Many body calculation for neutral system
my.py_trdm1(C0,0,1,norb=norb,nelec=nelec)
my.py_trdm1(C0,0,2,norb=norb,nelec=nelec)
my.py_trdm1(C0,0,3,norb=norb,nelec=nelec)
my.pyscf_no(C0[0], norb, nelec )                                   # Doing Natural orbitals 
# my.py_NTO( C0, h1e, h2e, norb, nelec, ecore = ecore)               # Doing NTO
# my.py_Dysons( C0, h1e, h2e, norb, nelec, ecore = ecore, nstat = 4 )    # Doing Dysons

