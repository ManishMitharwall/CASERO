import numpy as np
import sys
from print_level import *
from tools_dspmm import CAS_CI, RAS_CI, FCIDUMP_Hfile, one_and_two_body_Sfile

class CASERO_CODE:
    def __init__(self) -> None:
        Print_Title()
    def CI_cas(self,Nst,Nel,max_oc,Sz):
        return CAS_CI().CAS_construct_basis(Nst,Nel,max_oc,Sz)
    def CI_ras(self,CAS_el,CAS_orb,RAS1_orb,RAS3_orb,Sz,excitation=1):
        return RAS_CI().RAS_construct_basis(CAS_el,CAS_orb,RAS1_orb,RAS3_orb,Sz,excitation=excitation)
    def FCI_2_Hf(self,FCI_file):
        return FCIDUMP_Hfile(FCI_file)

    def Sfile_12(self,Nst):
        return one_and_two_body_Sfile(Nst)
    
    


