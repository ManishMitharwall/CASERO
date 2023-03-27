from NTO_orbitals import Natural_Transition_orbitals
from Dyson_orbitals import Dyson_Orbitals
from natural_orbitals import natural_Orbitals

class DSPMM_CODE:

    def Dysons(self,C,Cm,CM,B,Bm,BM,nstat):
        Dyson_Orbitals(C,Cm,CM,B,Bm,BM,nstat)