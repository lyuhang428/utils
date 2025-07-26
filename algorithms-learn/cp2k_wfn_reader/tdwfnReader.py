#!/home/lyh/miniconda3/envs/myenv2/bin/python3.12
# cp2k/src/qs_tddfpt2_restart.F
from scipy.io import FortranFile # works for fortran generated binary file (UNFORMATTED)
from ctypes import c_int,c_double,c_float # 4Bytes, 8Bytes, 4Bytes
from dataclasses import dataclass
import sys

@dataclass
class TDWFN:
    '''use it as struct, no methods inside'''
    nstates:int  = None
    nspin  :int  = None
    nao    :int  = None
    nmo_occ:list = None
    eigval :list = None
    coef   :dict = None # dict of dict

def tdwfnReader(filename:str)->TDWFN:
    if not filename.endswith('.tdwfn'):
        raise ValueError("Invalid cp2k .tdwfn file")
    f = FortranFile(filename,'r',)
    nstates, nspin, nao = f.read_ints(dtype=c_int) # or dtype='i4'
    nmo_occ = []
    nmo_occ.extend(f.read_ints(dtype=c_int))
    eigval = f.read_reals(c_double)

    coef = {i:{} for i in range(nspin)} # eigenvectors
    for key in coef.keys():
        coef[key] = {j:[] for j in range(nstates)}
    for i in range(nspin):
        for j in range(nstates):
            for k in range(nmo_occ[j]):
                coef[i][j].append(f.read_reals(c_double))

    tdwfn = TDWFN(nstates=nstates,nspin=nspin,nao=nao,nmo_occ=nmo_occ,eigval=eigval,coef=coef)
    return tdwfn


if __name__ == "__main__":
    filename = sys.argv[1]
    tdwfn = tdwfnReader(filename)
    print(f"num of excited states: {tdwfn.nstates}\n \
          num of spin: {tdwfn.nspin}\n \
        num of ao: {tdwfn.nao}\n \
        num of occupied mo: {tdwfn.nmo_occ}\n \
        eigenvalues: {tdwfn.eigval}\n")
    for ispin in range(tdwfn.nspin):
        for istate in range(tdwfn.nstates):
            print(f"spin{ispin},state{istate}\n{tdwfn.coef[ispin][istate]}")
        print()