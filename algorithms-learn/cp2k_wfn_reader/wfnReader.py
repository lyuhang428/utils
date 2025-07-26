#!/home/lyh/miniconda3/envs/myenv2/bin/python3.12

# cp2k wfn, integer(4)/c_int, real(8)/c_double
# scipy.io.FortranFile is vert handy, if data in wfn is an array it automatically read it
# struct.unpack also works, but needs more efforts. In fortran unformatted file(binary), each record has leading and tailing record length
# i.e., the data in wfn read by scipy is [1,2,3,4,5], but struct will read [20,1,2,3,4,5,20], where 20=5 x sizeof(c_int)
# struct of cp2k wfn file:
'''
Number of atoms, num of spin, num of ao, max num of set, max num of shell: 3 1 24 3 3
Number of sets for each atom: 3 2 2
Number of shells for each [atom] and set: [ 3 2 1 ] [ 2 1 0 ] [ 2 1 0 ] 
Number of atomic orbitals for each atom, [set] and shell:
atom 0:	[ 1 1 1 ] [ 3 3 0 ] [ 5 0 0 ] 
atom 1:	[ 1 1 0 ] [ 3 0 0 ] [ 0 0 0 ] 
atom 2:	[ 1 1 0 ] [ 3 0 0 ] [ 0 0 0 ] 

# Information for spin 0: # this is not part of wfn
Number of molecular orbitals (MOs), num of homo, num of lumo, num of electron: 5 5 6 10
Eigenvalue for each MO: [ 0.000 0.000 0.000 0.000 0.000 ]
Electron occupancy for each MO: [ 2.000 2.000 2.000 2.000 2.000 ]
Coefficients for the MO/AO matrix:
MO 0:	[ -0.998 +0.000 -0.010 -0.068 +0.008 +0.006 +0.002 +0.001 -0.001 -0.000 -0.000 -0.001 -0.000 -0.001 +0.039 -0.017 +0.000 +0.000 +0.005 +0.048 -0.020 +0.001 +0.000 -0.005 ]
MO 1:	[ +0.028 +0.887 -0.106 -0.056 -0.010 +0.123 +0.075 -0.001 -0.022 -0.004 +0.000 -0.004 -0.000 +0.002 +0.193 -0.108 +0.016 -0.000 +0.030 +0.383 -0.168 +0.026 -0.000 -0.035 ]
MO 2:	[ -0.006 +0.156 -0.015 +0.004 -0.006 -0.718 +0.012 -0.000 +0.127 +0.022 +0.000 -0.000 -0.000 +0.001 +0.594 -0.195 +0.033 -0.000 +0.022 -0.509 +0.153 -0.025 -0.000 +0.012 ]
MO 3:	[ +0.001 +0.029 +0.016 +0.068 +0.919 -0.002 -0.001 +0.067 +0.000 +0.000 -0.019 +0.001 -0.000 +0.001 -0.028 +0.010 +0.001 +0.033 -0.002 -0.031 +0.011 +0.001 +0.033 +0.002 ]
MO 4:	[ +0.091 -0.208 -0.203 -0.800 +0.079 -0.013 +0.021 +0.006 +0.002 +0.000 -0.002 -0.011 +0.000 -0.015 +0.394 -0.146 -0.011 +0.003 +0.031 +0.374 -0.140 -0.012 +0.003 -0.031 ]
cp2k/src/qs_mo_io.F,  subroutine write_mo_set_low
'''
# incomplete sample usage of struct
'''
# with open(filename,'rb') as f:
#     content = f.read()
# print(f"size of file: {len(content)}")
# for i in range(0,28,4):
#     print(struct.unpack("i",content[i:i+4])[0],end=' ')
# print()
# for i in range(28,48,4):
#     print(struct.unpack("i",content[i:i+4])[0],end=' ')
# print()
# print(struct.unpack("i",content[48:52])[0])
# print
# for i in range(0,len(content),4):
#     print(struct.unpack("i",content[i:i+4])[0])
'''


#import struct # struct is a way in python to read binary, a general way
#import numpy as np
from scipy.io import FortranFile # if binary file is gnerated other than fortran, this way does not work
from ctypes import c_int,c_double
from dataclasses import dataclass
import sys

@dataclass
class WFN:
    '''use it as struct, no methods inside'''
    natom    :int  = None
    nspin    :int  = None
    nao      :int  = None
    nmo      :list = None
    nhomo    :list = None
    nlumo    :list = None
    ne       :list = None
    occupancy:dict = None
    coef     :dict = None

def wfnReader(filename:str)->WFN:
    if not filename.endswith('wfn'):
        raise ValueError("Invalid cp2k .wfn file")
    f = FortranFile(filename,'r')
    nmo = []
    nhomo = []
    nlumo = []
    ne = []
    occupancy = {}
    coef = {}
    natom, nspin, nao, _, __ = f.read_ints(dtype=c_int)
    _ = f.read_ints(c_int)
    _ = f.read_ints(c_int)
    _ = f.read_ints(c_int)
    
    if nspin > 1:
        for ispin in range(nspin):
            occupancy[ispin] = []
            coef[ispin] = []
            _nmo, _nhomo, _nlumo, _ne = f.read_ints(c_int)
            nmo.append(_nmo)
            nhomo.append(_nhomo)
            nlumo.append(_nlumo)
            ne.append(_ne)
            tmp = f.read_ints(c_double)
            occupancy[ispin] = tmp[int(len(tmp)/2):]
            for _ in range(nmo[ispin]):
                coef[ispin].append(f.read_reals(c_double))
    elif nspin == 1:
            occupancy[0] = []
            coef[0] = []
            _nmo, _nhomo, _nlumo, _ne = f.read_ints(c_int)
            nmo.append(_nmo)
            nhomo.append(_nhomo)
            nlumo.append(_nlumo)
            ne.append(_ne)
            tmp = f.read_ints(c_double)
            occupancy[0] = tmp[int(len(tmp)/2):]
            for _ in range(nmo[0]):
                coef[0].append(f.read_reals(c_double))
    wfn = WFN(natom=natom,nspin=nspin,nao=nao,nmo=nmo,nhomo=nhomo,nlumo=nlumo,ne=ne,occupancy=occupancy,coef=coef)
    return wfn



if __name__ == "__main__":
    #filename = './rks.wfn'
    filename = sys.argv[1]
    wfn = wfnReader(filename=filename)
    print(f"natom: {wfn.natom}\n \
            nspin: {wfn.nspin}\n \
            nao: {wfn.nao}\n \
            nmo: {wfn.nmo}\n \
            nhomo: {wfn.nhomo}\n \
            nlumo: {wfn.nlumo}\n \
            occupancy: {wfn.occupancy}")
    for key in wfn.coef.keys():
        print(key)
        for i in wfn.coef[key]:
            print(i)
