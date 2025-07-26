This is home-made Hartree-Fock quantum chemistry code. 

Currently it only does closed-shell calculation. 

The molecular integral engine is also home-made, third party libraries such as Psi4 and PySCF are not used

The `McMurchie-Davidson` scheme is used to evaluate the `overlap`, `kinetic`, `nuclear attraction`, and `electron repulsion integral`. In principle can calculate any angular quantum number, but only up to f-orbital is tested. 

For `eri` calculation, `Rys` quadrature is also implemented, but can only evaluate up to d-orbital (i.e. $[dd|dd]$)

For better performance, Fortran code is used to do the core part, and `f2py` is used to provide Python wrapper. 

The whole code is easy to use, I try to make everything **Header only**, so no annoying compilation is needed.

Python `NumPy` module provides many handy APIs for numerical operations, but such convenient library is lacking in Fortran and C/C++. This code uses `gsl` to do some numerical calculations. 

In the SCF part, DIIS algorithm is used to accelerate convergence. Only very naive `multiprocessing` level parallelization is used. Overall, this code just computes single point energy of a closed shell molecule. Though not implemented, wavefunction of given molecule under given basis set can be obtained very straightford. 

Author: lyh