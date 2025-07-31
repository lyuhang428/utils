## Home-made minimalistic restricted Hartree-Fock calculation

\- utils.py

\- molecule.py

\- mmd2.py

\- hf.py

\- core.f90

are working code. 

This is home-made minimalistic Hartree-Fock demo, currently only works for restricted closed-shell molecules. The molecular integral (overlap, kinetic, nuclear-electron attraction, and electron repulsion integral) are implemented under `McMurchie-Davidson` scheme (in principle any arbitrary angular quantum number) and `Rys` quadrature (up to $d$-orbital). Fortran code does the computational heavy part, and interfaced to Python via `f2py` module. 

The code here is much slower than `Psi4` (one order of magnitude slower), but readability wins, and is very easy to be modified and implemented with new methods (i.e. wavefunction extraction). 

To use the code, Gnu scientific library `gsl` should be installed, and linked via `export LD_LIBRARY_PATH=/path/to/gsl_lib/`. 

`multiprocessing` module is used to do parallel computing, pretty straightford and handy. 


Author: lyh
