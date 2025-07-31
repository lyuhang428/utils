# Author: lyh

import numpy as np
import numpy.polynomial.legendre as L
import autograd.numpy as agnp
from typing import override


class Grid:
    '''
    Gauss-Legendre grid and quadrature
    '''
    def __init__(self, norder:int=5):
        self.norder = norder
        self.roots, self.weights = L.leggauss(self.norder)
        self.basis = L.legvander(self.roots, self.norder-1)    # Vandermonde matrix w.r.t. Legendre polynomial
        self._u, self._s, self._vt = np.linalg.svd(self.basis) # Penrose pseudo inverse, spectral method (coefficient space)
        self.roots_transformed = Grid.transform(self.roots)    # map self.roots ∈ [-1,1] to self.roots_transformed ∈ [0.big_number)
        self.jacobian = Grid.transform_der(self.roots)         # map x -> T(x), T'(x) = dT(x)/dx, this is Jacobian det under change of variable

    @staticmethod
    def transform(x: np.ndarray) -> np.ndarray:
        '''map x ∈ [-1, 1] to T(x) ∈ [0, Inf]'''
        left = 1e-3
        right = 1e4
        alpha = np.log(right / left)
        return left * (np.exp(alpha * (1 + x) / 2) - 1)

    @staticmethod
    def transform_der(x:np.ndarray) -> np.ndarray:
        '''derivative of `transform` function w.r.t. x'''
        left = 1e-3
        right = 1e4
        alpha = np.log(right / left)
        return left * alpha / 2. * np.exp(alpha * (1 + x) / 2.)

    def legquad(self, fnvals:np.ndarray) -> np.ndarray:
        r'''Gauss-Legendre quadrature
        
        Integrand functions are discretized on [0, inf). 
        
        A single row of `fnvals` is one integrand function at all points; a single col of `fnvals` is all integrand functions at a given point

        Note this method only works for definite integral in [-1,1]. For [0,inf), `fnvals` should be multiplied with jacobian det, as 

        \int f(y) dy = \int f[T(x)] d T(x) = \int f[T(x)] T^\prime (x) dx
        
        Parameters
        ----------
        fnvals
            f0   |f0(x0)   f0(x1)   f0(x2) ...     f0(xN-1)|

            f1   |f1(x1)   f1(x2)   f1(x3) ...     f1(xN-1)|

                 |...            ...                    ...|

            fm-1 |fm-1(x1) fm-1(x2) fm-1(x3) ... fm-1(xN-1)|

        Return
        ------
            If there's only one integrand function, `fnvals` has shape (N,), returns [val]. If there're many integrand functions, returns a ndarray with shape (M,), 
            elements correspond to the integrals of each integrand functions
        '''
        dims = fnvals.shape
        if len(dims) == 1:
            if dims[0] != len(self.weights):
                raise ValueError('Dimensions do not match')
            quad = np.sum(fnvals * self.weights)
        elif len(dims) == 2:
            if dims[1] != len(self.weights):
                raise ValueError('Dimensions do not match')
            quad = np.sum(fnvals * self.weights[None,:], axis=1)
        else:
            raise NotImplementedError('3D+ tensors not supported')
        
        return quad

    def tocoefs(self, fnvals:np.ndarray) -> np.ndarray:
        r'''
        Vander Coef = Fn <=> VC = F, C = V^-1 F
        
        `Vander` may have no inverse, using svd to do Penrose pseudo-inverse, numerically more stable than direct inverse

        f(x) = \sum_{n=0} c_n P_n(x). Compute coefficient {c_n}, turn integral-differential operation into arithmetric operation
        
        Working in the coefficient space is called `spectral method`
        '''
        coefs = (self._vt.T @ np.diag(1 / self._s) @ self._u.T) @ fnvals # 彭罗斯伪逆
        return coefs
    
    def tofnvals(self, coefs:np.ndarray) -> np.ndarray:
        '''
        Given coefficients of Legendre polynomial, return fnvals at Legendre roots. P C = Fn
        '''
        if len(coefs) > self.norder:
            import warnings
            warnings.warn("not enough basis to convert fnvals to coefs", category=Warning)
            _coefs = coefs[:self.norder]
            fnvals = self.basis @ _coefs
        elif len(coefs) < self.norder:
            width = self.norder - len(coefs)
            _coefs = np.pad(coefs, pad_width=(0, width), constant_values=(0., 0.), mode="constant")
            fnvals = self.basis @ _coefs
        else:
            fnvals = self.basis @ coefs
        
        return fnvals

    def derivative(self, fnvals:np.ndarray) -> np.ndarray:
        '''
        Project f onto Legendre polynomial basis functions, getting coefficients, then make use of the 
        three-term recursion relation of orthogonal polynomial to turn differential operation to arithmetric 
        operation (spectral method). Once new coefficients are known, f' values at Legendre roots will be known
        '''
        if self.norder != len(fnvals):
            raise ValueError("Lengths do not match")
        coefs = self.tocoefs(fnvals)
        coefs_der = np.pad(L.legder(coefs), pad_width=(0, 1), constant_values=(0, 0))
        fnvals_der = self.basis @ coefs_der
        
        return fnvals_der
    
    def antiderivative(self, fnvals:np.ndarray) -> np.ndarray:
        '''
        Spectral method to compute the indefinite integral of f and its values at Legendre roots. 
        Note this is indefinite integral, boundary conditions needed to uniquely determine the integral
        '''
        if self.norder != len(fnvals):
            raise ValueError('Lengths do not match')
        coefs = self.tocoefs(fnvals)
        coefs_int = L.legint(coefs)[:-1]
        fnvals_int = self.basis @ coefs_int

        return fnvals_int



class Basis:
    r'''e^{-\alpha r^2} is the true basis function. Multiplying r gives r e^{-\alpha r^2}. This change make evaluation 
    of overlap, kinetic, external matrices easier, as the r^2 term shown in integral under spherical coordinate is absorbed 
    into the modified basis function already. This is just the radial part of basis function'''

    def __init__(self, grid:Grid, alphamin:float=1e-6, alphamax:float=1e7, nbasis:int=80):
        r'''
        Different `alpha` gives different basis functions. Basis functions are discretized on transformed 
        Legendre roots (mapped from [-1,1] to [0,inf)
        '''
        self.nbasis = nbasis
        self.grid = grid
        self.alphas = np.power(10, np.linspace(np.log10(alphamin), np.log10(alphamax), nbasis)) # logarithmically even

        # x exp(-a x^2). But true basis function is exp(-a x^2), `x` here is for spherical coordinates integration
        self.fnvals = grid.roots_transformed * np.exp(-np.outer(self.alphas, grid.roots_transformed**2))
        self.fnvals *= (np.sqrt(np.sqrt(self.alphas))**3 * np.sqrt(np.sqrt(2 / np.pi) * 8))[:,None] # normalization
        if not np.allclose(grid.legquad(self.fnvals**2 * grid.jacobian), 1.):
            # test if basis functions are normalized
            raise ValueError("Basis functions are not normalized")
        
        alpha_sums  = np.add.outer(self.alphas, self.alphas) # matrix outer addition
        alpha_prods =     np.outer(self.alphas, self.alphas) # matrix outer product

        # overlap, external, kinetic operator are independent of electron density, they can be precomputed once basis set is given
        # Hartree and XC operator depends on e density, and will be updated in each SCF loop
        self.olp = (2 * np.sqrt(2)) * (alpha_prods) ** 0.75 / alpha_sums**1.5    # overlap
        self.kin_rad = np.sqrt(72) * alpha_prods**1.75 / alpha_sums**2.5         # kinetic radial part
        self.kin_ang = np.sqrt(32) * (alpha_prods) ** 0.75 / np.sqrt(alpha_sums) # kinetic angular part, l=1 only
        self.ext = -np.sqrt(32 / np.pi) / alpha_sums * alpha_prods**0.75         # external field (nuclear-electron attraction)

    @override
    def __len__(self) -> int:
        return self.nbasis
    

class LDA:
    '''LDA exchange functional. This class is not for instantiation

    TODO: include LDA correlation functional (VWN V for instance). GGA exchange-correlation functional is the next next step.
    '''
    def __init__(self, *args, **kwds):
        pass

    @staticmethod
    def exc(rho:np.ndarray) -> np.ndarray:
        r'''
        Exchange functional at LDA level: 

        Exc[n(r)] = -3/4 * (3/pi)^{1/3} \int n(r)^{4/3} r^2 \sin\theta dr d\theta d\phi. This is excahnge energy, scalar
        
        εxc(r) = -3/4 * (3/pi)^{1/3} n(r)^{4/3}. Exchange energy density on the grid, vector

        Returns integrand, quadrature is done outside here
        '''
        return -3/4 * (3/np.pi)**(1/3) * rho**(4/3)
    
    @staticmethod
    def vxc(rho:np.ndarray) -> np.ndarray:
        r'''
        Functional derivative gives exchange potential Vxc(r) = \delta Exc[n(r)] / \delta n(r). 
        But at LDA level, Vxc is simply common function derivative of εxc w.r.t. n(r)
        '''
        return -(3/np.pi)**(1/3) * rho**(1/3)


def build_rho(econfs:list[list[int]], e_orbitals:list[np.ndarray], grid:Grid, basis:Basis) -> np.ndarray:
    '''
    Arguments
    ---------
    econfs : list[list[int]]
        i.e. [[2,2], [2]], 1s2 2s2 2s2, stands for carbon atom
    e_orbitals : list[np.ndarray]
        truncated eigenvecs of scipy.linalg.eigh(a,b, subset_by_index). Only occupied orbitals are stored. 
        Each ndarray in `e_orbitals` correspond to each angular momentum quantum number
    '''
    rho = np.zeros((grid.norder, )) # electron density on the grid
    # Modified basis functions have prefactor r, now remove it and times angular part (only spherically symmetrical s-orbital spherical harmonics is considered)
    for econf, e_orbital in zip(econfs, e_orbitals):
        u = e_orbital.T @ basis.fnvals # use coefficients and basis functions to build orbitals
        r = u / grid.roots_transformed / np.sqrt(4 * np.pi) 
        rho += np.dot(econf, r**2)
    
    return rho


def poisson_solver(rho:np.ndarray, grid:Grid) -> np.ndarray:
    r'''
    Spherically symmetrical Poisson equation solver (radial part)

    Let u(r) = r V(r), Poisson equation is converted from \nabla^2 V(r) = -4 \pi \rho(r) to u''(r) = -4 \pi r rho(r)

    Only electron potential is considered, so at r=0 there's no singularity, electron density and potential at r=0 are smooth and differentiable
    
    If nuclear potential is also considered, things become complicated as it's ∞ at r=0
    '''

    ne = grid.legquad(4 * np.pi * grid.roots_transformed**2 * rho * grid.jacobian)  # number of electrons
    u_prime = grid.antiderivative(grid.roots_transformed * rho * grid.jacobian)     # u'(r) = r * rho(r)
    u_prime -= u_prime[0]                                                           # force u'(0) = 0, boundary condition
    u = grid.antiderivative(u_prime * grid.jacobian)                                # u(r) = \int_0^r u'(r') dr'
    u -= u[0]                                                                       # force u(0) = 0, boundary condition
    epot = -4 * np.pi * u
    epot += ((ne - epot[-1]) / grid.roots_transformed[-1]) * grid.roots_transformed # fulfill asymptotic behaviour at r->∞
    epot /= grid.roots_transformed                                                  # V(r) = u(r) / r

    return epot


class AtomDFT:
    def __init__(self, atom:str, grid:Grid, basis:Basis):
        from leconfig import SYMBOLS, E_CONFIG_ARRAY
        self.grid    : Grid  = grid
        self.basis   : Basis = basis
        self.atom    : str   = atom
        self.charge  : int   = SYMBOLS[atom]
        self.occups  : list  = E_CONFIG_ARRAY[atom]
        self.vol     : float = 4. * np.pi * self.grid.roots_transformed**2 # r^2 \sin\theta dr d\theta d\phi = 4\pi r^2
        self.maxqn   : int   = len(self.occups)-1                          # maximum angular momentum quantum number

        if SYMBOLS[atom] != np.concatenate(self.occups).sum():
            # test if it's neutral atom
            raise ValueError('number of electron and atomic charge do not match!')
        self.ne : int = SYMBOLS[atom]

    def scf(self, maxiter:int=50, mixing:float=0.5, eps_scf:float=1e-6):
        from scipy.linalg import eigh             # generalized hermitian eigendecomposition FC = (SC)ε
        rho      = np.zeros((self.grid.norder, )) # ectron density only on the grid, INITIAL GUESS
        vhartree = np.zeros((self.grid.norder, )) #   Coulomb potential on the grid, INITIAL GUESS
        vxc      = np.zeros((self.grid.norder, )) #        xc potential on the grid, INITIAL GUESS
        etot_old = float('inf')
        fock_old = []

        for iscf in range(maxiter):
            e_ext = 0.
            e_kin_rad = 0.
            e_kin_ang = 0.
            e_orbitals = [] # energy on each orbital (orbitals are linear combination of basis functions)
            jxc = (self.basis.fnvals * (vhartree+vxc) * self.grid.weights * self.grid.jacobian) @ self.basis.fnvals.T # Coulomb operator + xc operator, shape (nbasis,nbasis)

            for qn in range(self.maxqn+1):
                # kinetic angular part relies on angular quantum number, and so is fock
                fock = self.basis.kin_rad + self.basis.ext*self.charge + jxc # Fock operator, shape (nbasis,nbasis)
                fock += (qn * (qn + 1) / 2) * self.basis.kin_ang             # kinetic angular part, kin_ang is for l=1 only

                # mix in previous Fock matrix
                if iscf == 0:
                    fock_old.append(fock) # [l=0 old, l=1 old, l=2 old, ...]. It's meaningless to mix up Focks with different angular quantum number
                else:
                    fock = mixing * fock + (1 - mixing) * fock_old[qn]
                    fock_old[qn] = fock.copy()

                subset = (0, len(self.occups[qn]) - 1)                          # in principle nbasis orbitals can be formed, but only occupied orbitals are stored
                vals, vecs = eigh(fock, self.basis.olp, subset_by_index=subset) # FC = (SC)ε. F is a, S is b, C is vecs. Generalized eigendecomposition. vals sorted in descending order
                e_orbitals.append((vals, vecs))                                 # update for each angular quantum number, s=0, p=1, d=2, etc.
                e_kin_rad += (np.array(self.occups[qn]) * np.diag(vecs.T @ self.basis.kin_rad @ vecs)).sum()
                e_ext += self.charge * (np.array(self.occups[qn]) * np.diag(vecs.T @ self.basis.ext @ vecs)).sum()
                if qn > 0: # qn=0, s orbital has no contribution to angular part kinetic energy
                    e_kin_ang += (qn * (qn + 1) / 2) * (self.occups[qn] * np.diag(vecs.T @ self.basis.kin_ang @ vecs)).sum()

            # once orbitals are built, construct electron density on the grid
            rho = build_rho(self.occups, [eig[1] for eig in e_orbitals], self.grid, self.basis)
            if not np.isclose(self.grid.legquad(rho * self.vol * self.grid.jacobian), self.ne):
                raise ValueError('integrate over electron density on the grid does not give total number of electron!')

            # once rho is built, vhartree, vxc are updated and go to the next iter
            vhartree = poisson_solver(rho, self.grid)                                           # Coulomb potential on the grid, vector
            e_hartree = 0.5 * self.grid.legquad(vhartree * rho * self.vol * self.grid.jacobian) # Coulomb energy, scalar

            exc, vxc = LDA.exc(rho), LDA.vxc(rho)                         # exchange energy density and exchange potential on the grid, vector
            e_xc = self.grid.legquad(exc * self.vol * self.grid.jacobian) # exchange energy, scalar

            etot = e_kin_ang + e_kin_rad + e_ext + e_hartree + e_xc
            print(f"Iter: {iscf:>3} \t etot: {etot:>10.6f} \t kin_rad: {e_kin_rad:>10.6f} \t kin_ang: {e_kin_ang:>10.6f} \t e_hartree: {e_hartree:>10.6f} \t e_xc: {e_xc:>10.6f} \t e_ext: {e_ext:>10.6f}")

            if abs(etot - etot_old) <= eps_scf:
                print(f"SCF converged after {iscf:>3} iters")
                break
            etot_old = etot.copy()

            if iscf == maxiter-1:
                import warnings
                warnings.warn(f'SCF did not converge in {iscf:>3}  iters!', category=Warning)
            
        print(f"Total energy: {etot}")





def testjxc():
    '''how Coulomb-XC operator is constructed'''
    _pot = np.random.normal(0., 1., (256,))

    grid = Grid(norder=256)
    basis = Basis(grid)
    jxc = (basis.fnvals * _pot * grid.weights * grid.jacobian) @ basis.fnvals.T

    print(jxc.shape)
    

def main():
    grid = Grid(256)
    basis = Basis(grid, )
    for elem in  ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K"]:
        atom = AtomDFT(elem, grid, basis, )
        print(f"Atom: {elem}")
        atom.scf()
        print()




if __name__ == '__main__':
    main()


