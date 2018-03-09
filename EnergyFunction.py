from __future__ import print_function, division
import numpy as np
from numpy.linalg import eigvalsh as GetEigen
from scipy.special import gamma, roots_laguerre as GetRoots, roots_genlaguerre as GetGenRoots, eval_genlaguerre as LaguerreL


"""
Energy Eigenvalues.
Quark and meson masses from calculatemasses.py
"""

# Basis is dimensions of matrix.
#basis = int(raw_input('Basis: '))
basis = 30
size = range(0,basis)
dimensions = (len(size),len(size))

# Number of states; angular momentum.
state = 3
l = 0

# Number of sample points for Gauss-Laguerre integration.
points = 30

# Constants.
coupling = (4./3.)*0.2 # Fine structure
sigma = 1000. *197.#* 10**6. # String tension


"""
NORMALIZATION
"""

def norm(index):
	return np.sqrt(2*gamma(index + 1) / gamma(index + l + 3./2.))

"""
KINETIC ENERGY

"""

## Weights are w(x) = p^1.5 e^-p and f(x) = p^l
## So when l = 0
## Roots of LaguerreL[n,1.5,x]

#### but when l = 1
#### w(x) = p^2.5 e^-p and f(x) is still p^l
#### but now LaguerreL[n,2.5,x]


def p_int(n, m, var, p):
	return (-1)**n * (-1)**m * 0.5 * var**2. \
	* LaguerreL(n,l+0.5,p) * LaguerreL(m,l+0.5,p) * p**l

"""
POTENTIAL ENERGY 

"""

## Weights are w(x) = r^0 e^-x and f(x) = r^l
## So roots are of LaguerreL[n,0,x] (non-generalized)

#### At some point, l=1
#### Thus w(x) = r^1 e^-x and f(x) = r^0
#### So now find roots of LaguerreL[n,1.,x] (generalized)

def r_int(n, m, var, r):
	return var * 0.5 * LaguerreL(n,l+0.5,r) \
	* LaguerreL(m,l+0.5,r) * r**l

"""
LINEAR POTENTIAL ENERGY 

"""

## Weights are w(x) = r^1.0 e^-x and f(x) = r^l
## So roots are of LaguerreL[n,1.0,x]

def l_int(n, m, var, r):
	return var**(-1.) * 0.5 * LaguerreL(n,l+0.5,r) \
	* LaguerreL(m,l+0.5,r) * r**l



"""
Gaussian-Laguerre Quadrature
"""

# Gauss-Laguerre quadrature: on the interval [0,inf] and the weight function f(x) = exp(-x)
# x = sample points, w = weights
roots, weights = GetRoots(points,mu=False)

# Points for non-generalized Laguerre l = 0 so superscript = 0.5
genroots, genweights = GetGenRoots(points,1.5,mu=False)

# Points for non-generalized Laguerre l = 0 so superscript = 0.5
linroots, linweights = GetGenRoots(points,1.0,mu=False)

"""
CALCULATING ENERGY EIGENVALUES FUNCTION DECLARATION
"""

def Calculate_Energy(qm_1,qm_2,meson_m,alpha):

	Hmatrix = np.zeros(dimensions)
	p_sum = np.zeros(dimensions)
	r_sum = np.zeros(dimensions)
	l_sum = np.zeros(dimensions)

	for i in size:
		for j in size:
			for k in range(points):
			
				p_sum[i][j] += (p_int(i,j,alpha,genroots[k]) * genweights[k])
				r_sum[i][j] += (r_int(i,j,alpha,roots[k]) * weights[k])
				l_sum[i][j] += (l_int(i,j,alpha,linroots[k]) * linweights[k])
			
			Hmatrix[i][j] =  (norm(i)*norm(j))*(0.5*(1/meson_m) * p_sum[i][j] \
				- coupling * r_sum[i][j] + sigma * l_sum[i][j] + qm_1 + qm_2)


	EigenList = GetEigen(Hmatrix)
	
	# Sorting for eigenvalues
	#idx = u.argsort()[::1]   
	#u = u[idx]
	
	u = []
	for v in range(state):
		u.append(EigenList[v])
	
	return(u)