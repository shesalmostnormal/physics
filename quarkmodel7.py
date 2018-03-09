from __future__ import print_function, division
from scipy.special import gamma, roots_laguerre as lag, roots_genlaguerre as genlag, eval_genlaguerre as LaguerreL
from math import sqrt, exp, pi
from numpy import linspace, zeros,linalg as LA
from numpy.polynomial.laguerre import laggauss as lg
from matplotlib.pyplot import plot, show
import pylab as py

"""
Cern paper: http://cds.cern.ch/record/528258/files/0111325.pdf
Trying to match eigenvalues for bottom and top quarks.

"""

d = int(raw_input('Basis (d x d matrix): ')) #Basis
points = 20

state = int(raw_input('State: '))


row = range(0,d)

l = 0

#mass = 0.511 * 10**6.
#e = 1/137.

m_b =  4.5 #* 10**9. eV # Bottom quark
m_t = 170. #*10**9. eV # Top quark

mu = m_b/2. #* 10**6. # Reduced mass
coupling = 0.5 # Fine structure
sigma = 1. #* 10**9. # String tension

def norm(index):
	return sqrt(2*gamma(index + 1) / gamma(index + l + 3./2.))

"""
KINETIC ENERGY

"""

## Weights are w(x) = p^1.5 e^-p and f(x) = p^l
## So when l = 0
## Roots of LaguerreL[n,1.5,x]

#### but when l = 1
#### w(x) = p^2.5 e^-p and f(x) is still p^l
#### but now LaguerreL[n,2.5,x]


def p_integrand(n, m, var, p):
	return norm(n) * norm(m) * (-1)**n * (-1)**m * 0.5 * var**2. \
	* LaguerreL(n,l+0.5,p) * LaguerreL(m,l+0.5,p) * p**l

"""
POTENTIAL ENERGY 

"""

## Weights are w(x) = r^0 e^-x and f(x) = r^l
## So roots are of LaguerreL[n,0,x] (non-generalized)

#### At some point, l=1
#### Thus w(x) = r^1 e^-x and f(x) = r^0
#### So now find roots of LaguerreL[n,1.,x] (generalized)

def r_integrand(n, m, var, r):
	return norm(n) * norm(m) * var * 0.5 * LaguerreL(n,l+0.5,r) \
	* LaguerreL(m,l+0.5,r) * r**l

"""
LINEAR POTENTIAL ENERGY 

"""

## Weights are w(x) = r^1.0 e^-x and f(x) = r^l
## So roots are of LaguerreL[n,1.0,x]

def l_integrand(n, m, var, r):
	return norm(n) * norm(m) * var**(-1.) * 0.5 * LaguerreL(n,l+0.5,r) \
	* LaguerreL(m,l+0.5,r) * r**l



"""
Gaussian-Laguerre Quadrature
"""

# Gauss-Laguerre quadrature: on the interval [0,inf] and the weight function f(x) = exp(-x)
# x = sample points, w = weights
roots, weights = lag(points,mu=False)

# Points for non-generalized Laguerre l = 0 so superscript = 0.5
genroots, genweights = genlag(points,1.5,mu=False)

# Points for non-generalized Laguerre l = 0 so superscript = 0.5
linroots, linweights = genlag(points,1.0,mu=False)


def E(alpha):
	matrix = zeros((len(row),len(row)), dtype=float)
	p_sum = zeros((len(row),len(row)), dtype=float)
	r_sum = zeros((len(row),len(row)), dtype=float)
	l_sum = zeros((len(row),len(row)), dtype=float)


	# Note: If printing J, the matrix should not be Hermitian: varying the alphas from the original
	#	will cause the vectors to no longer be the eigenstates of the Hamiltonian

	for i in row:
		for j in row:
			for k in range(points):
				p_sum[i][j] += (p_integrand(i,j,alpha,genroots[k]) * genweights[k])
				r_sum[i][j] += (r_integrand(i,j,alpha,roots[k]) * weights[k])
				l_sum[i][j] += (l_integrand(i,j,alpha,linroots[k]) * linweights[k])
			matrix[i][j] =  (0.5*(1/mu) * p_sum[i][j] - coupling * r_sum[i][j] + sigma/5. * l_sum[i][j])


	# Diagonalizing the matrix will reduce off-diagonal 
	# elements to zero, iff alpha =/= alpha_0

	#Check (without parameter):
	#print(r_sum[0][0], p_sum[0][0]) Should equal
	#print(p_sum[0][1],r_sum[0][1]) 
	#print(norm(0,0)**2 * 3*sqrt(pi)/8) 

	u, v = LA.eig(matrix)

	# Sorting

	idx = u.argsort()[::1]   
	u = u[idx]
	v = v[:,idx]

	# v is represented by v[number of basis][state],
	# where state would remain consistent throughout 

	return(u[state])


"""
Golden Ratio Search
"""

accuracy = 1*10**(-6)	# Required accuracy
z = (1+sqrt(5))/2		# Golden ratio

# Initial positions of the four points
x1 = 1. * 10**-2.
x4 = 1. * 10**2.
x2 = x4 - (x4-x1)/z
x3 = x1 + (x4-x1)/z

f1 = E(x1)
f2 = E(x2)
f3 = E(x3)
f4 = E(x4)

# Main loop of search process
while x4-x1>accuracy:
	if f2<f3:
		x4,f4 = x3,f3
		x3,f3 = x2,f2
		x2 = x4 - (x4-x1)/z
		f2 = E(x2)
	else:
		x1,f1 = x2,f2
		x2,f2 = x3,f3
		x3 = x1 + (x4-x1)/z
		f3 = E(x3)


alphamin = 0.5*(x1+x4)

lam = 4*sigma/(m_b**2 * coupling**3) * (1/5.)
ep = 0.5*m_b*coupling**2*(-0.5+1.5*lam-1.5*lam**2+(27/4)*lam**3)

print('The minimum alpha falls at', alphamin)
print('The minimum eigenvalue is', E(alphamin))
print('Lambda=',lam)
print('Ground st. corr=',ep)
print('Ground st.=', m_b*2 + ep)