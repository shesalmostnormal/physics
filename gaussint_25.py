from __future__ import print_function, division
from scipy.special import gamma
from math import sqrt, exp, pi
from numpy import linspace, empty,zeros,linalg as LA, arange
from numpy.polynomial.laguerre import laggauss as lg
from matplotlib.pyplot import plot, show
import pylab as py
from sympy import *
from sympy.solvers import solve

"""
ugh

"""
d = int(raw_input('Basis (d x d): '))

row = range(0,d)

state = 1 # Principal quantum number
l = 0 # Angular momentum quantum number
m = 0 # Magnetic quantum number

### CONSTANTS (Natural Units) ###

# Rest energy calculation ~ 0.5109989461 MeV
M = 510998.9 # [Units: eV / c^2]
c = 1.0
const1 = M * c**2 

# Fine structure constant
const2 = 1973.3 # [Units: eV (length)] represents h_bar*c
const3 = 1/137. # [Units: dimensionless] represents e^2/const2

# Bohr radius
bohr = (const2/const1)*(1/const3) # [Units: Ang] ~ 0.55

# Rydberg, natural unit of energy, for measuring energy
	# levels of hydrogen:
Ry = (const1/2)*(const3)**2 # [Units: eV] ~ 13.3
#E_n = - Ry / n**2 # [Units:eV] ~ -13.6/n^2

#################
# NATURAL UNITS #
#################
#e = 8.5424546*10**(-2.)
e = 4*pi*const3
m_e = 0.51099906 #[eV]
a_0 = 2.6817268*10**(-4.) #[1/eV]
E = - e**2 / (2*a_0) #[eV]

##########################
# Variational Parameters #
##########################
Z = 1
alpha_0 = 1.
alpha = alpha_0 # Variational parameter

# Recurrence relaiton for Laguerre polynomials

def LaguerreL(a, b, x):
	if a==0:
		return 1
	if a==1:
		return 1 + b - x
	else:
		return ((2*(a-1) + 1 + b - x) * LaguerreL(a-1, b, x) -((a-1) + b) * LaguerreL(a-2, b, x)) / a

# Use Gamma functions for non-integer factorials

def norm(n, l):
	return sqrt(2*gamma(n + 1) / gamma(n + l + 3./2.))


"""
POTENTIAL ENERGY 

"""
#with substitution
def r_integrand(n_p, l, n, var, r):
	return norm(n,l) * norm(n_p,l) * var*0.5 * LaguerreL(n,l+0.5,r) \
	* LaguerreL(n_p,l+0.5,r) * r**(l) 

#without substitution
def r_integrand2(n_p, l, n, var, r):
	return norm(n,l) * norm(n_p,l) * var**(2.*l+3.) * exp(-var**2 * r**2) *\
	LaguerreL(n,l+0.5,var**2 * r**2) * LaguerreL(n_p,l+0.5,var**2 * r**2) * r**(2.*l+1.)


"""
KINETIC ENERGY

"""

#with substitution
def p_integrand(n_p, l, n, var, p):
	return norm(n,l) * norm(n_p,l) * (-1)**n * (-1)**n_p * 0.5 * var**2 *\
	 LaguerreL(n,l+0.5,p) * LaguerreL(n_p,l+0.5,p) #* p**(l + 1.5)


#without substitution
def p_integrand2(n_p, l, n, var, p):
	return norm(n,l) * norm(n_p,l) * (-1)**n * (-1)**n_p * var**(-2.*l-3.) *\
	exp(-p**2/var**2)*LaguerreL(n,l+0.5,p**2/var**2) \
	* LaguerreL(n_p,l+0.5,p**2/var**2) * p**(2.*l+4.)


"""
Gaussian-Laguerre Quadrature and Hamiltonian Matrix
"""
N = 20 ######### DO NOT CHANGE UNLESS POINTS COMPUTED AGAIN

# Gauss-Laguerre quadrature: on the interval [0,inf] and the weight function f(x) = exp(-x)
	# x = sample points, w = weights
x, w = lg(N)

# Points for non-integer exponent of x
	# NSolve[LaguerreL[20,1.5,x],x]

point = [0.237787, 0.704138, 1.40675, 2.3497, 3.53851, 4.98045, 6.68474, 8.66296, 10.9295, 13.5024, 16.4041, 19.6633, 23.3165, 27.4122, 32.0161, 37.2213, 43.1679, 50.0868, 58.4192, 69.2956]

#point = [0.0985091, 0.291263, 0.580559, 0.966699, 1.45006, 2.03109, 2.71037, 3.48856, 4.36641, 5.34479, 6.42445, 7.60508, 8.91048, 10.2568, 11.7942, 12.487, 14.0496, 15.1898, 15.9255, 17.1694, 19.5783, 19.9688, 23.4689, 23.9427, 26.1484, 29.1393, 31.8191, 34.1786, 38.3017, 42.3578, 44.4853, 48.7796, 50.8871, 53.1575, 59.7652, 64.2413, 69.1332, 78.0773, 84.8365, 95.1858, 100.9, 105.919, 117.283, 127.648, 139.615, 147.673, 160.86,  166.893, 177.136, 183.056]
weight = []

for s in range(0,N):
	weight.extend([gamma(N+1.5+1)*point[s]/(gamma(N+1.)*(N+1.)**2*\
		(LaguerreL(N+1, 1.5, point[s]))**2)])

varlist = arange(1.,10.,0.1)

matrix = empty((len(row),len(row)))
r_sum = empty((len(row),len(row)))
p_sum = empty((len(row),len(row)))

mass = 9.1*10**(-31)
q = 2.4*10**(-6)
planck = 1.054*10**(-27)

# Function for energy eigenvalues
def E(n,l,var,d):

	row = range(0, d) # Varying values of i and j, or size of matrix
	"""
	# Note: If printing J, the matrix should not be Hermitian: varying the alphas from the original
	#	will cause the vectors to no longer be the eigenstates of the Hamiltonian
	for i in row:
		for j in row:
			for k in range(0, len(w)):
				r_sum[i][j] += (r_integrand(j,l,i,alpha,x[k]) * w[k])
				p_sum[i][j] += (p_integrand(j,l,i,alpha,point[k]) * weight[k])
			matrix[i][j] = r_sum[i][j] + p_sum[i][j]

	print(p_sum[0][0],p_sum[1][1],p_sum[2][2])
	print(r_sum[0][0],r_sum[1][1],r_sum[2][2])
	"""

	matrix2 = empty((len(row),len(row)))
	r_sum2 = empty((len(row),len(row)))
	p_sum2 = empty((len(row),len(row)))

	x2 = 10.
	x1 = 0.01
	h = 0.001
	
	for i in row:
		for j in row:
			for k in range(0, 4000):
				p1 = p_integrand2(j,l,i,var,x1+k*h)
				p2 = p_integrand2(j,l,i,var,k*h+h+x1)
				r1 = r_integrand2(j,l,i,var,x1+k*h)
				r2 = r_integrand2(j,l,i,var,k*h+h+x1)
				#p_sum2[i][j] += ((p1+p2)*h/2.)
				#r_sum2[i][j] += ((r1+r2)*h/2.)
				p_sum2[i][j] += (0.5*(p1+p2)*h/2.)
				r_sum2[i][j] += (-(r1+r2)*h/2.)
			matrix2[i][j] = r_sum2[i][j] + p_sum2[i][j]

	print('\n',p_sum2[0][0],p_sum2[1][1],p_sum2[2][2])
	print(r_sum2[0][0],r_sum2[1][1],r_sum2[2][2])
	
	# Computes the eigenvalues and right eigenvectors of a square array.
	#print(r_sum)		
	#print(p_sum)
	u, v = LA.eigh(matrix2)

		# Note: Diagonalizing the matrix will reduce off-diagonal elements to zero, iff alpha =/= alpha_0

	# Reordering the eigenvalues: largest to smallest
	
	# Note: Smallest for ground state calculation

	idx = u.argsort()[::1]   
	u = u[idx]
	v = v[:,idx]

	# Diagonalizing the matrix will reduce off-diagonal 
		# elements to zero, iff alpha =/= alpha_0

	# v is represented by v[number of basis][state],
		# where state would remain consistent throughout
	
	return (u[state])


#E(state,0,1.,d)


E_list = empty(len(varlist))

for t, var in enumerate(varlist):
	E_list[t] += E(state,0,var,d)

print('\nState:',state,'\nEigenvalue:',mass*q**4/(2*planck**2),\
	'\nMin. Eigenvalue:',min(E_list),\
	'\nAlpha:',min(zip(E_list, varlist))[1])

y1, = py.plot(varlist,E_list)

py.xlabel('alpha')
py.ylabel('energy')
#py.xlim(0,2)
#py.ylim(-30,10)
py.show()


"""
cristina@cristina-X555LAB:~/Desktop/physics/coulomb$ python gaussint_21.py 
Basis (d x d): 3
1.500028006 3.50008751925 5.50009978091
1.1283791671 0.94031597258 0.836881215596

 1.49999171592 3.49863736298 5.44817214924
1.12826602961 0.940127408037 0.835957808782
cristina@cristina-X555LAB:~/Desktop/physics/coulomb$ python gaussint_21.py 
Basis (d x d): 7
1.500028006 3.50008751925 5.50009978091
1.1283791671 0.94031597258 0.836881215596

 1.49999171592 3.49863736298 5.44817214924
1.12826602961 0.940127408037 0.835957808782

"""