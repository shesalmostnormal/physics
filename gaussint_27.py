from __future__ import print_function, division
from scipy.special import gamma
from math import sqrt, exp, pi
from numpy import linspace, empty,zeros,linalg as LA, arange
from numpy.polynomial.laguerre import laggauss as lg
#from matplotlib.pyplot import plot, show
#import pylab as py
from sympy import *
from sympy.solvers import solve

import gaulag

"""
F2PY

"""
dimen = int(raw_input('Basis: '))
N = int(raw_input('Points: '))


row = range(0,dimen)

state = 1 # Principal quantum number
l = 0 # Angular momentum quantum number
m = 0 # Magnetic quantum number

# Recurrence relation for Laguerre polynomials

def LaguerreL(a, b, x):
	if a==0:
		return 1
	if a==1:
		return 1 + b - x
	else:
		return ((2*(a-1) + 1 + b - x) * LaguerreL(a-1, b, x) -((a-1) + b) * LaguerreL(a-2, b, x)) / a

# Use Gamma functions for non-integer factorials

def norm(index, l):
	return sqrt(2*gamma(index + 1.) / gamma(index + l + 3./2.))


"""
POTENTIAL ENERGY 

"""
#with substitution
def r_integrand(i, l, j, var, r):
	return norm(i,l) * norm(j,l) * var*0.5 * LaguerreL(i,l+0.5,r) \
	* LaguerreL(j,l+0.5,r) * r**(l) 

"""
KINETIC ENERGY

"""

#with substitution
def p_integrand(i, l, j, var, p):
	return norm(i,l) * norm(j,l) * (-1)**(i+j) * 0.5 * var**2 *\
	 LaguerreL(i,l+0.5,p) * LaguerreL(j,l+0.5,p) #* p**(l + 1.5)


"""
Gaussian-Laguerre Quadrature and Hamiltonian Matrix
"""

# Gauss-Laguerre quadrature: on the interval [0,inf] and the weight function f(x) = exp(-x)
	# x = sample points, w = weights
x, w = lg(N)

# Points for non-integer exponent of x

point = []

def point_list(N):
	if N == 2: # NSolve[LaguerreL[2,1.5,x],x]
		sample2 = [1.62917, 5.37083]
		return sample2
	if N == 3: # NSolve[LaguerreL[3, 1.5, x], x]
		sample3 = [1.2204, 3.80888, 8.47072]
		return sample3
	if N == 4: # NSolve[LaguerreL[4,1.5,x],x]
		sample4 = [0.978507, 2.99038, 6.3193, 11.7118]
		return sample4
	if N == 5: # NSolve[LaguerreL[5,1.5,x],x]
		sample5 = [0.817632, 2.47233, 5.11601, 9.04415, 15.0499]
		return sample5
	if N == 6: # NSolve[LaguerreL[6, 1.5, x], x]
		sample6 = [0.702602, 2.11141, 4.32071, 7.48505, 11.921, 18.4593]
		return sample6
	if N == 7: # NSolve[LaguerreL[7, 1.5, x], x]
		sample7 = [0.616146, 1.84433, 3.74877, 6.42247, 10.0327, 14.912, 21.9236]
		return sample7
	if N == 8: # NSolve[LaguerreL[8, 1.5, x], x]
		sample8 = [0.548742, 1.63818, 3.31504, 5.64031, 8.71566, 12.7177, 17.9924, 25.4319]
		return sample8
	if N == 9: # NSolve[LaguerreL[9, 1.5, x], x]
		sample9 = [0.494692, 1.47402, 2.97366, 5.03608, 7.72862, 11.1593, 15.5124, 21.1451, 28.9762]
		return sample9
	if N == 10: # NSolve[LaguerreL[10, 1.5, x], x]
		sample10 = [0.450371, 1.34008, 2.69741, 4.55325, 6.95492, 9.97444, 13.7248, 18.3967, 24.3575, 32.5505]
		return sample10
	if N == 20: # NSolve[LaguerreL[20, 1.5, x], x]
		sample20 = [0.237787, 0.704138, 1.40675, 2.3497, 3.53851, 4.98045, 6.68474, 8.66296, 10.9295, 13.5024, 16.4041, 19.6633, 23.3165, 27.4122, 32.0161, 37.2213, 43.1679, 50.0868, 58.4192, 69.2956]
		return sample20

#f2py -m gaulag -c gaulag.f
#print(gaulag(10,10,2,1.5))

weight = []
point = point_list(N)
s = 0
for s in range(0,N):
	weight.extend([gamma(N+1.5+1)*point[s]/(gamma(N+1.)*(N+1.)**2*\
		(LaguerreL(N+1, 1.5, point[s]))**2)])

print(weight)
# Function for energy eigenvalues

matrix = empty((len(row),len(row)))
r_sum = empty((len(row),len(row)))
p_sum = empty((len(row),len(row)))

e = 1.
h = 1.
m_e = 1.
"""
def E(n,l,var,d):

	row = range(0, dimen) # Varying values of i and j, or size of matrix


	# Note: If printing J, the matrix should not be Hermitian: varying the alphas from the original
	#	will cause the vectors to no longer be the eigenstates of the Hamiltonian
	for i in row:
		for j in row:
			for k in range(0, len(w)):
				
				p_sum[i][j] += (p_integrand(i,l,j,var,point[k]) * weight[k])
				r_sum[i][j] += (r_integrand(i,l,j,var,x[k]) * w[k])			
				
			matrix[i][j] = -e*r_sum[i][j] + 0.5*(h/m_e)*p_sum[i][j]

#	print(p_sum[0][0],p_sum[1][1],p_sum[2][2])
	
#	print(r_sum[0][0],r_sum[1][1],r_sum[2][2])

#	print(matrix[0][0],matrix[1][1],matrix[2][2])
	
	#print(matrix)

	# Computes the eigenvalues and right eigenvectors of a square array.
	#print(r_sum)		
	#print(p_sum)
	u, v = LA.eigh(matrix)

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

#	print(u[state])
	return (u[state])
"""

varlist = arange(0.1,0.3,0.01)

#E(state,0,0.2685,dimen)

"""
E_list = empty(len(varlist))

for t, var in enumerate(varlist):
	E_list[t] += E(state,0,var,dimen)

print('\nState:',state,'\nEigenvalue:',1.,\
	'\nMin. Eigenvalue:',min(E_list),\
	'\nAlpha:',min(zip(E_list, varlist))[1])

#y1, = py.plot(varlist,E_list)

py.xlabel('alpha')
py.ylabel('energy')
#py.xlim(0,2)
#py.ylim(-30,10)
#py.show()
"""