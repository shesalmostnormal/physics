from __future__ import print_function, division
from scipy.special import gamma
from numpy import pi, zeros, linalg as LA, exp, linspace, absolute as ab, sqrt, arange
from math import factorial
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend
import pylab as py

""" 
FINAL RESULT. WORKING FOR 10x10 and 30x30.

Graphing energy eigenvalues.
"""

#d = int(raw_input('Basis (d x d): '))
#n = int(raw_input('State (n): '))
#l = int(raw_input('Ang momentum (l): '))

K = 0.002 # Spring constant
M = 0.02 # Mass of particle
h_bar = 1

b0 = (M * K / h_bar) ** (1/4.)
b = 0.9*b0 # Vary by percent
omega = (K/M) ** (0.5)

print('True alpha =',b)

def LHS(i, j, b0, b, l): # LHS of Schrodinger equation

# Note: This function returns the matrix elements <j|H|i> 
	
	#z = complex(0,1)

	def norm(index, l):
		return sqrt(gamma(index+1) / gamma(index + l + 3./2.))

	def kron(i, j):
		if i == j:
			return 1.0  # If i = j
		return 0.0  # If i =/= j

    # Atomic units: hbar = 1

	normalization =  norm(i,l) * norm(j,l)

	coeff1 = 1/(2.* M) * b**2.0
	
	T = (((2*j) + l + (3./2.))*kron(i,j) - (j + l + (1/2.))*kron(i,j-1) - (j+1)*kron(i,j+1)) * (-1)**(i+j) #* z**(2.*l)
	
	coeff2 = b0**4.0 * 1/(2.*M) * b**(-2.)

	U = ( (2*j) + l + (3./2.))*kron(i,j) - (j + l + (1/2.))*kron(i,j-1) - (j+1)*kron(i,j+1)

	return  ((coeff1 * T) + (coeff2 * U)) #*normalization

# R matrix and P matrix diagonal elements should be the same test them individually.

# Function for energy eigenvalues
def E(n,l,b,d):

	row = range(0, d) # Varying values of i and j, or size of matrix
	Hmatrix = zeros((len(row), len(row)), dtype=float) # Initialization of matrix

	# Note: If printing J, the matrix should not be Hermitian: varying the alphas from the original
	#	will cause the vectors to no longer be the eigenstates of the Hamiltonian

	for i in row:
		for j in row:
       			Hmatrix[i][j] = LHS(i,j, b0, b, l)

	# Computes the eigenvalues and right eigenvectors of a square array.
	u, v = LA.eigh(Hmatrix)

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

	return (u[n])

d = int(raw_input('Basis (d x d): '))

# Eigenval = 2n + l + 3/2

statelist = [0,1,10]
varlist = arange(0.001,0.16,0.001)

E_list = zeros((len(statelist),len(varlist)),dtype=float)

for i, state in enumerate(statelist):
	for j, var in enumerate(varlist):
		E_list[i][j] += E(state,0,var,d)
	print('\nEnergy level:',state, '\nMin. Eigenvalue:',min(E_list[i]),'\nExact Eigenvalue:',omega*(2*state + 1.5))

y1, = py.plot(varlist,E_list[0])
y2, = py.plot(varlist,E_list[1])
y3, = py.plot(varlist,E_list[2])

py.legend([y1,y2,y3], ['E_0','E_1','E_10'])
py.xlabel('alpha')
py.ylabel('energy')
py.xlim(0,0.16)
py.ylim(0,10)
py.show()
