from __future__ import print_function, division
from scipy.special import gamma
from numpy.polynomial.laguerre import laggauss as lg
from sympy import *
from sympy.solvers import solve
from numpy import *

"""
Copying Gauss-Laguerre formula for **generalized**

WORKS.

"""

num = 6
alpha = 1.

x = Symbol('x')


# Recurrence relaiton for Laguerre polynomials

def LaguerreL(N, alpha, x):
	if N==-1:
		return 0.
	if N==0:
		return 1.
	if N==1:
		return 1 + alpha - x
	else:
		return ((2*(N-1) + 1 + alpha - x) * LaguerreL(N-1, alpha, x) -((N-1) + alpha) * LaguerreL(N-2, alpha, x)) / N


varlist = arange(1.,2.,1)


def lg(alpha):

	point = []	
	weight = []

	for i in range(1,num):

		R = LaguerreL(i,alpha,x)
		Laguerreprime = R.diff(x)
		point+= real(solve(x * Laguerreprime + ((i+alpha) * LaguerreL(i-1,alpha,x)), x))

	k = 0
	index = [1.,2.,2.,3.,3.,3.]

	for j in index:
		#print(k)
		weight.extend([(point[k]*gamma(j+alpha)) / (j * gamma(j) * (j+alpha)*(LaguerreL(j-1,alpha,point[k]))**2)])
		#print(weight)
		k += 1

	print(point,weight)

	return point,weight

point,weight = lg(1.)