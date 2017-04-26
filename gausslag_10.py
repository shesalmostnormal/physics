from __future__ import print_function, division
from scipy.special import gamma
from numpy.polynomial.laguerre import laggauss as lg
from sympy import *
from sympy.solvers import solve
from numpy import *

"""
Copying Gauss-Laguerre formula for **generalized**

WORKS FOR ONE ALPHA AND VARIOUS INDEX.

"""

num = 3
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

f = open("glfilegen_10.txt","w")

point = []
weight = []

for i, alpha in enumerate(varlist):
	for j in range(1,num):

		R = LaguerreL(j,alpha,x)
		Laguerreprime = R.diff(x)
		point += solve(x * Laguerreprime + ((j+alpha) * LaguerreL(j-1,alpha,x)), x)

k = 0
index = [1.,2.,2.]

for alpha in varlist:
	for j in index:

		R = LaguerreL(j,alpha,x)
		Laguerreprime = R.diff(x)
		S = Laguerreprime.subs(x,point[k])
		weight.append(- gamma(j+alpha)/ (j * gamma(j) * LaguerreL(j-1,alpha,point[k]) * S))
		
		print(j,alpha,point[k],weight[k])
		
		k += 1

f.write(str(point) + '\n')
f.write(str(weight) + '\n')

f.close()
