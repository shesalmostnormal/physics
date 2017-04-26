from __future__ import print_function, division
from scipy.special import gamma
from numpy.polynomial.laguerre import laggauss as lg
from sympy import *
from sympy.solvers import solve
from numpy import *
from math import sqrt

"""
Gauss-Laguerre formula GENERATE SAMPLE POINTS AND WEIGHTS WORKING!!!!!!!!!!!!!!!

"""

N = 4
alpha = 1.5

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


R = LaguerreL(N,alpha,x)
Laguerreprime = R.diff(x)

point = solve(x * Laguerreprime + ((N+alpha) * LaguerreL(N-1,alpha,x)), x)

pointlist = []
for i in range(len(point)):
	pointlist += [re(point[i])]

weightlist = []
for j in range(0,N):
	weightlist.extend([gamma(N+1.5+1)*pointlist[j]/(gamma(N+1.)*(N+1.)**2*\
		(LaguerreL(N+1, 1.5, pointlist[j]))**2)])

print(pointlist,weightlist)