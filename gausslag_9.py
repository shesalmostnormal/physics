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

N = 20
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


u,v = lg(20)

point = [0.237787, 0.704138, 1.40675, 2.3497, 3.53851, 4.98045, 6.68474, \
	8.66296, 10.9295, 13.5024, 16.4041, 19.6633, 23.3165, 27.4122, 32.0161,\
	37.2213, 43.1679, 50.0868, 58.4192, 69.2956]
weight = []

s = 0
for s in range(0,20):
	weight.extend([gamma(N+1.5+1)*point[s]/(gamma(N+1.)*(N+1.)**2*\
		(LaguerreL(21, 1.5, point[s]))**2)])

