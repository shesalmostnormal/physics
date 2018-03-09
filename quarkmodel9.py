from __future__ import print_function, division
import numpy as np
from numpy.linalg import eigvalsh as GetEigen
from scipy.special import gamma, roots_laguerre as GetRoots, roots_genlaguerre as GetGenRoots, eval_genlaguerre as LaguerreL

import time

"""
Jan. 25, 2018

For heavy quarkonium states the bound state energy is much lower
than the mass of the quarks (hence non-relativistic) so calculation
works best with heavy quarks (eg. less accurate for pi mesons). For
relativistic, must include state of gluons and use equation different
from Schrofdinger.

"""

#d = int(raw_input('Basis (d x d matrix): ')) #Basis

d = 30

state = 3

start = time.time()

points = 30

size = range(0,d)
dimensions = (len(size),len(size))

l = 0

coupling = (4./3.)*0.2 # Fine structure
sigma = 1000. *197.#* 10**6. # String tension



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


def Calculate_Energy(m,qm,alpha):

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
			
			Hmatrix[i][j] =  (norm(i)*norm(j))*(0.5*(1/m) * p_sum[i][j] \
				- coupling * r_sum[i][j] + sigma * l_sum[i][j] + qm )


	EigenList = GetEigen(Hmatrix)
	
	# Sorting for eigenvalues
	#idx = u.argsort()[::1]   
	#u = u[idx]
	
	u = []
	for v in range(state):
		u.append(EigenList[v])
	
	return(u)


"""
Golden Ratio Search
"""
accuracy = 1*10**(-6)		# Required accuracy
z = (1+np.sqrt(5))/2		# Golden ratio

# Initial positions of the four points
x1 = 1. * 10**1.
x4 = 1. * 10**3.
x2 = x4 - (x4-x1)/z
x3 = x1 + (x4-x1)/z

def Golden(m,qm,s):
	global x1, x2, x3, x4

	f1 = Calculate_Energy(m,qm,x1)
	f2 = Calculate_Energy(m,qm,x2)
	f3 = Calculate_Energy(m,qm,x3)
	f4 = Calculate_Energy(m,qm,x4)

	# Main loop of search process
	while x4-x1>accuracy:
		if f2[s]<f3[s]:
			x4,f4[s] = x3,f3[s]
			x3,f3[s] = x2,f2[s]
			x2 = x4 - (x4-x1)/z
			f2[s] = Calculate_Energy(m,qm,x2)[s]
		else:
			x1,f1[s] = x2,f2[s]
			x2,f2[s] = x3,f3[s]
			x3 = x1 + (x4-x1)/z
			f3[s] = Calculate_Energy(m,qm,x3)[s]

	return (0.5*(x1+x4))

"""
Mesons
"""

m_u = m_d = 250. 	# MeV
m_s = 400. 			# MeV
m_c = 1.5 *10**3 	# MeV
m_b = 5.2 *10**3	# MeV

m_ud = m_u*m_d/(m_u+m_d)
m_uu = m_u*m_u/(m_u+m_u)
m_us = m_u*m_s/(m_u+m_s)
m_uc = m_u*m_c/(m_u+m_c)
m_ub = m_u*m_b/(m_u+m_b)

m_sc = m_s*m_c/(m_s+m_c) 
m_sb = m_s*m_b/(m_s+m_b)
m_ss = m_s*m_s/(m_s+m_s)

m_bc = m_b*m_c/(m_b+m_c)
m_bb = m_b*m_b/(m_b+m_b)
 
m_cc = m_c*m_c/(m_c+m_c)

mesons = ['ud', 'uu', 'us', 'uc', 'ub', 'sc', 'sb', 'ss', 'bc', 'bb', 'cc']
quark = Point(250., 400, 1.5*10**3, 5.2*10**3)
mass = [m_ud ,m_uu, m_us, m_uc, m_ub, m_sc, m_sb, m_ss, m_bc, m_bb, m_cc]
qmass = [m_u+m_d ,m_u+m_u, m_u+m_s, m_u+m_c, m_u+m_b, m_s+m_c, m_s+m_b, m_s+m_s, m_b+m_c, m_b+m_b, m_c+m_c]


mesonlist = []

for j in range(len(mass)):
	alphalist = []
	for i in range(state):
		alphalist.append(Golden(mass[j],qmass[j],i))
	mesonlist.append(Calculate_Energy(mass[j],qmass[j],alphalist[i]))

stop = time.time()

"""
WRITING INTO PYTHON FILE
"""

import csv

with open('mesonsplusmasses.csv', 'wb') as f:

	writer = csv.writer(f, delimiter='\t')

	writer.writerow(['Meson','	E_0','	E_1','	E_2'])

	record_list = [ list(item) for item in list(zip(mesons, qmass, mesonlist)) ]
	writer.writerows(record_list)


#print('Time of calculation:', stop-start)