"""
Python file to calculate masses of mesons.
Writes into file: Quark 1, Quark 2, Mass, Chi-Square.
"""
from project1 import groundstate

"""
Quark Masses

From literature:
m_u = 2.3 +0.7/-0.5 MeV
m_d = 4.8 +0.7/-0.3 MeV
m_s = 95 +/- 5 MeV
m_c = 1.275 +/- 0.025 GeV
m_b = 4.18 _ +/0.03 GeV

"""

m_u = m_d = 250. 	# MeV
m_s = 400. 			# MeV
m_c = 1.5 *10**3 	# MeV
m_b = 5.2 *10**3	# MeV

qm = [m_u,m_s,m_c,m_b]


"""
Meson Masses

ud = 139.57018 +/- 0.00035 MeV
us = 493.677 +/- 0.016 MeV
uc = 1869.62 +/- 0.15 MeV
ub = 5279.25 +/- 0.17 MeV

ss
sc = 1968.49 +/- 0.32 MeV
sb = 5366.77 +/- 0.24 MeV

cc = 2983.6 +/- 0.7 MeV
cb = 6277. +/- 6.0 MeV

bb = 9460.30 +/- 0.26 MeV
ds = 497.164 +/- 0.024 MeV
"""

Q1 = [0,0,0,0,1,1,1,2,2,3]
Q2 = [0,1,2,3,1,2,3,2,3,3]
Mesonlist = [0.139,0.493,1.869,5.279,0,1.968,5.366,2.983,6.277,9.460,0.497]


# Reduced mass of meson.
def MesonMass(n1,n2):
	return (qm[n1]*qm[n2])/(qm[n1]+qm[n2])

# Chi-square compares observed mass (m_o) and expected mass (m_e).
def ChiSquare(m_o,m_e):
	chisum = 0
	for i in range(1):
		chisum += (m_o - m_e)**2/m_e

	return chisum


"""
WRITING INTO PYTHON FILE
"""

with open('mesonmasses.txt', 'w') as f:
	
	f.write("Q1\tQ2\tMass (GeV)\t\tEnergy\n")

	for i in range(0,10):
		f.write("%d\t" % Q1[i])
		f.write("%d\t" % Q2[i])
		f.write("%f\t" % Mesonlist[i])
		f.write("%f\n" % groundstate())
		#f.write("%f\n" % MesonMass(i,j))

f.close()
