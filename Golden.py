from __future__ import print_function, division
import numpy as np

from EnergyFunction import Calculate_Energy


"""
Golden Ratio Search
"""
state = 3

accuracy = 1*10**(-6)		# Required accuracy
z = (1+np.sqrt(5))/2		# Golden ratio

# Initial positions of the four points
x1 = 1. * 10**1.
x4 = 1. * 10**3.
x2 = x4 - (x4-x1)/z
x3 = x1 + (x4-x1)/z

def Golden(qm1,qm2,mm,s):
	global x1, x2, x3, x4

	f1 = Calculate_Energy(qm1,qm2,mm,x1)
	f2 = Calculate_Energy(qm1,qm2,mm,x2)
	f3 = Calculate_Energy(qm1,qm2,mm,x3)
	f4 = Calculate_Energy(qm1,qm2,mm,x4)

	# Main loop of search process
	while x4-x1>accuracy:
		if f2[s]<f3[s]:
			x4,f4[s] = x3,f3[s]
			x3,f3[s] = x2,f2[s]
			x2 = x4 - (x4-x1)/z
			f2[s] = Calculate_Energy(qm1,qm2,mm,x2)[s]
		else:
			x1,f1[s] = x2,f2[s]
			x2,f2[s] = x3,f3[s]
			x3 = x1 + (x4-x1)/z
			f3[s] = Calculate_Energy(qm1,qm2,mm,x3)[s]

	return (0.5*(x1+x4))

"""
WRITING INTO PYTHON FILE
"""

with open('mesonmasses.txt', 'r') as file1, open('alphalist.txt', 'w') as file2:

	#num_lines = sum(1 for line in f1)

	for line in file1:
		value = map(float,line.split())

		for st in range(state):
			file2.write("%d\t" % value[0])
			file2.write("%d\t" % value[1])
			file2.write("%f\t" % value[2])
			file2.write("%f\t" % value[3])
			file2.write("%f\t" % value[4])
			file2.write("%d\t" % st)			
			file2.write("%f\n" % Golden(value[2],value[3],value[4],st))

file1.close()
file2.close()