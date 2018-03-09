from __future__ import print_function, division
import numpy as np
from numpy.linalg import eigvalsh as GetEigen
from scipy.special import gamma, roots_laguerre as GetRoots, roots_genlaguerre as GetGenRoots, eval_genlaguerre as LaguerreL


"""
Energy Eigenvalues written into file.
Golden.py generates an alphalist.txt
Quark and meson masses from mesonmasses.txt
Call Calculate_Energy for eigenvalues.
"""

from EnergyFunction import Calculate_Energy

"""
WRITING INTO PYTHON FILE
"""

with open('alphalist.txt', 'r') as file1, open('eigenlist.txt', 'w') as file2:

	#num_lines = sum(1 for line in f1)

	for line in file1:
		value = map(float,line.split())
		
		for i in range(3):	
			file2.write("%d " % value[0])
			file2.write("%d " % value[1])
			file2.write("%f\t" % value[2])
			file2.write("%f\t" % value[3])
			
			file2.write("%f\n" % Calculate_Energy(value[2],value[3],value[4],value[6])[i])
		
file1.close()
file2.close()