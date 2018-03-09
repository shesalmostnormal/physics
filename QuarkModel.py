from __future__ import print_function, division
import numpy as np
from numpy.linalg import eigvalsh as GetEigen
from scipy.special import gamma, roots_laguerre as GetRoots, roots_genlaguerre as GetGenRoots, eval_genlaguerre as LaguerreL


"""
Read files for output.
"""

q1 = int(raw_input('Quark 1 (0=u/d, 1=s, 2=c, 3=d): '))
q2 = int(raw_input('Quark 2 (0=u/d, 1=s, 2=c, 3=d): '))


"""
READING FROM PYTHON FILE
"""

with open('eigenlist.txt', 'r') as f:

	for line in f:
		value = map(str,line.split())
	
	if (q1 == value[0] and q2 == value[1]):
		print(value[5])
		
f.close()
