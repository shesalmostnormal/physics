from __future__ import print_function, division
from scipy.special import gamma
from numpy import pi, zeros, linalg as LA, exp, linspace, absolute as ab, sqrt, arange
from math import factorial
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend
import pylab as py

import proj12

""" 
Cython to calculate eigenvalues.
"""

#d = int(raw_input('Basis (d x d): '))
#n = int(raw_input('State (n): '))
#l = int(raw_input('Ang momentum (l): '))

K = 0.002 # Spring constant
M = 0.02 # Mass of particle
h_bar = 1

b0 = (M * K / h_bar) ** (1/4.)
b = b0 # Vary by percent
omega = (K/M) ** (0.5)

print('True alpha =',b)

d = int(raw_input('Basis (d x d): '))

x, y0, y1, y2 = proj12.main(M, omega, d, b, b0)

plt1, = py.plot(x,y0)
plt2, = py.plot(x,y1)
plt3, = py.plot(x,y2)

py.legend([plt1,plt2,plt3], ['E_0','E_1','E_10'])
py.xlabel('alpha')
py.ylabel('energy')
py.xlim(0,0.08)
py.ylim(0,50)
py.show()
