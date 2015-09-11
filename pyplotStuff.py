'''Ploting results'''

import matplotlib.pyplot as plt
import numpy as np

def readResults():
    rfname = 'Results.dat'
    ptname = 'target.dat'
    
    data = np.loadtxt(rfname)
    nx   = data[0]
    
    lbound = 1
    
    x    = data[lbound:lbound+nx]
    lbound += nx
    pres = data[lbound:lbound+nx]
    lbound += nx
    rho  = data[lbound:lbound+nx]
    lbound += nx
    Mach = data[lbound:lbound+nx]
    lbound += nx
    
    xhalf= data[lbound:lbound+nx+1]
    lbound += nx+1
    S    = data[lbound:lbound+nx+1]
    lbound += nx+1

def readTargetP():
    pass


readResults
figurePath='./Figures'
plt.figure(1)
plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.5)
plt.plot(x,pres,'-o',markerfacecolor='None', markeredgecolor='b')


plt.draw()

plt.show()
