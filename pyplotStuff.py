'''Ploting results'''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

rfname = 'Results.dat'
ptname = 'targetP.dat'

# Read results and assign to variables
data = np.loadtxt(rfname)
nx   = data[0]

lbound = 1
x = data[lbound:lbound+nx]
lbound += nx
pres = data[lbound:lbound+nx]
lbound += nx
rho  = data[lbound:lbound+nx]
lbound += nx
Mach = data[lbound:lbound+nx]
lbound += nx

xhalf= data[lbound:lbound+nx+1]
lbound += nx+1
S = data[lbound:lbound+nx+1]
lbound += nx+1

# Read target pressure
data = np.loadtxt(ptname)
nx = data[0]

targetp=data[1:nx+1]

pp=PdfPages('Figures.pdf')
# Plot Pressure vs Target Pressure
plt.figure(1)
plt.title('Pressure Distribution')
PressureCurve = plt.plot(x,pres,'-ob',markerfacecolor='None', markeredgecolor='b', label='Current')
TargetPCurve  = plt.plot(x,targetp,'-or',markerfacecolor='None', markeredgecolor='r', label='Target')
plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.5)
plt.legend(loc='upper right')
pp.savefig()

# Plot
plt.figure(2)
plt.title('Mach Distribution')
plt.plot(x,Mach,'-ob',markerfacecolor='None', markeredgecolor='b')
plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.5)
pp.savefig()

# Close PDF
pp.close()

