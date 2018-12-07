#!/usr/local/bin/python3
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
col=('b', 'g', 'r', 'c', 'm', 'y', 'k')
symb = ('o', 'v', '<', 's', 'p', '*', 'h', 'H', 'D', 'd')

case_name = sys.argv[1]
filename = 'flow_convergence.dat'
full_path = '../Results/'+case_name+'_'+filename
flow_conv = np.loadtxt(full_path, dtype={'names': ('text', 'time', 'residual'),
									'formats': ('S', 'f4', 'f4')})

filename = 'gradient_convergence.dat'
full_path = '../Results/'+case_name+'_'+filename
opt_conv = np.loadtxt(full_path, dtype={'names': ('text', 'time', 'gradient'),
									'formats': ('S', 'f4', 'f4')})


pdfName = './'+case_name+'_convergence.pdf'
pp=PdfPages(pdfName)

plt.figure()
print(flow_conv['time'])
print(flow_conv['residual'])
plt.semilogy(flow_conv['time'],flow_conv['residual'],marker='o',ms=1, lw=1, label=r"flow\_residual")
plt.semilogy(opt_conv['time'],opt_conv['gradient'],marker='o',ms=2, label=r"gradient\_norm")
plt.tight_layout()
plt.legend(loc='upper right')
plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.2)
#plt.axis([10, 35, 1, 1e-12])
plt.xlabel(r'CPU time $[s]$')

pp.savefig(bbx_inches='tight')

pp.close()
