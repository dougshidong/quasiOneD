#!/usr/bin/python
'''Ploting results'''
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

f = open("input.in", "r")
header = f.readline()
header = f.readline()
header = f.readline()
header = f.readline()
line   = f.readline()
line   = line.strip()
columns = line.split()
nx = int(columns[0])

dx = 1.0 / nx

x = np.linspace(dx/2.0, 1.0 - dx/2.0, nx)

## Plot Channel
#plt.figure()
#plt.title('Channel Shape')
#cshape1 = plt.plot(xhalf,S,'-ob',markerfacecolor='None', markeredgecolor='b')
#plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.5)
#pp.savefig()

def readPMachConv(filename, nx):
    flfname = './Results/'+filename+'Flow.dat'
    fcfname = './Results/'+filename+'FlowConv.dat'
    ftfname = './Results/'+filename+'FlowTime.dat'

    # Read flow results
    data = np.loadtxt(flfname)
    lbound = 0
    pres = data[lbound:lbound+nx]
    lbound += nx
    rho  = data[lbound:lbound+nx]
    lbound += nx
    Mach = data[lbound:lbound+nx]
    lbound += nx
    # Read convergence information
    data = np.loadtxt(fcfname)
    lbound = 0
    conv = data[lbound:len(data)]
    # Read time convergence
    data = np.loadtxt(ftfname)
    lbound = 0
    time = data[lbound:len(data)]

    return pres, Mach, conv, time


def plotPressure(pressure, names, pp):
    # Plot Pressure Distribution
    plt.figure(figsize=(8,4))
    plt.title('Pressure Distribution')

    for n, name in enumerate(names):
        nx = len(pressure[n][:])
        dx = 1.0 / nx
        x = np.linspace(dx/2.0, 1.0 - dx/2.0, nx)

        plt.plot(x, pressure[n][:], color = col[n],
                 marker = symb[n], mec = col[n], mfc = 'None', ms = 4,
                 label=name)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$p/p_{t_{in}}$')
    plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.2)
    plt.legend(loc='lower left')
    plt.tight_layout()
    pp.savefig(bbx_inches='tight')

def plotMach(Mach, names, pp):
    plt.figure(figsize=(8,4))
    plt.title('Mach Distribution')
    for n, name in enumerate(names):
        nx = len(Mach[n][:])
        dx = 1.0 / nx
        x = np.linspace(dx/2.0, 1.0 - dx/2.0, nx)
        plt.plot(x, Mach[n][:], color = col[n],
                 marker = symb[n], mec = col[n], mfc = 'None', ms = 6,
                 label=name)
    plt.xlabel(r'$x$')
    plt.ylabel(r'Mach')
    plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.2)
    plt.legend(loc='upper left')
    plt.tight_layout()
    pp.savefig(bbx_inches='tight')

def plotConv(conv, names, pp):
    plt.figure(figsize=(8,4))
    plt.title('Residual Convergence')
    step = 1
    for n, name in enumerate(names):
        plt.semilogy(range(1, len(conv[n][:]), step),
                     conv[n][0:-1:step],
                     color = col[n],
                     marker = symb[n], mec = col[n], mfc = 'None', ms = 0,
                     label=name)

    plt.xlabel(r'Iterations')
    plt.ylabel(r'Density Residual')
    plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.2)
    plt.legend(loc='upper right')
    plt.tight_layout()
    pp.savefig(bbx_inches='tight')

def plotTime(time, conv, names, pp):
    plt.figure(figsize=(8,4))
    plt.title('Residual Convergence')
    step = 1
    for n, name in enumerate(names):
        plt.semilogy(time[n][0:-1:step],
                     conv[n][0:-1:step],
                     color = col[n],
                     marker = symb[n], mec = col[n], mfc = 'None', ms = 0,
                     label=name)

    plt.xlabel(r'Times [s]')
    plt.ylabel(r'Density Residual')
    plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.2)
    plt.legend(loc='upper right')
    plt.tight_layout()
    pp.savefig(bbx_inches='tight')

def question1():
    filenames = ['q1']
    for n, name in enumerate(filenames):
        p, m, c, t = readPMachConv(name, nx)
        pres = [p]
        mach = [m]
        conv = [c]

    curvenames = ['Test']
    pdfName = './report/figures/q1p.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q1m.pdf'
    pp=PdfPages(pdfName)
    plotMach(mach, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q1c.pdf'
    pp=PdfPages(pdfName)
    plotConv(conv, curvenames, pp)
    pp.close()

def question2():
    filenames = ['q2_76', 'q2_72', 'q2_68', 'q2_60']
    for n, name in enumerate(filenames):
        p, m, c, t = readPMachConv(name, nx)
        if n == 0:
            pres = [p]
            mach = [m]
            conv = [c]
        else:
            pres.append(p)
            mach.append(m)
            conv.append(c)

    curvenames = ['$p_{exit} = 0.76p_{t_{in}}$',
                  '$p_{exit} = 0.72p_{t_{in}}$',
                  '$p_{exit} = 0.68p_{t_{in}}$',
                  '$p_{exit} = 0.60p_{t_{in}}$']

    pdfName = './report/figures/q2p.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q2m.pdf'
    pp=PdfPages(pdfName)
    plotMach(mach, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q2c.pdf'
    pp=PdfPages(pdfName)
    plotConv(conv, curvenames, pp)
    pp.close()

def question3():
    filenames = ['q3_25', 'q3_50', 'q3_100', 'q3_200']
    nxs = [25, 50, 100, 200]
    for n, name in enumerate(filenames):
        p, m, c, t = readPMachConv(name, nxs[n])
        if n == 0:
            pres = [p]
            mach = [m]
            conv = [c]
        else:
            pres.append(p)
            mach.append(m)
            conv.append(c)

    curvenames = ['$nx = 25$',
                  '$nx = 50$',
                  '$nx = 100$',
                  '$nx = 200$']

    pdfName = 'Figures.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    plotMach(mach, curvenames, pp)
    plotConv(conv, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q3p.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q3m.pdf'
    pp=PdfPages(pdfName)
    plotMach(mach, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q3c.pdf'
    pp=PdfPages(pdfName)
    plotConv(conv, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q3loss.pdf'
    pp=PdfPages(pdfName)
    plt.figure(figsize=(8,4))
    plt.title('Total Pressure Loss vs Grid Size')
    plt.plot([25, 50, 100, 200], [0.0367182,0.0379852,0.0385511,0.0388222],
             color = col[0],
             marker = symb[0], mec = col[0], mfc = 'None', ms = 6,
             label=name)
    plt.xlabel(r'Grid Size')
    plt.ylabel(r'Total Pressure Loss')
    plt.grid(b=True, which='major', color='black', linestyle='-',alpha=0.2)
    plt.tight_layout()
    pp.savefig(bbx_inches='tight')
    pp.close()

def question4():
    filenames = ['q4_SD', 'q4_SW', 'q4_CMSW', 'q4_Roe']
    for n, name in enumerate(filenames):
        p, m, c, t = readPMachConv(name, nx)
        if n == 0:
            pres = [p]
            mach = [m]
            conv = [c]
            time = [t]
        else:
            pres.append(p)
            mach.append(m)
            conv.append(c)
            time.append(t)

    curvenames = ['SD',
                  'SW',
                  'CMSW',
                  'Roe']

    pdfName = 'Figures.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    plotMach(mach, curvenames, pp)
    plotConv(conv, curvenames, pp)
    plotTime(time, conv, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q4p.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q4m.pdf'
    pp=PdfPages(pdfName)
    plotMach(mach, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q4c.pdf'
    pp=PdfPages(pdfName)
    plotConv(conv, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q4t.pdf'
    pp=PdfPages(pdfName)
    plotTime(time, conv, curvenames, pp)
    pp.close()

def question5():
    filenames = ['q5_Euler', 'q5_jrk4', 'q5_eulerI']
    for n, name in enumerate(filenames):
        p, m, c, t = readPMachConv(name, nx)
        if n == 0:
            pres = [p]
            mach = [m]
            conv = [c]
            time = [t]
        else:
            pres.append(p)
            mach.append(m)
            conv.append(c)
            time.append(t)

    curvenames = ['EE',
                  'JRK4',
                  'EI']

    pdfName = 'Figures.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    plotMach(mach, curvenames, pp)
    plotConv(conv, curvenames, pp)
    plotTime(time, conv, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q5p.pdf'
    pp=PdfPages(pdfName)
    plotPressure(pres, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q5m.pdf'
    pp=PdfPages(pdfName)
    plotMach(mach, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q5c.pdf'
    pp=PdfPages(pdfName)
    plotConv(conv, curvenames, pp)
    pp.close()

    pdfName = './report/figures/q5t.pdf'
    pp=PdfPages(pdfName)
    plotTime(time, conv, curvenames, pp)
    pp.close()

#print 'Plotting Question 1. Don\'t forget to remove legend'
#question1()
#print 'Done Q1'

#print 'Plotting Question 2'
#question2()
#print 'Done Q2'

#print 'Plotting Question 3'
#question3()
#print 'Done Q3'

#print 'Plotting Question 4'
#question4()
#print 'Done Q4'

print 'Plotting Question 5'
question5()
print 'Done Q5'

