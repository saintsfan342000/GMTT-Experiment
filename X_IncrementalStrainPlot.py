import numpy as n
from numpy import array, nanmean, nanstd, sqrt
from numpy import hstack, vstack, dstack
import matplotlib.pyplot as p
p.style.use('mysty')
import figfun as f
import os

'''
Plots published values and ONLY Kelin's incremental values, pt-to-pt and averaged
'''

expts = n.array([1, 2, 3, 4, 5, 7, 9, 10])#, 12, 13, 14, 15])
#expts = n.array([2,12,13,3,14,15])
# [0]Expt No., [1]Expt Type, [2]Tube No., [3]Alpha, [4]Alpha-True, 
# [5]Mean Radius, [6]Thickness, [7]Eccentricity
# Type:  0=Radial, 1=Sigma-Tau Corner, 2=Tau-Sigma Corner
key = n.genfromtxt('../ExptSummary.dat', delimiter=',')

for k, x in enumerate(expts):

    # [0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = 0.657804 inch
    rot = n.genfromtxt('../TTGM-{}_FS19SS6/disp-rot.dat'.format(x), delimiter=',', usecols=(5))
    # [0]Stage, [1]Kelin VM, [2]Kelin H8, [3]Alternate VM, [4]Alternate H8
    e = n.genfromtxt('../TTGM-{}_FS19SS6/IncrementalStrain.dat'.format(x), delimiter=',')
    # # [0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx 
    # [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq
    sig, tau, dmean = n.genfromtxt('../TTGM-{}_FS19SS6/mean.dat'.format(x), 
            delimiter=',', usecols=(3,4,11), unpack=True)
    
    expinfo = key[ key[:,0] == x ].ravel()
    if n.isnan(expinfo[3]):
        masterlabel = '{:.0f}; $\\infty$ ; {:.0f}'.format(x,expinfo[2])
    else:
        masterlabel = '{:.0f}; {}; {:.0f}'.format(x,expinfo[3],expinfo[2])
    
    triax = ((sig/2)/n.sqrt(.75*(sig**2+4*tau**2)))[-1]

    # Fig0:  Eeq vs rotation
    if k == 0:
        fig0 = p.figure()
        ax0 = fig0.add_subplot(111)
    master, = ax0.plot(rot, dmean, label=masterlabel)
    line2, = ax0.plot(rot, e[:,1],'--', color=master.get_color())
    line3, = ax0.plot(rot, e[:,-2], ':', color=master.get_color())
    if x == expts[-1]:
        ax0.axis(xmin=0, ymin=0)
        ax0.set_title('Mean Strain vs Rotation')
        ax0.set_xlabel('$\\phi^\\circ$')
        ax0.set_ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$')
        lines = [master, line2, line3]
        leg02 = ax0.legend(lines, ['Normal', 'Increm/P2P', 'Increm/Avg'], loc='lower right')
        p.setp(leg02.get_lines(), color='k')
        f.ezlegend(ax0)
        ax0.add_artist(leg02)
        f.myax(ax0)

    # Fig1:  Failure Strain vs Triax
    if k == 0:
        fig1 = p.figure()
        ax1 = fig1.add_subplot(111)
    lines = []
    lines.append( ax1.plot(triax, dmean[-1], 'o', color='C0', mfc='none', label='Traditional')[0] )
    lines.append(ax1.plot(triax, e[-1,1], 's', color='C1', label='Increm/P2P')[0])
    lines.append(ax1.plot(triax, e[-1,-2], '^', color='C2', label='Increm/Avg')[0])
    if n.isnan(expinfo[3]):
        tx = ax1.text(triax,1.7,'\n$\\infty$',
                color=master.get_color(),ha='center',va='top',size=10)
    else:
        tx = ax1.text(triax,1.7,'\n{}'.format(expinfo[3]),
                color=master.get_color(),ha='center',va='top',size=10)
    if x == expts[-1]:
        ax1.axis([0,0.6,0,1.7])
        ax1.set_xlabel('$\\sigma_{\\mathsf{m}}/\\sigma_{\\mathsf{e}}$')
        ax1.set_ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$')
        leg1 = ax1.legend(lines, [i.get_label() for i in lines],
                loc='center left', bbox_to_anchor=(1.01,0.5),
                bbox_transform=ax1.transAxes,
                handlelength=0)
        [i.set_color(leg1.get_lines()[k].get_mec()) for k,i in enumerate(leg1.get_texts())]
        f.myax(ax1)

fig0.savefig('../8 - New Incremental eeq-rot.png', bbox_inches='tight', dpi=175)
fig1.savefig('../9 - New Incremental eeq-triax.png', bbox_inches='tight', dpi=175)
p.show('all')
