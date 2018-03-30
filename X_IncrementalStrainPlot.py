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
max_or_mean = 'max'
expts = n.array([ 31,   35, 20, 22,  24, 17, 27,   34,  32,   9,  30,  16, 15, 18])#, 12, 13, 14, 15])
alpha = n.array([.25, .375, .5, .5, .75,  1,  1, 1.25, 1.5, 2.0, 3.0, 3.5, 4.0, n.nan])
#expts = n.array([2,12,13,3,14,15])
# [0]Expt No., [1]Expt Type, [2]Tube No., [3]Alpha, [4]Alpha-True, 
# [5]Mean Radius, [6]Thickness, [7]Eccentricity
# Type:  0=Radial, 1=Sigma-Tau Corner, 2=Tau-Sigma Corner

for k, (x,alp) in enumerate(zip(expts,alpha)):

    # [0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = 0.657804 inch
    rot = n.genfromtxt('../TT2-{}_FS19SS6/disp-rot.dat'.format(x), delimiter=',', usecols=(5))
    # [0]Stage, [1]AvgF-All in Box-VM, [2]H8, [3]AvgF-Passing-VM, [4]H8, [5]Passing P2P-VM, [6]H8, [7]MaxPt-VM, [8]H8
    e = n.genfromtxt('../TT2-{}_FS19SS6/IncrementalStrain.dat'.format(x), delimiter=',')
    # # [0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx 
    # [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq
    sig, tau, dmean = n.genfromtxt('../TT2-{}_FS19SS6/mean.dat'.format(x), 
            delimiter=',', usecols=(3,4,11), unpack=True)
    dmax = n.genfromtxt('../TT2-{}_FS19SS6/MaxPt.dat'.format(x))[:,-1]
    
    if n.isnan(alp):
        masterlabel = '{:.0f}; $\\infty$ '.format(x)
    else:
        masterlabel = '{:.0f}; {}'.format(x,alp)
    
    triax = ((sig/2)/n.sqrt(.75*(sig**2+4*tau**2)))[-1]

    # Fig0:  Eeq vs rotation
    if k == 0:
        fig0 = p.figure()
        ax0 = fig0.add_subplot(111)
    
    if max_or_mean == 'max':
        data = dmax
        col = 7
    elif max_or_mean == 'mean':
        data = dmean
        col = 5

    master, = ax0.plot(rot, data, label=masterlabel)
    line2, = ax0.plot(rot, e[:,col],'--', color=master.get_color())
    line3 = 'spare'
    if max_or_mean == 'mean':
        line3, = ax0.plot(rot, e[:,3], ':', color=master.get_color())

    if x == expts[-1]:
        ax0.axis(xmin=0, ymin=0)
        ax0.set_title('M{} Strain vs Rotation'.format(max_or_mean[1:]))
        ax0.set_xlabel('$\\phi^\\circ$')
        ax0.set_ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$')
        lines = [master, line2, line3]
        labels = ['Normal', 'Increm/P2P', 'Increm/Avg']
        if max_or_mean == 'max':
            lines = lines[:2]
            labels = labels[:2]
        leg02 = ax0.legend(lines, labels, loc='lower right')
        p.setp(leg02.get_lines(), color='k')
        f.ezlegend(ax0)
        ax0.add_artist(leg02)
        f.myax(ax0)

    # Fig1:  Failure Strain vs Triax
    if k == 0:
        fig1 = p.figure()
        ax1 = fig1.add_subplot(111)
    lines = []
    lines.append( ax1.plot(triax, data[-1], 'o', color='C0', mfc='none', label='Traditional')[0] )
    lines.append(ax1.plot(triax, e[-1,col], 's', color='C1', label='Increm/P2P')[0])
    if max_or_mean == 'mean':
        lines.append(ax1.plot(triax, e[-1,3], '^', color='C2', label='Increm/Avg')[0])
    
    if n.isnan(alp):
        tx = ax1.text(triax,1.8,'\n$\\infty$',
                color=master.get_color(),ha='center',va='top',size=10)
    else:
        tx = ax1.text(triax,1.8,'{}'.format(alp),
                color=master.get_color(),ha='center',va='top',size=10)
    if x == expts[-1]:
        ax1.axis([0,0.6,0,1.8])
        ax1.set_xlabel('$\\sigma_{\\mathsf{m}}/\\sigma_{\\mathsf{e}}$')
        ax1.set_ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$')
        leg1 = ax1.legend(lines, [i.get_label() for i in lines],
                loc='center left', bbox_to_anchor=(1.01,0.5),
                bbox_transform=ax1.transAxes,
                handlelength=0)
        [i.set_color(leg1.get_lines()[k].get_mec()) for k,i in enumerate(leg1.get_texts())]
        f.eztext(ax1, 'M{} values'.format(max_or_mean[1:]), 'bl')
        f.myax(ax1)

fig0.savefig('../8 - Incremental {} eeq-rot.png'.format(max_or_mean), bbox_inches='tight', dpi=175)
fig1.savefig('../9 - Incremental {} eeq-triax.png'.format(max_or_mean), bbox_inches='tight', dpi=175)
p.show('all')
