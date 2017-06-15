import numpy as n
L = n.genfromtxt
abs=n.abs
import matplotlib.pyplot as p
import os
from sys import argv
from scipy.signal import detrend
import figfun as ff
from pandas import read_excel

try:
    worthless, expt, FS, SS, path = argv
    expt = int(expt)
    FS = int(FS)
    SS = int(SS)
    savefigs = True
except ValueError:
    expt = 9
    FS = 19
    SS = 6
    path = '../TTGM-{}_FS{}SS{}'.format(expt,FS,SS)
    savefigs = True


try:
   os.chdir(path)
except FileNotFoundError:
    os.chdir('TYPE PATH HERE')

key = n.genfromtxt('{}/../ExptSummary.dat'.format(path), delimiter=',')
alpha, Rm,thickness, tube, material = key[ key[:,0] == expt, [3,5,6,2,1]].flatten()

if n.isnan(alpha):
    alpha = '$\\infty$'
    


#########
#Max.dat
#   '[0]Stage [1]Time [2]AxSts [3]ShSts [4]NEx [5]NEy [6]Gamma [7]F11-1 [8]F22-1 [9]atan(F12/F22) [10]epeq [11]AramX [12]AramY'
#mean.dat
#   [0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
#like_Scott
#   [0]Stage [1]Time [2]SizeAveragingZone(in) [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
# profStgs
#   'Stages at which profiles were generated'
# profUr
#   'First Row Stage number.  Second row begin data.\n[0]Ycoord [1:] Stage Ur/Ro'
#MaxPt.dat
#   [0]Stage [1]Time [2]AxSts [3]ShSts [4]NEx [5]NEy [6]Gamma [7]F11-1 [8]F22-1 [9]atan(F12/F22) [10]epeq'

STF = L('STF.dat',delimiter=',')
dmax = L('max.dat',delimiter=',')
dmaxPt = L('MaxPt.dat',delimiter=',')
dmean = L('mean.dat',delimiter=',')
dscot = L('like_Scott.dat',delimiter=',')
DR = L('disp-rot.dat',delimiter=',')
#'[0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = {:.6f} inch'
profStg = L('prof_stages.dat',delimiter=',').astype(int)
profLEp = L('StrainProfiles.dat',delimiter=',')[1:]
profUr = L('RadialContraction.dat',delimiter=',')[1:]

##################################################
# Figure 1 - AxSts-Delta and ShearSts-Rot
##################################################
p.style.use('mysty-sub')
fig1 = p.figure(1,facecolor='w',figsize=(8,12))
p.subplot(2,1,1)

for j in profStg:
    p.plot(DR[int(j),4],DR[int(j),2],'o',markersize=7,mec='none')
p.plot(DR[:,4],DR[:,2],'b',zorder=0)
p.xlabel('$\\delta/\\mathsf{L}$')
p.ylabel('$\\Sigma$\n$(\\mathsf{ksi})$')
titlestring = 'TTGM-{:.0f}, $\\alpha$ = {}.  FS{:.0f}SS{:.0f}. Tube {:.0f}-{:.0f}'.format(expt,alpha,FS,SS,material,tube)
p.title(titlestring,size=18)
p.gcf().gca().set_ylim([0,1.2*n.max(STF[:,2])])
p.gcf().gca().set_xlim(left=0)

ff.myax(p.gcf(),ff.ksi2Mpa,'$\\Sigma$\n$(\\mathsf{MPa})$')

p.subplot(2,1,2)
for j in profStg:
    p.plot(DR[j,5],DR[j,3],'o',markersize=7,mec='none')
p.plot(DR[:,5],DR[:,3],'b',zorder=0)
p.xlabel('$\\phi^{\\circ}$')
p.ylabel('$\\mathcal{T}$\n$(\\mathsf{ksi})$')
p.gcf().gca().set_ylim([0,1.2*n.max(STF[:,3])])
p.gcf().gca().set_xlim(left=0)

ff.myax(p.gcf(),ff.ksi2Mpa,'$\\mathcal{T}$\n$(\\mathsf{MPa})$')

##################################################
# Figure 2 - Epsilon-Gamma
##################################################
p.style.use('mysty')
fig2 = p.figure(2,facecolor='w')
p.plot(abs(dmax[:,6]),dmax[:,5],'ro',mec='r',ms=4,label='Max',alpha=0.5)
p.plot(abs(dmean[:,7]),dmean[:,6],'o',mec='b',mfc='b',ms=4,label='Mean',alpha=0.5)
l1 = p.plot(abs(dscot[:,7]),dscot[:,6],'ko',mec='k',ms=4,label='Macro',alpha=0.5)[0]
p.plot(abs(dmax[:,9]),dmax[:,8],'s',mfc='none', mec='r',ms=4,alpha=0.5)
p.plot(abs(dmean[:,10]),dmean[:,9],'s',mec='b',mfc='none',ms=4,alpha=0.5)
l2 = p.plot(abs(dscot[:,10]),dscot[:,9],'s', mfc='none',mec='k',ms=4,alpha=0.5)[0]
L1 = p.legend([l2,l1],["Haltom 2013","Aramis"],numpoints=1,handletextpad=.000,title="Stn. Defn.",loc='upper left',bbox_to_anchor=(.99,1),borderpad=0,borderaxespad=0,frameon=False)
p.setp(L1.get_title(),fontsize=L1.get_texts()[1].get_fontsize())
L2 = p.legend(loc='center left',bbox_to_anchor=(.99,.5),numpoints=1,handletextpad=.000,borderpad=0,borderaxespad=0,frameon=False)
p.gca().add_artist(L1)
p.xlabel('$\\gamma$')
p.ylabel('$\\epsilon_{\\mathsf{y}}$')
p.title(titlestring,size=18)
p.gcf().gca().set_ylim(bottom=0)
p.gcf().gca().set_xlim(left=0)

ff.myax(p.gcf())

##################################################
# Figure 3 - Strain Profile thru Max Pt
##################################################
p.style.use('mysty')
fig3 = p.figure(3,facecolor='w',figsize=(12,6) )
p.gcf().add_axes([.12,.12,.8,.78])
for i in range(len(profStg)):
    p.plot(profLEp[:,3*i]/thickness,profLEp[:,3*i+2],lw=1.5)
#p.gcf().gca().set_ylim(bottom=0)
p.title(titlestring,size=18)
p.xlabel('y$_{\\mathsf{o}}$/t$_{\\mathsf{o}}$')
p.ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}}$')
p.gca().set_xlim([-8,8])

ff.myax(p.gcf(),TW=.0025,HW=.3,HL=.05,OH=.2)

##################################################
# Figure 4 - Radial contraction profile thru x = 0
##################################################
profUr = profUr[~n.any(n.isnan(profUr),axis=1), :]  # Detrend the data!
profUr[:,1:] = detrend(profUr[:,1:],axis=0)
profUr[:,1:] -= n.nanmax( profUr[:,1:], axis=0 )
p.style.use('mysty')
fig4 = p.figure(4,facecolor='w',figsize=(12,6) )
p.gcf().add_axes([.12,.12,.8,.78])    
p.plot(2*profUr[:,0]/0.62,profUr[:,1:],lw=1.5)
#p.gcf().gca().set_ylim(bottom=0)
p.title(titlestring,size=18)
p.xlabel('2y$_{\\mathsf{o}}$/L$_{\\mathsf{g}}$')
p.ylabel('$\\frac{\\mathsf{u}_{\\mathsf{r}}}{\\mathsf{R}_{\\mathsf{o}}}$')
p.gca().set_xlim([-1,1])
p.grid(True,which='both',axis='y',linestyle='--',alpha=0.5)

ff.myax(p.gcf(),TW=.0025,HW=.3,HL=.05,OH=.2)

##################################################
# Figure 5 - LastStgMax vs each Stage Max
##################################################    
p.style.use('mysty')
fig5 = p.figure(5)
p.plot(abs(dmax[:,6]),dmax[:,5],'ro',label='Each Stage Max',mec='r',ms=4,alpha=0.5)
p.plot(abs(dmax[:,9]),dmax[:,8],'s',mfc='none', mec='r',ms=4,alpha=0.5)
l1, = p.plot(abs(dmaxPt[:,6]),dmaxPt[:,5],'bo',label='Last Stage Max',mec='b',ms=4,alpha=0.5)
l2, = p.plot(abs(dmaxPt[:,9]),dmaxPt[:,8],'s',mfc='none', mec='b',ms=4,alpha=0.5)
p.xlabel('$\\gamma$')
p.ylabel('$\\epsilon_\\mathsf{y}$')
p.title(titlestring + '\nComparing Last Stg Max to Each Stage Max',size=14)
p.gcf().gca().set_ylim(bottom=0)
p.gcf().gca().set_xlim(left=0)
L1 = p.legend([l1,l2],["Aramis","Haltom 2013"],numpoints=1,handletextpad=.000,title="Stn. Defn.",loc='upper left',bbox_to_anchor=(.99,1),borderpad=0,borderaxespad=0,frameon=False)
p.setp(L1.get_title(),fontsize=L1.get_texts()[0].get_fontsize())
L2 = p.legend(loc='center left',bbox_to_anchor=(.99,.5),numpoints=1,handletextpad=.000,borderpad=0,borderaxespad=0,frameon=False)
p.gca().add_artist(L1)
ff.myax(p.gcf())

##################################################
# Figure 6 - Machine force-disp, torque-rot respnses
##################################################  
# [0]Stg [1]Time [2]AxSts [3]ShSts [4]AxForce [5]Torque [6]MTS Disp [7]MTS Rot
STF[:,[6,7]]-=STF[0,[6,7]]
p.style.use('mysty-sub')
fig6, ax61, ax62 = ff.make21()
ax61.grid(axis='x',alpha=0.5)
ax61.plot(STF[:,7],STF[:,5],'b')
loc = profStg[2]
ax61.plot(STF[loc,7],STF[loc,5],'r^',mec='r',mfc='r',ms=8)
ax61.set_xlabel('$\\phi^{\\circ}$')
ax61.set_ylabel('T\n(lb-in)')
ax61.set_title(titlestring)
ax61.axis(xmin=0,ymin=0)
LLstr = '$_{\\mathsf{LL}}$'
Fastr = '$_{\\mathsf{f}}$'
annot = ''.join(['T{:s} = {:.0f} lb-in\n'.format(LLstr,STF[loc,5]),
        '$\\phi${:s} = {:.1f}$^{{\\circ}}$\n'.format(LLstr,STF[loc,7]),
        'T{:s} = {:.0f} lb-in\n'.format(Fastr,STF[-1,5]),
        '$\\phi${:s} = {:.1f}$^{{\\circ}}$\n'.format(Fastr,STF[-1,7]),
        'Time: {:.0f} sec'.format(STF[-1,1])])
ff.eztext(ax61,annot,'br')
ff.myax(ax61)
ax62.grid(axis='x',alpha=0.5)
ax62.plot(STF[:,6],STF[:,4],'b')
ax62.plot(STF[loc,6],STF[loc,4],'r^',mec='r',mfc='r',ms=8)
ax62.set_xlabel('$\\delta$',fontsize=22)
ax62.set_ylabel('F\n(lb)')
ax62.axis(xmin=0,ymin=0)
annot = ''.join(['F{:s} = {:.0f} lb\n'.format(LLstr,STF[loc,4]),
        '$\\delta${:s} = {:.3f} in\n'.format(LLstr,STF[loc,6]),
        'F{:s} = {:.0f} lb\n'.format(Fastr,STF[-1,4]),
        '$\\delta${:s} = {:.3f} in\n'.format(Fastr,STF[-1,6]),
        'Time: {:.0f} sec'.format(STF[-1,1])])        
ff.eztext(ax62,annot,'br')
ff.myax(ax62)
atrue = n.polyfit(STF[:,3],STF[:,2],1)[0]
p.text(.5,.48,'$\\alpha_\\mathsf{{true}}$ = {:.2f}'.format(atrue),transform=fig6.transFigure,ha='center',va='center')

if savefigs:
    fig5.savefig('5 - MaxComparison.png',dpi=125)
    fig4.savefig('4 - RadialContraction.png',dpi=125)
    fig3.savefig('3 - StrainProfile.png',dpi=125)
    fig2.savefig('2 - StrainPath.png',dpi=125)
    fig1.savefig('1 - Sts-Delta-Rot.png',dpi=125)
    fig6.savefig('BinderFig1.pdf',bbox_inches='tight')
    p.close('all')
else:
    p.show('all')

