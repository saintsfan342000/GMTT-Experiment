import numpy as n
L = n.genfromtxt
fit = n.polyfit
abs=n.abs
from numpy import nanmax, nanmin
import os
import pandas as pd
import shutil
import matplotlib.pyplot as p
p.style.use('mysty-quad')
import figfun as f

makeplot = True

exporthead = 'Stage, AxSts(nom), ShSts(nom), Delta/L, Phi, eemax, eemean'
exporthead = exporthead.split(',')
llhead = ['Exp', 'Alpha']
[llhead.append(Z) for Z in exporthead]
#So I just need to open LV_DIC and InSGZone

expts = [1,2,3,4,5,7,9,10,11]

# [0]Expt No., [1]Material, [2]Tube No., [3]Alpha, [4]Alpha-True, [5]Mean Radius, [6]Thickness, [7]Eccentricity
key = n.genfromtxt('../ExptSummary.dat', delimiter=',')
key = key[ n.in1d(key[:,0],expts) ]
key = key[ key[:,3].argsort() ]
expts = key[:,0].astype(int)

fid = pd.ExcelWriter('../TTGM_SetData.xlsx')
LL = n.empty((len(expts),9)) #Store the LL properties

for k,X in enumerate(expts):
    print(X)
    
    if X in [8,11]:
        path = '../TTGM-{}_FS32SS8'.format(X)
    else:
        path = '../TTGM-{}_FS19SS6'.format(X)
        
    thickness = key[ k, 6 ]
    
    profstg = L(path+'/prof_stages.dat', delimiter=',').astype(int)
    
    # [0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = 0.686910 inch
    dr = L(path+'/disp-rot.dat', usecols=(0,2,3,4,5),
                          delimiter=',')
                          
    # Open file and read off the Lg
    with open(path+'/disp-rot.dat','r') as tempfid:
        line0 = tempfid.readline()
        tempfid.close()
    Lg = line0.split(' ')[-2]

    emean = L(path+'/mean.dat', delimiter=',', usecols=(11))
    emax = L(path+'/MaxPt.dat', delimiter=',', usecols=(10))

    d = n.c_[dr,emax,emean]
    d_stg = d[profstg].copy()
    LL[k,:2] = X, key[k,3]
    LL[k,2:] = d_stg[2]
    
    '''    
    LEprofs = L(path + '/StrainProfiles.dat', delimiter=',', skip_header=2)
    # I want col 0,2, 3,5, 6,8, 9,11, etc
    allcol = n.arange(LEprofs.shape[1])
    elim = n.arange(1,allcol[-2]+3,3)
    keepcol = allcol[ ~n.in1d(allcol, elim) ]
    LEprofs = LEprofs.take(keepcol, axis=1)
    LEprofs[:,0::2] *= (1/thickness)
    LEprofs = LEprofs[ n.abs(LEprofs[:,0])<=5 ]
    
    fig2 = p.figure()
    for k in range(int(LEprofs.shape[1]/2)):
        p.plot(LEprofs[:,2*k], LEprofs[:,2*k+1])
    
    '''
    
    if makeplot:
        if k == 0:
            fig, *ax = f.makequad()
        ax[0].plot(d[:,4],d[:,2])
        [ax[0].plot(d[k,4], d[k,2], 'o') for k in d_stg[:,0].astype(int)]
        ax[1].plot(d[:,3],d[:,1])
        [ax[1].plot(d[k,3], d[k,1], 'o') for k in d_stg[:,0].astype(int)]
        ax[2].plot(d[:,4], d[:,5])
        ax[2].plot(d[:,4], d[:,6], label=str(int(X)))
        ax[2].legend(loc='lower right')
        ax[3].plot(d[:,1],d[:,2])
        [ax[3].plot(d[k,1],d[k,2]) for k in d_stg[:,0].astype(int)]    
        if X == expts[-1]:
            ax[0].plot(LL[:,4+2],LL[:,2+2],'ko')
            ax[1].plot(LL[:,3+2],LL[:,1+2],'ko')
            ax[2].plot(LL[:,4+2],LL[:,6+2],'ko')
            ax[3].plot(LL[:,1+2],LL[:,2+2],'ko')
    
    exporthead[-1] += '   Lg={} in'.format(Lg)
    pd.DataFrame(d).to_excel( fid, sheet_name='Exp{} _ {}'.format(int(X),key[k,3]), index=False, header=exporthead)
    
    pd.DataFrame(d_stg).to_excel( fid, sheet_name='Exp{} _ {}'.format(int(X),key[k,3]), startcol = 10, index=False, header=exporthead)
    
    pd.DataFrame(LL).to_excel( fid, sheet_name='LL', index=False, header=llhead )
    
fid.save()
p.show('all')

