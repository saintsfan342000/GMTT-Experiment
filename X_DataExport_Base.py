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

exporthead = 'Stage, AxSts(nom), ShSts(nom), Delta/L, Phi, eemax, eemean'
exporthead = exporthead.split(',')
#So I just need to open LV_DIC and InSGZone

expts = ['3_FS19SS6']

fid = pd.ExcelWriter('../TTGM-{}_ExptData.xlsx'.format(expts[0][0]))
misc = n.empty(len(expts)) #Store the misc properties

for k,X in enumerate(expts):
    print(X)
    
    path = '../TTGM-{}'.format(X)
    
    key = n.genfromtxt('../ExptSummary.dat', delimiter=',')
    thickness = key[ key[:,0] == int(X[0]), :].ravel()[-2]
    
    profstg = L(path+'/prof_stages.dat', delimiter=',').astype(int)
    
    stage,time,sig,tau,delta,phi = L(path+'/disp-rot.dat',
                          delimiter=',', unpack=True)
                          
    # Open file and read off the Lg
    with open(path+'/disp-rot.dat','r') as tempfid:
        line0 = tempfid.readline()
        tempfid.close()
    Lg = line0.split(' ')[-2]

    emean = L(path+'/mean.dat', delimiter=',', usecols=(11))
    emax = L(path+'/MaxPt.dat', delimiter=',', usecols=(10))

        
    LEprofs = L(path + '/StrainProfiles.dat', delimiter=',', skip_header=2)
    # I want col 0,2, 3,5, 6,8, 9,11, etc
    allcol = n.arange(LEprofs.shape[1])
    elim = n.arange(1,allcol[-2]+3,3)
    keepcol = allcol[ ~n.in1d(allcol, elim) ]
    LEprofs = LEprofs.take(keepcol, axis=1)
    LEprofs[:,0::2] *= (1/thickness)
    LEprofs = LEprofs[ n.abs(LEprofs[:,0])<=5 ]
    
    d = n.c_[stage,sig,tau,delta,phi,emax,emean]
    d_stg = d[profstg].copy()

    fig, *ax = f.makequad()
    ax[0].plot(d[:,4],d[:,2])
    [ax[0].plot(d[k,4], d[k,2], 'o') for k in d_stg[:,0].astype(int)]
    ax[1].plot(d[:,3],d[:,1])
    [ax[1].plot(d[k,3], d[k,1], 'o') for k in d_stg[:,0].astype(int)]
    ax[2].plot(d[:,4], d[:,5], label='max')
    ax[2].plot(d[:,4], d[:,6], label='mean')
    ax[2].legend(loc='lower right')
    ax[3].plot(d[:,1],d[:,2])
    [ax[3].plot(d[k,1],d[k,2]) for k in d_stg[:,0].astype(int)]
    
    fig2 = p.figure()
    for k in range(int(LEprofs.shape[1]/2)):
        p.plot(LEprofs[:,2*k], LEprofs[:,2*k+1])
        
    exporthead[-1] += '   Lg={} in'.format(Lg)
    pd.DataFrame(d).to_excel( fid, sheet_name='Exp{}'.format(X.split('_')[0]), index=False, header=exporthead)
    
    pd.DataFrame(d_stg).to_excel( fid, sheet_name='Exp{}_stages'.format(X.split('_')[0]), index=False, header=exporthead)
    
    header=['Stg{} yo/to,Stg{} ee'.format(k,k) for k in profstg]
    header=','.join(header)
    header=header.split(',')
    pd.DataFrame(LEprofs).to_excel( fid, sheet_name='Exp{}_LEpProfs'.format(X.split('_')[0]), index=False, header=header )
    
fid.save()
p.show('all')
