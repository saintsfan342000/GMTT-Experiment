import numpy as n
from numpy import arccos, pi, sqrt
from numpy import nanmax, nanmin, abs, array
L = n.genfromtxt
fit = n.polyfit
acos = arccos
import os
import pandas as pd
import shutil
import matplotlib.pyplot as p
p.style.use('mysty-quad')
import figfun as f

makeplot = True

def stress_measures(sig, tau):
    triax = (sig/2)/sqrt(.75*(sig**2+4*tau**2))
    qbar = -2*acos(27*sig*tau**2/(4*(3*sig**2/4 + 3*tau**2)**(3/2)))/pi + 1
    Sp = array([3*sig/4 - sqrt(sig**2 + 16*tau**2)/4, 
                    3*sig/4 + sqrt(sig**2 + 16*tau**2)/4, 0])
    mu = (2*n.median(Sp) - Sp.max() - Sp.min())/(Sp.max()-Sp.min())
    return triax, qbar, mu

# Header for each expt    
exporthead = 'Stage, AxSts(nom), ShSts(nom), Delta/L, Phi, eemax, eemean'
exporthead = exporthead.split(',')
# Header for each expt at limit load
llhead = ['Exp', 'Alpha']
[llhead.append(Z) for Z in exporthead]
# Header for tables
t1head = 'Expt, Tube, alp, Rm, t, Ecc(%), SigLL, TauLL, eemnLL'.split(', ')
t2head = 'Expt, alp, SigF, TauF, Triax, Q-bar, Mu, efmn, efmx'.split(', ')
# Header for failure summary
failhead = 'Expt, Alpha, Triax, VM-Mean, H8-Mean, Anis-Mean, VM-Max, H8-Max, Anis-Max'
failhead = failhead.split(', ')

expts = [1,2,3,4,5,7,9,10,11]

# [0]Expt No., [1]Material, [2]Tube No., [3]Alpha, [4]Alpha-True, [5]Mean Radius, [6]Thickness, [7]Eccentricity
key = n.genfromtxt('../ExptSummary.dat', delimiter=',')
key = key[ n.in1d(key[:,0],expts) ]
key = key[ key[:,3].argsort() ]
expts = key[:,0].astype(int)

# Initialize empties
LL = n.empty((len(expts),9)) #Store the LL properties
table1 = n.zeros((len(expts), len(t1head)))
table1 = pd.DataFrame(table1, index=expts, columns=t1head)
table2 = n.zeros((len(expts), len(t2head)))
table2 = pd.DataFrame(table2, index=expts, columns=t2head)
fails = n.zeros((len(expts)-1, len(failhead)))

fid = pd.ExcelWriter('../TTGM_SetData.xlsx')

for k,X in enumerate(expts):
    print(X)
    
    if X in [8,11]:
        path = '../TTGM-{}_FS32SS8'.format(X)
    else:
        path = '../TTGM-{}_FS19SS6'.format(X)
        
    alpha, tube, Rmean, thickness, exx = key[ k, [ 4,2,5,6,7] ]
    
    profstg = L(path+'/prof_stages.dat', delimiter=',').astype(int)
    
    # [0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = 0.686910 inch
    dr = L(path+'/disp-rot.dat', usecols=(0,2,3,4,5),
                          delimiter=',')
                          
    # Open file and read off the Lg
    with open(path+'/disp-rot.dat','r') as tempfid:
        line0 = tempfid.readline()
        tempfid.close()
    Lg = line0.split(' ')[-2]

    if X == 11:
        emean = L(path+'/mean.dat', delimiter=',', usecols=(11))
        emax = L(path+'/MaxPt.dat', delimiter=',', usecols=(10))
    else:
        # [0-5]Mean VM-H8-Anis-de00-01-00, [6-11]Max VM-H8-Anis-de00-01-00, [12-13]Mean, max Classic LEp
        stndata = L(path+'/IncrementalAnalysis/NewFilterResults_3med.dat', 
                        delimiter=',')
        emean, emax =  stndata[:,[0,6]].T

    d = n.c_[dr,emax,emean]
    d_stg = d[profstg].copy()
    LL[k,:2] = X, key[k,3]
    LL[k,2:] = d_stg[2]

    if X != 11:
        fails[k-1,:3] = X, alpha, stress_measures(*d_stg[-1,1:3])[0]
        fails[k-1,3:6] = stndata[-1, :3]
        fails[k-1,6:] = stndata[-1, 6:9]
        
    ## Table 1
    # 'Expt, Tube, alp, Rm, t, Ecc(%), SigLL, TauLL'
    data = [X, tube, '{:g}'.format(alpha)]
    data.append('{:.3f}\n({:.1f})'.format(Rmean, Rmean*25.4))
    data.append('{:.4f}\n({:.2f})'.format(thickness, thickness*25.4))
    data.append('{:.1f}'.format(exx))
    data.extend(['{:.1f}\n({:.1f})'.format(j,j*6.98) for j in LL[k,3:5]])
    data.append('{:.3f}'.format(LL[k,-1]))
    table1.loc[X] = data
    ## Table 2
    # 'Expt, alpha, SigF, TauF, Triax, Q-bar, Mu, efmn, efmx'
    data = [X, '{:g}'.format(alpha)]
    data.extend(['{:.1f}\n({:.1f})'.format(j,j*6.98) for j in d_stg[-1,1:3]])
    data.extend( ['{:.3f}'.format(j) for j in stress_measures(*d_stg[-1,1:3]) ])
    data.extend( ['{:.3f}'.format(j) for j in d_stg[-1,[-1,-2]] ] )
    table2.loc[X] = data
    
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
        
        if k > 0:
            p.close(fig2)
        fig2, ax2 = p.subplots()
        ax2.plot(fails[:,2], fails[:,3], 'o')
        ax2.plot(fails[:,2], fails[:,6], '^')
    
    exporthead[-1] += '   Lg={} in'.format(Lg)
    pd.DataFrame(d).to_excel( fid, sheet_name='Exp{} _ {}'.format(int(X),key[k,3]), index=False, header=exporthead)
    
    pd.DataFrame(d_stg).to_excel( fid, sheet_name='Exp{} _ {}'.format(int(X),key[k,3]), startcol = 10, index=False, header=exporthead)
    
pd.DataFrame(LL).to_excel( fid, sheet_name='LL', index=False, header=llhead )
pd.DataFrame(fails).to_excel( fid, sheet_name='FailureStrain', index=False, 
                    header=failhead)

table1.to_excel( fid, sheet_name = 'Table 1', index=False)
table2.to_excel( fid, sheet_name = 'Table 2', index=False)
    
    
fid.save()
p.show('all')

