# Most recently based on the Jupyter notebook located in TTGM-2_FS19SS6/Jupyter

import numpy as n
from numpy import array, nanmean, nanstd, sqrt
from numpy import hstack, vstack, dstack
n.set_printoptions(linewidth=300,formatter={'float_kind':lambda x:'{:g}'.format(x)})
import matplotlib.pyplot as p
from mysqrtm import mysqrtm
import os
from tqdm import tqdm, trange
from sys import argv
from scipy.io import savemat, loadmat

expt = argv[1]
proj = 'TTGM-{}_FS19SS6'.format(expt)
print('\n')
print(proj)

expt = int( proj.split('_')[0].split('-')[1])
FS = int( proj.split('_')[1].split('SS')[0].split('S')[1] )
SS = int( proj.split('_')[1].split('SS')[1] )
arampath = '../{}/AramisBinary'.format(proj)
prefix = '{}_'.format(proj)
savepath = '../{}'.format(proj)

key = n.genfromtxt('../ExptSummary.dat', delimiter=',')
alpha = key[ key[:,0] == int(expt), 3 ]

os.chdir(savepath)

# [0]Stg [1]Time [2]AxSts [3]ShSts [4]AxForce [5]Torque [6]MTS Disp [7]MTS Rot
STF = n.genfromtxt('STF.dat', delimiter=',')
last = int(STF[-1,0])
maxi, maxj = n.genfromtxt('max.dat', delimiter=',', usecols=(-2,-1), unpack=True)
Xmin,Xmax,Ymin,Ymax = n.genfromtxt('box_limits.dat', delimiter=',')
profStg = n.genfromtxt('prof_stages.dat', delimiter=',', dtype=int)
LL = profStg[2]

def increm_strains(A00,A01,A10,A11,B00,B01,B10,B11):
    de00 = (2*((A00 - B00)*(A11 + B11) - (A01 - B01)*(A10 + B10))/
            ((A00 + B00)*(A11 + B11) - (A01 + B01)*(A10 + B10))
           )
    de01 = ((-(A00 - B00)*(A01 + B01) + (A00 + B00)*(A01 - B01) + 
            (A10 - B10)*(A11 + B11) - (A10 + B10)*(A11 - B11))/
            ((A00 + B00)*(A11 + B11) - (A01 + B01)*(A10 + B10))
           )
    de11 = (2*((A00 + B00)*(A11 - B11) - (A01 + B01)*(A10 - B10))/
            ((A00 + B00)*(A11 + B11) - (A01 + B01)*(A10 + B10))
           )
    if type(de00) is n.ndarray:
        return de00.mean(), de01.mean(), de11.mean()
    else:
        return de00, de01, de11

def deeq(de00, de01, de11, sig00, sig01, sig11, sigeq):
    return (sig00*de00 + sig11*de11 - 2*sig01*de01)/sigeq
    
# Averaging F  for passing points
deeq2 = n.empty((last+1,2))
# Point to point for passing points
deeq3 = n.empty((last+1,2))
# Max point from last stage traced all the way back
deeq4 = n.empty((last+1,2))
# Max point each stage
deeq5 = n.empty((last+1,2))

d = loadmat('TTGM-{}_FS19SS6_PassingPts.mat'.format(expt))
for k in trange(1,last+1):
    
    #if k%25 == 0:  print(k, end=',', flush=True)

    sig00 = STF[k,2]/2  # Hoop sts (assumed 1/2 axial)
    sig11 = STF[k,2]   # Ax sts
    sig01 = STF[k,3]   # Sh sts
    # Mises equivalent sts
    sigvm = n.sqrt(sig11**2 + sig00**2 - sig00*sig11 + 3*sig01**2)
    # Principle stresses
    s1 = sig00/2 + sig11/2 + sqrt(sig00**2 - 2*sig00*sig11 + 4*sig01**2 + sig11**2)/2
    s2 = sig00/2 + sig11/2 - sqrt(sig00**2 - 2*sig00*sig11 + 4*sig01**2 + sig11**2)/2
    # Hosford eq. sts
    sigh8 = (((s1-s2)**8+(s2-0)**8+(0-s1)**8)/2)**(1/8)

    A = d['stage_{}'.format(k)]
    if k == 1:
        B = n.zeros_like(A)
        B[:,2], B[:,-1] = 1, 1
        B[:,:2] = A[:,:2].copy()
    else:
        B = d['stage_{}'.format(k-1)]
   
    ##  deeq2: Average F for all passing!
    if k == 1:
        B00t, B11t = 1, 1
        B01t, B10t = 0, 0
    A00t, A01t, A10t, A11t = [i.mean() for i in A[:,2:].T]
    de00, de01, de11 = increm_strains(A00t,A01t,A10t,A11t,B00t,B01t,B10t,B11t)
    # For next stage, keep A and assign it to B
    B00t, B01t, B10t, B11t = [i.copy() for i in (A00t, A01t, A10t, A11t)]           
    deeq2[k,0] = (sig00*de00 + sig11*de11 - 2*sig01*de01)/sigvm    
    deeq2[k,1] = (sig00*de00 + sig11*de11 - 2*sig01*de01)/sigh8        
    
    # deeq3:  Point to Point comparison between passing vals
    InB = n.ones(len(A[:,0]), dtype=bool)
    Btemp = n.empty_like(A)
    for w, (arX, arY) in enumerate(zip(A[:,0], A[:,1])):
        Brow = B[ (B[:,0]==arX) & (B[:,1]==arY) ].ravel()
        if Brow.shape[0] == 0:
            # Then Brow is empty and those aramX,Y points from A aren't in B
            InB[w] = False
        else:
            Btemp[w] = Brow

    de00, de01, de11 = increm_strains(*[i.compress(InB) 
                                        for i in ( 
                                        (*A[:,2:].T, *Btemp[:,2:].T))])
    deeq3[k,0] = deeq(de00, de01, de11, sig00, sig01, sig11, sigvm)
    deeq3[k,1] = deeq(de00, de01, de11, sig00, sig01, sig11, sigh8)
    
    ## deeq4:  Max Point traced all the way back
    maxptA = A[ (A[:,0] == maxi[-1]) & (A[:,1] == maxj[-1]) ].ravel()
    if k == 1:
        maxptB = n.array([0,0,1,0,0,1])
    else:
        maxptB = B[ (B[:,0] == maxi[-1]) & (B[:,1] == maxj[-1]) ].ravel()
    if (len(maxptA) == 0) or (len(maxptB) == 0):
        deeq4[k] = 0
    else:
        de00, de01, de11 = increm_strains(*(*maxptA[2:], *maxptB[2:]))
        deeq4[k,0] = deeq(de00, de01, de11, sig00, sig01, sig11, sigvm)
        deeq4[k,1] = deeq(de00, de01, de11, sig00, sig01, sig11, sigh8)
    
    # deeq5:  Max point each stage
    maxptA = A[ (A[:,0] == maxi[k]) & (A[:,1] == maxj[k]) ].ravel()
    if k == 1:
        maxptB = n.array([0,0,1,0,0,1])
    else:
        maxptB = B[ (B[:,0] == maxi[k-1]) & (B[:,1] == maxj[k-1]) ].ravel()
    if (len(maxptA) == 0) or (len(maxptB) == 0):
        deeq5[k] = 0
    else:
        de00, de01, de11 = increm_strains(*(*maxptA[2:], *maxptB[2:]))
        deeq5[k,0] = deeq(de00, de01, de11, sig00, sig01, sig11, sigvm)
        deeq5[k,1] = deeq(de00, de01, de11, sig00, sig01, sig11, sigh8)

deeq2[0] = 0
deeq3[0] = 0
deeq4[0] = 0
deeq5[0] = 0

for i in [deeq2, deeq3, deeq4, deeq5]:
    i[ n.any(n.isnan(i), axis=1) ] = 0

headerline = ('[0]Stage, [1-2]AvgF-Passing-VM-H8, [3-4]PassingP2P-VM-H8, [5-6]MaxPtTracedBack-VM-H8, [7-8]MaxPtEachStage-VM-H8')
X = n.c_[STF[:,0], deeq2.cumsum(axis=0), deeq3.cumsum(axis=0), deeq4.cumsum(axis=0), deeq5.cumsum(axis=0)]
n.savetxt('IncrementalStrain.dat', X=X, header=headerline,
        fmt = '%.0f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f'
        )
