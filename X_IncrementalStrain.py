# Most recently based on the Jupyter notebook located in TTGM-2_FS19SS6/Jupyter

import numpy as n
from numpy import array, nanmean, nanstd, sqrt
from numpy import hstack, vstack, dstack
n.set_printoptions(linewidth=300,formatter={'float_kind':lambda x:'{:g}'.format(x)})
import matplotlib.pyplot as p
from mysqrtm import mysqrtm
import os
from sys import argv
from scipy.io import savemat

expt = argv[1]
proj = 'TT2-{}_FS19SS6'.format(expt)
print('\n')
print(proj)

if expt == '18':
    alpha = n.nan
else:
    alpha = 2.0

expt = int( proj.split('_')[0].split('-')[1])
FS = int( proj.split('_')[1].split('SS')[0].split('S')[1] )
SS = int( proj.split('_')[1].split('SS')[1] )
arampath = '../{}/AramisBinary'.format(proj)
prefix = '{}_'.format(proj)
savepath = '../{}'.format(proj)

os.chdir(savepath)

# [0]Stg [1]Time [2]AxSts [3]ShSts [4]AxForce [5]Torque [6]MTS Disp [7]MTS Rot
STF = n.genfromtxt('STF.dat', delimiter=',')
last = int(STF[-1,0])
maxi, maxj = n.genfromtxt('Max10.dat', delimiter=',')[0, 1:3]
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
    
# Average F for all points in box
deeq1 = n.empty((last+1,2))
# Averaging F  for passing points
deeq2 = n.empty((last+1,2))
# Point to point for passing points
deeq3 = n.empty((last+1,2))
# Max point (no choice but point to point)
deeq4 = n.empty((last+1,2))
# A dictionary for passing points to be saved
data = {}
for k in range(1,last+1):
    if k%25 == 0:  print(k, end=',', flush=True)

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

    # Stage k
    A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,k))
    # Stage k-1
    B = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,k-1))
    # Initialize
    LEp, NEx, NEy, NExy, gamma, aramX, aramY = ( [] for _ in range(7) )
    A00,A01,A10,A11 = ( [] for _ in range(4) )
    Alex, Aley, Alexy = ( [] for _ in range(3) )    
    # Downsize A
    Abox = A[ (A[:,2]>=Xmin) & (A[:,2]<=Xmax) & (A[:,3]>=Ymin) & (A[:,3]<=Ymax), :]
    Bbox = B[ (B[:,2]>=Xmin) & (B[:,2]<=Xmax) & (B[:,3]>=Ymin) & (B[:,3]<=Ymax), :]

    ## deeq1:  Average all points in the box
    de00, de01, de11 = increm_strains(Abox[:,8].mean(), Abox[:,9].mean(), Abox[:,10].mean(), Abox[:,11].mean(),
                                   Bbox[:,8].mean(), Bbox[:,9].mean(), Bbox[:,10].mean(), Bbox[:,11].mean())
    deeq1[k,0] = deeq(de00, de01, de11, sig00, sig01, sig11, sigvm)
    deeq1[k,1] = deeq(de00, de01, de11, sig00, sig01, sig11, sigh8)

    ## deeq4:  Max Point
    maxptA = A[ (A[:,0] == maxi) & (A[:,1] == maxj) ].ravel()
    maxptB = B[ (B[:,0] == maxi) & (B[:,1] == maxj) ].ravel()
    if (len(maxptA) == 0) or (len(maxptB) == 0):
        deeq4[k] = 0
    else:
        de00, de01, de11 = increm_strains(*(*maxptA[8:], *maxptB[8:]))
        deeq4[k,0] = deeq(de00, de01, de11, sig00, sig01, sig11, sigvm)
        deeq4[k,1] = deeq(de00, de01, de11, sig00, sig01, sig11, sigh8)
    
    # Break it down
    F=Abox[:,-4:].reshape(len(Abox[:,0]),2,2)   # A "stack" of 2x2 deformation gradients 
    FtF = n.einsum('...ji,...jk',F,F)
    U = mysqrtm( FtF )     #Explicit calculation of sqrtm
    eigU = n.linalg.eigvalsh(U)  #Each row is the vector of eigenvalues for that row's matrix
    LE = n.log(eigU) #Element-wise
    LE0,LE1 = LE[:,0], LE[:,1]
    colLEp = ( 2/3 * ( LE0**2 + LE1**2 + (-LE0-LE1)**2 ) )**0.5
    colNEx = U[:,0,0] - 1
    colNEy = U[:,1,1] - 1
    colNExy = U[:,0,1]
    colG = n.arctan(colNExy/(1+colNEx)) + n.arctan(colNExy/(1+colNEy))    
    eigU, V = n.linalg.eigh(U)
    LEprin = n.log(eigU)
    # Now rotate the principle log strains back to x, y using the eigenvectors
    # It turns out that LEx, and LEy are very close to just taking n.log(U)
    # And the resulting shears are very close to the off-diagonals of U
    LExy = n.einsum('...ij,...j,...jk',V,LEprin,V)
    LEx, LEy = LExy[:,0,0], LExy[:,1,1]     
    LExy = LExy[:,0,1]    
    for j in n.unique( Abox[:,0] ):
        rng = (Abox[:,0] == j)
        Yind = Abox[ rng, 1]
        if len(Yind) == (n.max(Yind) - n.min(Yind) +1):
            locLEp = n.argmax(colLEp[rng])
            LEp.append( colLEp[rng][locLEp] )
            NEy.append( colNEy[rng][locLEp] )
            NEx.append( colNEx[rng][locLEp] )
            gamma.append( colG[rng][locLEp] )      
            aramX.append( Abox[rng,0][locLEp] )
            aramY.append( Abox[rng,1][locLEp] )
            A00.append( Abox[rng,8][locLEp])
            A01.append( Abox[rng,9][locLEp])
            A10.append( Abox[rng,10][locLEp])
            A11.append( Abox[rng,11][locLEp])
            Alex.append( LEx[rng][locLEp] )
            Aley.append( LEy[rng][locLEp] )            
            Alexy.append( LExy[rng][locLEp] )         
    # Convert these lists to arrays so I can boolean index
    for varname in 'LEp,NEy,NEx,gamma,aramX,aramY,A00,A01,A10,A11,Alex,Aley,Alexy'.split(','):
        exec('{0} = array({0})'.format(varname))
    if not n.isnan(alpha):
        ratio = NEy / gamma
    else:
        ratio = NEx / NEy
    ratioAvg = nanmean(ratio)
    ratioSDEV = nanstd(ratio)    
    passed = (ratio >= ratioAvg - 0.5 * ratioSDEV) & (ratio <= ratioAvg + 0.5 * ratioSDEV)
    # Boolean index all the arrays I'm carrying
    for varname in 'LEp,NEy,gamma,aramX,aramY,A00,A01,A10,A11,Alex,Aley,Alexy'.split(','):
        exec('{0} = {0}[passed]'.format(varname))
    
    data['stage_{}'.format(k)] = n.c_[aramX,aramY,A00,A01,A10,A11]
    
    locmax = n.argmax( LEp )
    # So now I have LEp, NEy, ..., A11 for all the passing points.  And the loc of max
    # I need to get the same aramX and Y points from stage k-1
    B = B[ n.in1d(B[:,0],aramX) & n.in1d(B[:,1],aramY) ] # This just reduced the size of B
    B = B[ (B[:,2]>=Xmin) & (B[:,2]<=Xmax) & (B[:,3]>=Ymin) & (B[:,3]<=Ymax), :]
    # First I need to make sure all the aramX,Ys in A are also in B    
    InB = n.ones(len(aramX), dtype=bool)
    B00, B01, B10, B11, aramX_B, aramY_B = ( n.empty(len(aramX)) for _ in range(6) )
    for w, (arX, arY) in enumerate(zip(aramX,aramY)):
        Brow = B[ (B[:,0]==arX) & (B[:,1]==arY) ].ravel()
        if Brow.shape[0] == 0:
            # Then Brow is empty and those aramX,Y points from A aren't in B
            InB[w] = False
        elif Brow.ndim > 1:
            raise ValueError('Your logic is bad!')
        else:
            # Then the point is also in B, so get its F components and index
            B00[w], B01[w], B10[w], B11[w] = Brow[8:]
            aramX_B[w], aramY_B[w] = Brow[:2]
    # Now get rid of all the InB=False components
    for varname in 'aramX,aramY,aramX_B,aramY_B,A00,A01,A10,A11,B00,B01,B10,B11,Alex,Aley,Alexy'.split(','):
        exec('{0} = {0}[InB]'.format(varname))
    # Symbolic calculation of dF*F.inv()
    # B is F_i-1, A is F
    # dF = F_i - F_i-1 = A - B
    # F.inv() = inv( (F_i + F_i-1)/2 ) = inv((A+B)/2)
    de00, de01, de11 = increm_strains(A00,A01,A10,A11,B00,B01,B10,B11)
    deeq3[k,0] = deeq(de00, de01, de11, sig00, sig01, sig11, sigvm)
    deeq3[k,1] = deeq(de00, de01, de11, sig00, sig01, sig11, sigh8)

    ## Last method:  Average F!
    if k == 1:
        B00t, B11t = 1, 1
        B01t, B10t = 0, 0

    A00, A01, A10, A11 = [i.mean() for i in (A00, A01, A10, A11)]

    de00, de01, de11 = increm_strains(A00,A01,A10,A11,B00t,B01t,B10t,B11t)

    B00t, B01t, B10t, B11t = [i for i in (A00, A01, A10, A11)]           
           
    deeq2[k,0] = (sig00*de00 + sig11*de11 - 2*sig01*de01)/sigvm    
    deeq2[k,1] = (sig00*de00 + sig11*de11 - 2*sig01*de01)/sigh8        

stage_1 = data['stage_1'].copy()
stage_1[:,2:] = 1, 0, 0, 1    
data['stage_0'] = stage_1
        
deeq1[0] = 0
deeq2[0] = 0
deeq3[0] = 0
deeq4[0] = 0
savemat('{}_PassingPts'.format(proj), data)

for i in [deeq1, deeq2, deeq3, deeq4]:
    i[ n.any(n.isnan(i), axis=1) ] = 0

headerline = ('[0]Stage, [1]AvgF-All in Box-VM, [2]H8, [3]AvgF-Passing-VM, [4]H8, [5]Passing P2P-VM, [6]H8, [7]MaxPt-VM, [8]H8')
X = n.c_[STF[:,0], deeq1.cumsum(axis=0), deeq2.cumsum(axis=0), deeq3.cumsum(axis=0), deeq4.cumsum(axis=0)]
n.savetxt('IncrementalStrain.dat', X=X, header=headerline,
        fmt = '%.0f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f'
        )
