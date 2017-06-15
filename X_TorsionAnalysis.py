import numpy as n
from numpy import (dstack, vstack, hstack,
                   linspace, array, nanmean, nanstd)
from numpy.linalg import eigvalsh
abs = n.abs
import matplotlib.pyplot as p
from scipy.interpolate import griddata, interp1d, LinearNDInterpolator
from scipy.spatial.distance import pdist
import numexpr as ne
from pandas import read_excel, read_csv
from mysqrtm import mysqrtm
from CircleFitByPratt import CircleFitByPratt as CF
from CircleFitByPratt import PrattPlotter
import os, glob, shutil, sys

'''
For analysis of pure shear DIC Data (TTGM-8)
Quite a bit different in that there is no localization to study here.
Averages over three ranges:  +/- .05", .1", and .2" (saving into d_05, d_10, d_20)
Also interpolates and saves the interpolated circumferential surface at mid-length (for buckling analysis)
'''

proj = 'TTGM-8_FS32SS8'
    
BIN = True
makecontours = True
saveAram = True                      # Whether to save the missing removed array to a .npy binary.  Only does so if BIN = False

analyze_basic_results = 1 # Whether to overwrite basic results
analyze_whole_field = 1 # Whether to take avg epsX and epsQ over whole surface
analyze_LEpProfs = 1 # Whether to make LEpProf thru max point
analyze_urProfs = 1   # Whether to make urProf thru max point and using BF circ

print(" Have you added this experiment to the summary yet?\n"*4)
print(" And are the headerlines of 'ST.dat' commented out?\n"*4)

expt = int( proj.split('_')[0].split('-')[1])
FS = int( proj.split('_')[1].split('SS')[0].split('S')[1] )
SS = int( proj.split('_')[1].split('SS')[1] )

if BIN:
    arampath = '../{}/AramisBinary'.format(proj)
    prefix = '{}_'.format(proj)
    #last = len( glob.glob( '{}/{}*.npy'.format(arampath,prefix) ) ) - 1
    # calculation of last moved to STPF[-1,0] rather than last aramfile
else:
    arampath = 'D:/Users/user/Documents/AAA Scales/{}/AllPts'.format(proj)
    # If that path don't exist then try Martin Deep drive
    if not os.path.exists(arampath):
        arampath = 'F:/GMPT/{}/AllPts'.format(proj)
    if not os.path.exists(arampath):
        raise FileNotFoundError('\nThe arampath does not exist.\n{}\n\n'.format(arampath))
    prefix = '{}-Stage-0-'.format(proj)
    #last = len( glob.glob( '{}/{}*.txt'.format(arampath,prefix) ) ) - 1

savepath = '../{}'.format(proj)
saveprefix = '{}_'.format(proj)          #Prefix for the npy files

key = n.genfromtxt('../ExptSummary.dat', delimiter=',')
if not (expt in key[:,0]):
    raise ValueError('You need to update the experiment summary!')
else:
    alpha, Rm,thickness, tube, material = key[ key[:,0] == expt, [3,5,6,2,1]].flatten()

os.chdir(savepath)

###########################################################################
# Make STF file
###########################################################################
if False:
    ArSTF = read_csv('ArSTF.dat'.format(arampath),sep=',',header=None,comment='#',skiprows=3).values
    last = int(ArSTF[-1,0])
    ## Read csv, flexible enough for spaces and tabs
    LV = read_csv('./../LVFiles/TTGM-{:.0f}_LV.dat'.format(expt),header=None,skiprows=1,comment='#',index_col=None,skipinitialspace=True,sep='\s*',engine='python').values
    STF = n.empty( (len(ArSTF[:,0]),8) )
    timecol = 2
    STF[:,[0,1]] = ArSTF[:,[0,timecol]]
    STF[0,1] = LV[0,-1] # In case LV writes its first just after 0 sec
    TRFD = interp1d(LV[:,-1],LV[:,[0,1,2,3]],axis=0).__call__(STF[:,1])
    #STF[:,[4,5]] = interp1d(LV[:,-1],LV[:,[2,0]],axis=0).__call__(STF[:,1])
    STF[:,[4,5,6,7]] = TRFD[:,[2,0,3,1]]
    STF[:,2] = STF[:,4]/(2*n.pi*Rm*thickness)/1000
    STF[:,3] = STF[:,5]/(2*n.pi*Rm*Rm*thickness)/1000
    STF[:,-2]-=STF[0,-2]    #Initialize MTS disp and rot to zero
    STF[:,-1]-=STF[0,-1]
    headerline='[0]Stg [1]Time [2]AxSts [3]ShSts [4]AxForce [5]Torque [6]MTS Disp [7]MTS Rot'
    n.savetxt('STF.dat',X=STF,fmt='%.0f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
else:
    STF = n.loadtxt('STF.dat',delimiter=',')    
    last = int(STF[-1,0])

if BIN == True:
    A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,last))
else:
    A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,last),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
    A = A[ ~n.any(n.isnan(A),axis=1), :]
    #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)

###########################################################################
# Make the contour in which we select the area we analyze for localization
###########################################################################

# Calculate the facet size and spacing relative to thickness
if not os.path.exists('facetsize.dat'):
    Atemp = A[ (abs(A[:,2])<=0.2) & (abs(A[:,3])<=0.2) , :]
    min_disp = n.min( pdist(Atemp[:,[2,3,4]]) )
    SS_th = (min_disp/thickness)
    FS_th = FS / SS * (min_disp/thickness)
    n.savetxt('facetsize.dat', X=[SS_th, FS_th], fmt='%.6f', header='[0]Step/thickness, [1]Size/thickness')
else:
    SS_th, FS_th = n.genfromtxt('facetsize.dat', delimiter=',')

LL = n.max(n.argmax(STF[:,2:4], axis=0))


######################
##### Initialize #####
#### Empty Arrays ####
######################
if True:
    # Results:  [0]Stage, [1,2,3]Eps_x,y,xy [4]Gamma, [5,6,7]Alternate Eps_x,y,gamma, [8]epeq, [9]delta/L, [10]Gamma_alt stdev, [11]Total Pts in rng, [12]Pts after filter
    # Results +/- .05"
    d_05 = n.empty((last+1,13))
    # Results +/- .1"
    d_10 = n.empty_like(d_05)
    # Results +/- 0.19"
    d_20 = n.empty_like(d_05)
    # For xzinterpolation
    xzinterp = [0]*(last+1)

def d_analsis(Lg, d):
    '''
    Define the things I do for d_05, d_10, d_20 in a function so I don't have to rewrite code twice
    '''
    rng = (abs(A[:,3]) < Lg) & (Q <= q_mid+30) & (Q >= q_mid-30)
    # Filter based on axial-shear strain ratio, alternate definitions
    ratmean = (NEy_alt/G_alt)[rng].mean()
    ratdev = (NEy_alt/G_alt)[rng].std()
    keeps = (((NEy_alt/G_alt)[rng]<=ratmean+.5*ratdev) & ((NEy_alt/G_alt)[rng]>=ratmean-.5*ratdev) )
    # Assign
    d[k,0] = k
    d[k, 1:9] = data[rng].mean(axis=0)
    d[k, 10] = G_alt.std()
    d[k, 11:13] = len(rng), rng.sum()
    # Virtual extensometer
    rng = (
            ((Q <= q_mid+20) & (Q >= q_mid-20)) 
            &
             (
                ((A[:,3]>Lg-.015) & (A[:,3] < Lg+.015)) | 
                ((A[:,3]>-Lg-.015) & (A[:,3] < -Lg+.015))
              )
           )
    Atemp = A[rng]
    xspace = n.linspace(A[rng,2].min()+.01,A[rng,2].max()-.01,len(n.unique(A[rng,0])))
    dtop, dbot = griddata( Atemp[:,[2,3]], Atemp[:,6], (xspace[None,:],array([[Lg],[-Lg]])), method='linear')
    d[k, 10] = (dtop.mean()-dbot.mean())/(2*Lg) - 1
  
######################
# Iterate Thru Stage #
######################
for itcount, k in enumerate(range(last,-1,-1)):
    print('{}. Stage {}'.format(proj,k))
    if BIN:
        A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,k))
    else:
        A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,k),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
        A = A[ ~n.any(n.isnan(A),axis=1), :]
        if saveAram:
            if not os.path.exists('./AramisBinary'):
                os.mkdir('./AramisBinary')
            n.save('./AramisBinary/{}{:.0f}'.format(saveprefix,k),A)
    #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
    if k > 0:
        Q = n.arctan2(A[:,2], A[:,4])*180/n.pi
        q_rng = Q.max()-Q.min()
        q_mid = Q.min()+q_rng/2
        F = A[:,8:12].reshape(len(A[:,0]),2,2)   # A "stack" of 2x2 deformation gradients
        FtF = n.einsum('...ji,...jk',F,F)
        U = mysqrtm( FtF )     #Explicit calculation of sqrtm
        eigU = eigvalsh(U)  #Each row is the vector of eigenvalues for that row's matrix
        LE = n.log(eigU) #Element-wise
        LE0,LE1 = LE[:,0], LE[:,1]
        LEp = ne.evaluate('( 2/3 * ( LE0**2 + LE1**2 + (-LE0-LE1)**2 ) )**0.5')
        NEx = U[:,0,0] - 1 #Hoop
        NEy = U[:,1,1] - 1 #Axial
        NExy = U[:,0,1]
        G = abs(n.arctan(NExy/(1+NEx)) + n.arctan(NExy/(1+NEy)))
        NEx_alt = F[:,0,0]-1
        NEy_alt = F[:,1,1]-1
        G_alt = abs(n.arctan(F[:,0,1]/F[:,1,1]))
        data = n.vstack((NEx, NEy, NExy, G, NEx_alt, NEy_alt, G_alt, LEp)).T
        Ro = n.sqrt(A[:,2]**2 + A[:,4]**2)
        R =  n.sqrt(A[:,5]**2 + A[:,7]**2)

        ######################
        #### Basic Results ###
        ## Lots of griddata ##
        ######################
        d_analsis(.05, d_05)
        d_analsis(.10, d_10)
        d_analsis(.20, d_20)
    else:
        d_05[k] = 0
        d_05[k,11] = d_05[k+1,11]
        d_10[k] = 0
        d_10[k,11] = d_10[k+1,11]
        d_20[k] = 0
        d_20[k,11] = d_20[k+1,11]   
        
    ######################
    #### Ur/R profiles ###
    ######################
    #  Interp z coord at a given y using a xspace and y = 0, and save these
    rng = n.abs(A[:,3])<.01
    Atemp = A[rng]
    if k == last:
        y = array([[0]])
        xlen = 2*len(n.unique(Atemp[:,0]))
        xzinterp = n.empty((xlen,2*(last+1)))
    xspace = linspace(Atemp[:,2].min()+.04, Atemp[:,2].max()-.04, xlen )
    xzinterp[:,2*k:2*k+2] = griddata( Atemp[:,[2,3]], Atemp[:,[5,7]], (xspace[None,:], y), method='linear')[0]

#########################
###### END OF LOOP ######
#########################

header0 = ("[0]Stage, [1,2,3]Eps_x,y,xy [4]Gamma, " + 
          "[5,6,7]Alternate Eps_x,y,gamma, [8]epeq, [9]delta/L, " +
          "[10]Gamma_alt stdev, [11]Total Pts in rng, [12]Pts after filter")
header1 = 'Results +/-.05 inches\n'
n.savetxt('d_05.dat', X=d_05, 
          fmt='%.0f'+', %.6f'*10+', %.0f'*2, 
          header=header1+header0)
header1 = 'Results +/-.1 inches\n'
n.savetxt('d_10.dat', X=d_10, 
          fmt='%.0f'+', %.6f'*10+', %.0f'*2, 
          header=header1+header0)
header1 = 'Results +/-.19 inches\n'
n.savetxt('d_20.dat', X=d_20, 
          fmt='%.0f'+', %.6f'*10+', %.0f'*2, 
          header=header1+header0)

# xzinterps could all be different lengths
minlen = min([len(i[:,0]) for i in xzinterp])
for i in range(last+1):
    xzinterp[i] = xzinterp[i][:minlen]
xzinterp = array(xzinterp).reshape(-1, 2*(last+1))
header = ('Interpolated deformed xz coords along y = 0 for buckling study\n' +
          'Column 2*k is kth stage deformed x-coord\n' + 
          'Column 2*k+1 is the kth stage deformed z-coord')
n.savetxt('XZInterp.dat', X=xzinterp, fmt='%.6f', delimiter=', ')

