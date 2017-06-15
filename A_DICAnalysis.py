import numpy as n
dstack, vstack, hstack = n.dstack, n.vstack, n.hstack
linspace = n.linspace
array = n.array
nanmean, nanstd = n.nanmean, n.nanstd
eigvalsh = n.linalg.eigvalsh
abs = n.abs
import matplotlib.pyplot as p
from scipy.interpolate import griddata, interp1d
from scipy.spatial.distance import pdist
import numexpr as ne
from pandas import read_excel, read_csv
from mysqrtm import mysqrtm
import os, glob, shutil

'''

cd 'G:\Martin_Experiments\AAA_TensionTorsion\AA_PyScripts'
python 'G:\Martin_Experiments\AAA_TensionTorsion\AA_PyScripts\TT2_DIC_Analysis.py'

Designed to be the complete TT analysis code.
Does the same analysis as Nico's code, but more:
-Techincal strains are calcualted as done by Nico, but also uses Haltom 2013 definitions.  
-Record the aramis x,y indices of the max point in every stage; written to max.dat
-Profiles pass thru the last stage max point; to accomodate this the stages are processed in reverse order
-Computes rotation and delta/L
-profStgs are determined from a plot of torque vs MTS rotation, and
taken as equally-spaced increments of rotation rather than stage number.
**
Requires a file ArSTF, which is the Aramis stage, force, time output. (or stage, time, torque)
This is used for the interpolation of the Labview data.
**
'''
proj = 'TTGM-9_FS19SS6'
BIN = False
makecontours = True
saveAram = True                      # Whether to save the missing removed array to a .npy binary.  Only does so if BIN = False

ArSTF_timecol = 2  # The column in the ArSTF containing stage time...either 1 or 2 depending on whether F is force or torque...sometimes it changes.

print("\n\nHave you added this experiment to the summary yet?")
print(" Have you added this experiment to the summary yet?\n\n")

expt = int( proj.split('_')[0].split('-')[1])
FS = int( proj.split('_')[1].split('SS')[0].split('S')[1] )
SS = int( proj.split('_')[1].split('SS')[1] )

if BIN:
    arampath = '../{}/AramisBinary'.format(proj)
    prefix = '{}_'.format(proj)
    #last = len( glob.glob( '{}/{}*.npy'.format(arampath,prefix) ) ) - 1
    # calculation of last moved to STF[-1,0] rather than last aramfile
else:
    arampath = 'D:/Users/user/Documents/AAA Scales/TTGM/{}/AllPts'.format(proj)
    prefix = '{}-Stage-0-'.format(proj)
    #last = len( glob.glob( '{}/{}*.txt'.format(arampath,prefix) ) ) - 1       

savepath = '../{}'.format(proj)
saveprefix = '{}_'.format(proj)          #Prefix for the npy files

key = n.genfromtxt('{}/../ExptSummary.dat'.format(savepath), delimiter=',')
if not (expt in key[:,0]):
    raise ValueError('You need to update the experiment summary!')
else:
    alpha, Rm,thickness, tube, material = key[ key[:,0] == expt, [3,5,6,2,1]].flatten()

os.chdir(savepath)
Size_Scott = 1/16 #[in]

###########################################################################
# Make STF file
###########################################################################
if not os.path.exists('STF.dat'):
    ArSTF = read_csv('ArSTF.dat'.format(arampath),sep=',',header=None,comment='#',skiprows=3).values
    last = int(ArSTF[-1,0])
    ## Read csv, flexible enough for spaces and tabs
    LV = read_csv('./../LVFiles/TTGM-{:.0f}_LV.dat'.format(expt),header=None,skiprows=1,comment='#',index_col=None,skipinitialspace=True,sep='\s*',engine='python').values
    STF = n.empty( (len(ArSTF[:,0]),8) )
    STF[:,[0,1]] = ArSTF[:,[0,ArSTF_timecol]]
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

###########################################################################
# Make the contour in which we select the area we analyze
###########################################################################

if BIN == True:
    A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,last))
else:
    A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,last),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
    A = A[ ~n.any(n.isnan(A),axis=1), :]

    #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
A = A[ n.abs(A[:,3])<0.2, :]        # Reduces the amount we'll have to plot

if not os.path.exists('box_limits.dat'):
    x=A[:,2]
    y=A[:,3]

    #Calculation for each facet of the Logarithmic cumulative plastic strain
    F=A[:,-4:].reshape(len(A[:,0]),2,2) 
    #FtF = n.einsum('...ij,...jk',F.transpose((0,2,1)),F)
    FtF = n.einsum('...ji,...jk',F,F) #Same as that commented out above. Kept the comment to recall how to transpose a stack of matrices
    U = mysqrtm( FtF )
    eigU = eigvalsh(U) #dimension is len(A[:,0]) x 2.  Each row is the vector of eigenvalues for that row's matrix
    LE = n.log(eigU)
    LE0,LE1 = LE[:,0], LE[:,1]
    LEp = ne.evaluate('( 2/3 * ( LE0**2 + LE1**2 + (-LE0-LE1)**2 ) )**0.5')

    #xspace=linspace(n.min(x),n.max(x),len( n.unique(A[:,0]) )/2 )
    #yspace=linspace(n.min(y),n.max(y),len( n.unique(A[:,1]) )/2 )

    #LEP = griddata( vstack((x,y)).T , LEp, (xspace[None,:],yspace[:,None]), method='linear')
    #meanLEP = nanmean(LEP.flatten())
    #stdLEP = nanstd(LEP.flatten())
    stdLEP = nanstd(LEp)
    meanLEP = nanmean(LEp)

    laststgcontour = p.figure(1,facecolor='w',figsize=(16,12))
    #p.contourf(xspace,yspace,LEP,linspace(meanLEP-0*stdLEP,meanLEP+3*stdLEP,256),extend='both',cmap='viridis')
    p.tricontourf(x,y,LEp,linspace(meanLEP-1*stdLEP,meanLEP+5*stdLEP,256),extend='both',cmap='viridis')
    p.axis([n.min(x)*1.1, 1.1*n.max(x), n.min(y), n.max(y)])
    #p.axis('equal')
    p.grid(True,linestyle='-',linewidth=2)
    p.xlabel('Undeformed X (in)')
    p.ylabel('Undeformed Y (in)')

    crop = n.asarray( p.ginput(2) )
    Xmin = n.min( crop[:,0] )
    Xmax = n.max( crop[:,0] )
    Ymin = n.min( crop[:,1] )
    Ymax = n.max( crop[:,1] )
    
    headerline='Xmin, Xmax, Ymin, Ymax (all in inches!)'
    n.savetxt('box_limits.dat',X=array([Xmin, Xmax, Ymin, Ymax])[None,:],fmt='%.6f',delimiter=', ',header=headerline)
    
    box_x = array([Xmin , Xmax , Xmax , Xmin , Xmin])
    box_y = array([Ymin , Ymin , Ymax , Ymax , Ymin])
    p.plot(box_x,box_y,'w',linewidth=2.5)
    p.grid(False)
    p.draw()
else:
    Xmin,Xmax,Ymin,Ymax = n.loadtxt('box_limits.dat',delimiter=',')
    box_x = array([Xmin , Xmax , Xmax , Xmin , Xmin])
    box_y = array([Ymin , Ymin , Ymax , Ymax , Ymin])

A = A[ (abs(A[:,2])<=0.2) & (abs(A[:,3])<=0.1) , :]
min_disp = n.min( pdist(A[:,[2,3,4]]) )
SS_th = (min_disp/thickness)
FS_th = FS / SS * (min_disp/thickness)

#Definition of the size of the averaging zone for grid method
size_av = Size_Scott - thickness * FS_th

#Stage Time Force data...Create the 10 stations at which we'll plot profiles
LL = n.max(n.argmax(STF[:,2:4], axis=0))
# For TTGM-4
if expt == 4:
    LL = n.nonzero(STF[:,0]==519)[0][0]

if not os.path.exists('prof_stages.dat'):
    #[0]Stg [1]Time [2]AxSts [3]ShSts [4]AxForce [5]Torque [6]MTS Disp [7]MTS Rot
    TRFD = read_csv('STF.dat',sep=',',comment='#',header=None,index_col=None).values[:,[5,7,4,6]]
    p.figure()
    p.plot(TRFD[:,1],TRFD[:,0])
    p.plot(TRFD[LL,1],TRFD[LL,0],'ro',ms=10)
    p.title('Click where the first\ncalc stage will be', fontsize=30)
    #p.axis([0,1.1*TRFD[-1,1],0,1.2*n.max(TRFD[:,0])])
    p.gca().set_ylim([0,1.2*n.max(TRFD[:,0])])
    p.gca().set_xlim(left=0)
    yldrot = n.asarray(p.ginput(1)).flatten()[0]
    p.close()

    # Make evenly spaced increments
    for i in range(len(TRFD[:,1])):
        if TRFD[i,1] >= yldrot:
            break

    #Two before the limit load, seven after
    x = hstack( (n.linspace(yldrot,TRFD[LL,1],3)[0:2], n.linspace(TRFD[LL,1],TRFD[-1,1],8)) )
    profStg = n.empty(len(x))
    for i in range(len(x)):
        profStg[i] = n.where(TRFD[:,1]>=x[i])[0][0]
        
    headerline = 'Stages at which profiles were generated'
    n.savetxt('prof_stages.dat',X=profStg[None,:],fmt='%.0f',delimiter=', ',header=headerline)
else:
    profStg = n.loadtxt('prof_stages.dat',delimiter=',')

if not os.path.exists('disp-rot_limits.dat'):
    ## Now identify the upper and lower ranges for calculating delta and ph
    
    if BIN == True:
        A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,0))
    else:
        A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,0),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
        A = A[ ~n.any(n.isnan(A),axis=1), :]
    #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
    p.figure(figsize=(16,12))
    p.title('Click the four points that bound the cusp of the thick edges',size=20)
    p.plot(A[:,3],A[:,4],'.',alpha=0.5,ms=4,markevery=2)
    p.gcf().gca().set_ylim([.98*1.826/2,1.02*1.9685/2])
    p.gcf().gca().set_xlim(right=-.1)
    rdlim = n.sort(n.asarray(p.ginput(2))[:,0])
    p.close()
    p.figure(figsize=(16,12))
    p.plot(A[:,3],A[:,4],'.',alpha=0.5,ms=4,markevery=2)
    p.gcf().gca().set_ylim([.98*1.826/2,1.02*1.9685/2])
    p.gcf().gca().set_xlim(left=.1)
    rdlim = n.sort( hstack( (rdlim, n.asarray(p.ginput(2))[:,0]) ) )
    p.close()
    headerline='Lower sxn Ymin, Lower sxn Ymax, Upper sxn Ymin, Upper sxn Ymax'
    n.savetxt('disp-rot_limits.dat',X=rdlim[None,:],fmt='%.6f',delimiter=', ',header=headerline)
else:
    rdlim = n.sort(n.loadtxt('disp-rot_limits.dat',delimiter=','))

#Initialize a few things before looping and calculating every stage
export_max=n.zeros( (last+1,13) )   #MaxPt data
export_mean=n.zeros( (last+1,12) )  #MeanPt data
export_stdv=n.zeros( (last+1,12) )  #Std Deviation data
export_Scott=n.zeros( (last+1,12) ) #DIC Macro data
export_MaxPt=n.zeros( (last+1,11) ) #The last stage's max point traced thru all stages
prof_count=0                        #For indexing the profLEp list        
profLEp = [0]*len(profStg)          #List of proLEp arrays
profUr = n.empty( (201,len(profStg)+1) )    # First col. for y-coord, first row for stg no.
L = n.zeros(last+1)                   #For delta/L
up_th = n.zeros(last+1)               #For phi
lo_th = n.zeros(last+1)               #For phi

#Cycle through the stages
z=1
for k in range(last,-1,-1):
    print('{}. Stage {}'.format(proj,k))
    if k == 0:
        if BIN:
            A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,k))
        else:
            A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,k),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
            A = A[ ~n.any(n.isnan(A),axis=1), :]
            if saveAram:
                if not os.path.exists('./AramisBinary'):
                    os.mkdir('./AramisBinary')
                n.save('./AramisBinary/{}{:.0f}'.format(saveprefix,k),A)
        export_max[k] =  [STF[k,0], STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0,0,0]
        export_MaxPt[k] = [STF[k,0], STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0]
        export_mean[k] = [STF[k,0], 0, STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0]
        export_stdv[k] = [STF[k,0], 0, STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0]
        export_Scott[k] = [STF[k,0], 0, STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0]
    else:
        LEp, NEx, NEy, NExy, gamma, xcoord, ycoord, aramX, aramY, NEx_alt, NEy_alt, gamma_alt = ( [] for _ in range(12) )
        scot_count=0                #Count the number of points in the Scott method
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
        Abox = A[ (A[:,2]>=Xmin) & (A[:,2]<=Xmax) & (A[:,3]>=Ymin) & (A[:,3]<=Ymax), :]
        F=Abox[:,-4:].reshape(len(Abox[:,0]),2,2)   # A "stack" of 2x2 deformation gradients 
        FtF = n.einsum('...ji,...jk',F,F)
        U = mysqrtm( FtF )     #Explicit calculation of sqrtm
        eigU = eigvalsh(U)  #Each row is the vector of eigenvalues for that row's matrix
        LE = n.log(eigU) #Element-wise
        LE0,LE1 = LE[:,0], LE[:,1]
        colLEp = ne.evaluate('( 2/3 * ( LE0**2 + LE1**2 + (-LE0-LE1)**2 ) )**0.5')
        colNEx = U[:,0,0] - 1
        colNEy = U[:,1,1] - 1
        colNExy = U[:,0,1]
        colG = n.arctan(colNExy/(1+colNEx)) + n.arctan(colNExy/(1+colNEy))
        colNEx_alt = F[:,0,0]-1
        colNEy_alt = F[:,1,1]-1
        colG_alt=n.arctan(F[:,0,1]/F[:,1,1]);
                
        # What we have effectively done is created an additional column for each data point.
        # Now I just need to sort though those which are in passable columns
        
        for j in n.unique( Abox[:,0] ):
            rng = (Abox[:,0] == j)
            Yind = Abox[ rng, 1]
            if len(Yind) == (n.max(Yind) - n.min(Yind) +1):
                locLEp = n.argmax(colLEp[rng])                #Location of...
                LEp.append( colLEp[rng][locLEp] )             #Max LEp in the current column
                NEx.append( colNEx[rng][locLEp] )          
                NEy.append( colNEy[rng][locLEp] )
                gamma.append( colG[rng][locLEp] )      
                xcoord.append( Abox[rng,2][locLEp] )
                ycoord.append( Abox[rng,3][locLEp] )
                aramX.append( Abox[rng,0][locLEp] )
                aramY.append( Abox[rng,1][locLEp] )
                NEx_alt.append( colNEx_alt[rng][locLEp] )
                NEy_alt.append( colNEy_alt[rng][locLEp] )
                gamma_alt.append( colG_alt[rng][locLEp] )

        LEp, NEx, NEy, gamma, xcoord, ycoord, aramX, aramY = map(array,[LEp, NEx, NEy, gamma, xcoord, ycoord, aramX, aramY])    #Convert lists to arrays
        NEx_alt, NEy_alt, gamma_alt = map(array,[NEx_alt, NEy_alt, gamma_alt])
        ratio = NEy / gamma
        ratioAvg = nanmean(ratio)
        ratioSDEV = nanstd(ratio)
        passed = (ratio >= ratioAvg - 0.5 * ratioSDEV) & (ratio <= ratioAvg + 0.5 * ratioSDEV)
        LEp=LEp[passed]
        NEx=NEx[passed]
        NEy=NEy[passed]
        gamma=gamma[passed]
        NEx_alt=NEx_alt[passed]
        NEy_alt=NEy_alt[passed]
        gamma_alt=gamma_alt[passed]        
        xcoord=xcoord[passed]
        ycoord=ycoord[passed]
        aramX = aramX[passed]
        aramY = aramY[passed]
                
        ## Max 10 stages and identify I,J of max point in last
        locmax = n.argmax( LEp )
        if k == last:
            aramXmaxlast, aramYmaxlast = aramX[locmax], aramY[locmax]   #Save for making strain profiles
            MaxTen = n.flipud(array( vstack( (LEp,aramX,aramY,xcoord,ycoord) ) ).T[n.argsort(LEp),:][-10:,:])
        
        export_max[k] =  [ STF[k,0], STF[k,1], STF[k,2], STF[k,3], NEx[locmax], NEy[locmax], abs(gamma[locmax]), 
                                NEx_alt[locmax], NEy_alt[locmax], abs(gamma_alt[locmax]), LEp[locmax], aramX[locmax], aramY[locmax] ]
        export_mean[k] = [ STF[k,0],  STF[k,1], sum(passed), STF[k,2], STF[k,3], nanmean(NEx), nanmean(NEy), abs(nanmean(gamma)), 
                                nanmean(NEx_alt), nanmean(NEy_alt), abs(nanmean(gamma_alt)), nanmean(LEp)]
        export_stdv[k] = [ STF[k,0], STF[k,1], sum(passed), STF[k,2], STF[k,3], nanstd(NEx), nanstd(NEy), abs( nanstd(gamma) ), 
                                nanstd(NEx_alt), nanstd(NEy_alt), abs(nanstd(gamma_alt)), nanstd(LEp)]
        ## Store Last stage max point data
        ## Find the row in Abox (or depth-layer in F, or index in LEp, etx.) where the max point is
        rowmax = n.nonzero( (Abox[:,0]==aramXmaxlast) & (Abox[:,1]==aramYmaxlast) )[0]
        if len(rowmax) != 0:
            export_MaxPt[k] = [ STF[k,0], STF[k,1], STF[k,2], STF[k,3], colNEx[rowmax], colNEy[rowmax], colG[rowmax], 
                                colNEx_alt[rowmax], colNEy_alt[rowmax], colG_alt[rowmax], colLEp[rowmax] ]
        else:
            export_MaxPt[k] = [ STF[k,0], STF[k,1], STF[k,2], STF[k,3], n.nan, n.nan, n.nan, 
                                n.nan, n.nan, n.nan, n.nan ]

        ######################################################
        #Average using Scot-size box around the max point only
        ######################################################
        rgn = (Abox[:,2] <= (xcoord[locmax]+size_av/2)) & (Abox[:,2] >= (xcoord[locmax]-size_av/2)) & (Abox[:,3] <= (ycoord[locmax]+size_av/2)) &  (Abox[:,3] >= (ycoord[locmax]-size_av/2))
        
        OurSize_Scott = FS_th * thickness + .5*( (n.max(Abox[rgn,2]) - n.min(Abox[rgn,2])) +( n.max(Abox[rgn,3]) - n.min(Abox[rgn,3])) )
        export_Scott[k] = [STF[k,0], STF[k,1], nanmean(OurSize_Scott), STF[k,2], STF[k,3], nanmean(colNEx[rgn]), nanmean(colNEy[rgn]),
                            n.abs(nanmean(colG[rgn])), nanmean(colNEx_alt[rgn]), nanmean(colNEy_alt[rgn]), n.abs(nanmean(colG_alt[rgn]))
                            , nanmean(colLEp[rgn])]                                
                                
        ######################################################
        # Profile of ep, radial contraction, contours for different stations
        ######################################################
        if k in profStg:
            ### Profiles
            Acol = A[ A[:,0] == aramXmaxlast, :]
            ystn = n.empty( (len(Acol[:,0]),3) )    # 2-columns [0]def-ycoord [1]LEp
            ystn[:,[0,1]] = Acol[:,[3,6]] #Defored y-coords
            F=Acol[:,-4:].reshape(len(Acol[:,0]),2,2)
            FtF = n.einsum('...ji,...jk',F,F) 
            U = mysqrtm( FtF )
            eigU = eigvalsh(U)  #Each row is the vector of eigenvalues for that row's matrix
            LE = n.log(eigU)
            LE0,LE1 = LE[:,0], LE[:,1]
            ystn[:,2] = ne.evaluate('( 2/3 * ( LE0**2 + LE1**2 + (-LE0-LE1)**2 ) )**0.5')
            ystn = n.flipud(ystn[ n.argsort(ystn[:,0]) ])
            profLEp[-(1+prof_count)] = ystn #Want the profiles put in ascendng order
               
            ### Radial Contraction, interpolating based on undeformed coord, will be plotted against undeformed
            #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
            if prof_count == 0:
                yspace = linspace( -.4, .4, 200 )[:,None]    #Test sxt including rad is ~0.62" long
                profUr[1:,[0]] = yspace #First column is the yspace
            profUr[0,-(prof_count+1)] = k
            # Work on just a part of the point cloud to speed up griddata
            Acol = A[ (A[:,2]>=-.25) & (A[:,2]<=.25) , :]
            xo, zo, x, z = Acol[:,[2,4,5,7]].T
            er = ne.evaluate('(sqrt( x**2 + z**2 ) - sqrt(xo**2 + zo**2 ))/sqrt(xo**2 + zo**2 )')
            erint = griddata( Acol[:,[2,3]], er, (n.array([[0]]),yspace), method='linear')
            profUr[1:,[-(1+prof_count)]] = erint
            
            prof_count=prof_count+1    
            ######################################################
            # Strain contours
            ######################################################        
            if makecontours:
                if k == last:
                    try:
                        p.sca(laststgcontour.gca())
                    except NameError:
                        Acol = A[ abs(A[:,3])<=0.2, :]
                        F=Acol[:,-4:].reshape(-1,2,2) 
                        FtF = n.einsum('...ji,...jk',F,F) #Same as that commented out above. Kept the comment to recall how to transpose a stack of matrices
                        U = mysqrtm( FtF )
                        eigU = eigvalsh(U) #dimension is len(A[:,0]) x 2.  Each row is the vector of eigenvalues for that row's matrix
                        LE = n.log(eigU)
                        LE0,LE1 = LE[:,0], LE[:,1]
                        LEp = ne.evaluate('( 2/3 * ( LE0**2 + LE1**2 + (-LE0-LE1)**2 ) )**0.5')
                        stdLEP = nanstd(LEp)
                        meanLEP = nanmean(LEp)
                        laststgcontour = p.figure(1,facecolor='w',figsize=(16,12))
                        p.tricontourf(Acol[:,2],Acol[:,3],LEp,linspace(meanLEP-1*stdLEP,meanLEP+5*stdLEP,256),extend='both',cmap='viridis')
                        p.axis([n.min(Acol[:,2]), n.max(Acol[:,2]), n.min(Acol[:,3]), n.max(Acol[:,3])])
                        p.grid(True,linestyle='-',linewidth=2)
                        p.xlabel('Undeformed X (in)')
                        p.ylabel('Undeformed Y (in)')
                    z=1
                    for m in range(len(MaxTen[:,0])):
                        p.plot(MaxTen[m,3],MaxTen[m,4],marker='${}$'.format(z),ms=8,mfc='k',mec='k')
                        z+=1
                    cbar = p.colorbar()
                    cbar.set_label(r'     $\mathsf{e_{e}^{p}}$',rotation=0,fontsize=25)
                    p.title('$\mathsf{{e^{{p}}_{{e}}}}$: Station {} (Last Stage)'.format(n.where(profStg==k)[0][0]+1))
                    p.savefig('last_stage.png')
                else:
                    Acol = A[ abs(A[:,3])<=0.2, :]
                    F = Acol[:,8:].reshape(len(Acol[:,0]),2,2)
                    FtF = n.einsum('...ji,...jk',F,F) 
                    U = mysqrtm( FtF )
                    eigU = eigvalsh(U)  #Each row is the vector of eigenvalues for that row's matrix
                    LE = n.log(eigU)
                    LE0,LE1 = LE[:,0], LE[:,1]
                    LEp = ne.evaluate('( 2/3 * ( LE0**2 + LE1**2 + (-LE0-LE1)**2 ) )**0.5')
                    #xspace=linspace(n.min(x),n.max(x),len( n.unique(Acol[:,0]) )/2 )
                    #yspace=linspace(n.min(y),n.max(y),len( n.unique(Acol[:,1]) )/2 )
                    #LEP = griddata( vstack((x,y)).T , LEp, (xspace[None,:],yspace[:,None]), method='linear')
                    meanLEP = nanmean(LEp)
                    stdLEP = nanstd(LEp)
                    p.figure(facecolor='w',figsize=(16,12))
                    p.tricontourf(Acol[:,2],Acol[:,3],LEp,linspace(meanLEP-1*stdLEP,meanLEP+3*stdLEP,256),extend='both',cmap='viridis')
                    #p.axis('equal')
                    p.xlabel('Undeformed X (in)')
                    p.ylabel('Undeformed Y (in)')
                    p.axis([n.min(Acol[:,2]), n.max(Acol[:,2]), n.min(Acol[:,3]), n.max(Acol[:,3])])
                    p.plot(xcoord,ycoord,'o',mec='k',mfc='k',ms=6)
                    p.plot(box_x,box_y,'w',linewidth=2.5)
                    p.xlabel('Undeformed X (in)')
                    p.ylabel('Undeformed Y (in) (axial)')
                    if k == LL:
                        p.title('$\mathsf{{e^{{p}}_{{e}}}}$: Station {} (Limit Load)'.format(n.where(profStg==k)[0][0]+1))
                    else:
                        p.title('$\mathsf{{e^{{p}}_{{e}}}}$: Station {}'.format(n.where(profStg==k)[0][0]+1))
                    cbar = p.colorbar()
                    cbar.set_label(r'     $\mathsf{e_{e}^{p}}$',rotation=0,fontsize=25)
                if not os.path.exists('./Contours'):
                    os.mkdir('./Contours')
                p.savefig('./Contours/Stg{:.0f}'.format(k))
                p.close()
            else:
                p.close()
            
        ###############################################
        # Calculate delta and phi
        ###############################################

    lower = A[ (A[:,3] >= rdlim[0]) & (A[:,3] <= rdlim[1]) & (n.abs(A[:,2]) <= 0.3), :]
    upper = A[ (A[:,3] >= rdlim[2]) & (A[:,3] <= rdlim[3]) & (n.abs(A[:,2]) <= 0.3), :]
    L[k] = n.mean(upper[:,6]) - n.mean(lower[:,6])
    up_th[k] = n.mean( n.arctan2(upper[:,7], upper[:,5]) - n.arctan2(upper[:,4], upper[:,2]) )*180/n.pi
    lo_th[k] = n.mean( n.arctan2(lower[:,7], lower[:,5]) - n.arctan2(lower[:,4], lower[:,2]) )*180/n.pi


#########
# End of Loop
#########

# delta/L and phi
dL = (L - L[0]) / L[0]
dth = -(lo_th - up_th)

# Save data files!
#########
#Max.dat
headerline='[0]Stage [1]Time [2]AxSts [3]ShSts [4]NEx [5]NEy [6]Gamma [7]F11-1 [8]F22-1 [9]atan(F12/F22) [10]epeq [11]AramX [12]AramY'
n.savetxt('max.dat', X=export_max, fmt='%.0f, %.2f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.0f, %.0f',header=headerline)
#mean.dat
headerline='[0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
n.savetxt('mean.dat', X=export_mean, fmt='%.0f, %.2f, %.0f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
#std.dat
headerline='[0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
n.savetxt('std.dat', X=export_stdv, fmt='%.0f, %.2f, %.0f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
#like_Scott
headerline='[0]Stage [1]Time [2]SizeAveragingZone(in) [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
n.savetxt('like_Scott.dat', X=export_Scott, fmt='%.0f, %.2f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
#MaxPt.dat
headerline='Last Stage Last Point, traced thru whol exp.\n[0]Stage [1]Time [2]AxSts [3]ShSts [4]NEx [5]NEy [6]Gamma [7]F11-1 [8]F22-1 [9]atan(F12/F22) [10]epeq'
n.savetxt('MaxPt.dat', X=export_MaxPt, fmt='%.0f, %.2f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
#disp-rot
headerline='[0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = {:.6f} inch'.format(L[0])
n.savetxt('disp-rot.dat', X=hstack( (STF[:,[0,1,2,3]],dL[:,None],dth[:,None]) ), fmt='%.0f, %.2f, %.6f, %.6f, %.6f, %.6f',header=headerline)
# Save the Max10
headerline = '[0]Epeq [1]AramXIndex [2]AramYIndex [3]UndefXCoord [4]UndefYCoord'
n.savetxt('Max10.dat',X=MaxTen,fmt='%.6f, %.0f, %.0f, %.6f, %.6f',header=headerline)
# Save the radial contraction
headerline = 'First Row Stage number.  Second row begin data.\n[0]Ycoord [1:] Stage Ur/Ro'
n.savetxt('RadialContraction.dat',X=profUr,fmt = '%.6f',delimiter=',',header=headerline)

# Save profiles (takes some work)
maxprolen = 0
for i in range(len(profStg)):
    prolen = len(profLEp[i][:,0])   
    # If you get the "int object not subscriptable" warning, it likely means not every stage in profStg was iterated over
    # This typically occurs when a stage in profStg exceeds "last", which happens from time to time when I have to deviate from the norm
    if prolen > maxprolen:
        maxprolen = prolen
allprofs = n.zeros( (1+maxprolen,len(profStg)*3) )
allprofs[allprofs == 0] = n.nan
for i in range(len(profStg)):
    allprofs[0,3*i] = profStg[i]
    allprofs[1:len(profLEp[i][:,0])+1,[3*i,3*i+1,3*i+2]] = profLEp[i]
headerline = 'First Row: Stage number. 2nd Row Begin Column Triplets of [0]Undef y, [1]Def y, [2] vs epeq'
n.savetxt('StrainProfiles.dat',X=allprofs,fmt='%.6f',delimiter=',',header=headerline) 

os.system('python ../AA_PyScripts/B_Plots.py {} {} {} {}'.format(int(expt), int(FS), int(SS), savepath))



'''
### Some automation stuff used for when I recomputed on a different FSSS
### Paste this at the top and comment out the existing expt through saveprefix
proj = 'TT2-31_FS32SS8'
expt = int( proj.split('_')[0].split('-')[1])
FS = int( proj.split('_')[1].split('SS')[0].split('S')[1] )
SS = int( proj.split('_')[1].split('SS')[1] )

arampath = 'F:/TensTors/{}/AllPts'.format(proj)
prefix = '{}-Stage-0-'.format(proj)
last = len( glob.glob( '{}/{}*.txt'.format(arampath,prefix) ) ) - 1

BIN = False
if BIN:
    arampath = '../{}/AramisBinary'.format(proj)
    prefix = '{}_'.format(proj)
    last = len( glob.glob( '{}/{}*.npy'.format(arampath,prefix) ) ) - 1

savepath = '../{}'.format(proj)
makecontours = True
saveAram = True                         # Whether to save the missing removed array to a .npy binary
saveprefix = '{}_'.format(proj)          #Prefix for the npy files

existing = glob.glob('../TT2-{}_*'.format(expt))[0]

if os.path.exists(savepath):
    # Throw and error because I don't wanna overwrite anything!
    AAAAA += 1
    pass
else:
    os.mkdir(savepath)
    
shutil.copy( '{}/ArSTF.dat'.format(existing), '{}/ArSTF.dat'.format(savepath) )
shutil.copy( '{}/prof_stages.dat'.format(existing), '{}/prof_stages.dat'.format(savepath) )
'''
