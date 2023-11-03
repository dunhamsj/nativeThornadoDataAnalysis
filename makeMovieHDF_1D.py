#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( './publication.sty' )

from UtilitiesModuleHDF import ReadFieldsHDF, GetNorm

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---
THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True)
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

############################ User Input ############################

# --- Define root path for data ---

RootPath = THORNADO_DIR + 'SandBox/YahilCollapse_XCFC/'
suffix = 'Output/'

#Problem = 'Advection'
#Problem = 'RiemannProblem'
#Problem = 'RiemannProblem2d'
#Problem = 'RiemannProblemSpherical'
#Problem = 'SphericalSedov'
#Problem = 'KelvinHelmholtz_Relativistic'
#Problem = 'KelvinHelmholtz'
##Problem = '2D_Jet'
#Problem = 'StandingAccretionShock'
#Problem = 'GravitationalCollapse'
#Problem = 'DynamicTOV'
Problem = 'YahilCollapse'

CoordinateSystem = 'SPHERICAL'
UsePhysicalUnits = True
RunTime = 10.0 # seconds

Field             = 'PF_V1'
Dimension         = ['X1']#,'X2']
UseSemiLogXScale  = True
UseSemiLogYScale  = False
PlotTwoVariables  = False
cmap              = 'Reds'
PlotTranspose     = False
UseGeometryFields = True
SnapshotRange     = [0,542]

FigTitle = suffix[:-1]
FigTitle = Problem

############################

nFiles = SnapshotRange[1] - SnapshotRange[0] + 1

Snapshots \
  = np.linspace( SnapshotRange[0], SnapshotRange[1], nFiles, \
                 dtype = np.int64 )

# --- Define where to look for data ---
PathToData = RootPath + suffix + Problem

nDimsX = len( Dimension )

SaveFileAs = 'mov.{:}_{:}_notdiff.mp4'.format( Problem, Field )

if( PlotTwoVariables and nDimsX > 1 ):
    exit('Cannot plot two variables in two dimensions. Exiting...')

Names = ReadFieldsHDF \
          ( PathToData, Snapshots, \
            CoordinateSystem = CoordinateSystem, \
            UsePhysicalUnits = UsePhysicalUnits, \
            UseGeometryFields = UseGeometryFields )
TimeUnit = Names['Time'][0]
Time     = Names['Time'][1]

if( np.any( np.array( Dimension ) == 'X1' ) ):
    LengthUnit = Names[Dimension[0]][0]
    x1    = np.array( Names[Dimension[0]][1] )
    x1lim = ( np.min(x1), np.max(x1) )
if( np.any( np.array( Dimension ) == 'X2' ) ):
    x2    = np.array( Names[Dimension[1]][1] )
    x2lim = ( np.min(x2), np.max(x2) )
if( np.any( np.array( Dimension ) == 'X3' ) ):
    x3    = np.array( Names[Dimension[2]][1] )
    x3lim = ( np.min(x3), np.max(x3) )

Y1 = Names[Field][1]
y1Min = +np.inf
y1Max = -np.inf
for t in range( nFiles ):
    y1Min = min( y1Min, Y1[t][0,:,:].min() )
    y1Max = max( y1Max, Y1[t][0,:,:].max() )
ylim = [ y1Min, y1Max ]

if PlotTwoVariables:

    PathToData = RootPath + 'Output_Yahil_GravitationalMass/' + Problem
    Names = ReadFieldsHDF \
              ( PathToData, Snapshots, \
                CoordinateSystem = CoordinateSystem, \
                UsePhysicalUnits = UsePhysicalUnits, \
                UseGeometryFields = UseGeometryFields )

    Y2 = Names[Field][1]
    y2Min = +np.inf
    y2Max = -np.inf
    for t in range( nFiles ):
        y2Min = min( y2Min, Y2[t][0,0,:].min() )
        y2Max = max( y2Max, Y2[t][0,0,:].max() )

    ylim = [ min( ylim[0], y2Min ), max( ylim[1], y2Max ) ]

################ Plotting information

if( nDimsX == 1 ):

    # Animation program adapted from
    # https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
    fig, ax = plt.subplots()
    ax.set_title( r'$\texttt{{{:}}}$'.format( FigTitle ) )

    ax.set_xlim( x1lim )
    ax.set_ylim( [ ylim[0]-0.1, ylim[1]+0.1 ]  )

    ax.set_xlabel( '{:} {:}'.format( Dimension[0], LengthUnit ) )
    ax.set_ylabel( '{:}'.format( Field ) )

    if UseSemiLogYScale:
        if( y1Min < 0.0 ):
            ax.set_yscale( 'symlog' )
        else:
            ax.set_yscale( 'log' )
        Height = np.log10( ylim[1] ) - np.log10( ylim[0] )
        ytext = 10**( np.log10( ylim[0] ) + 0.9 * Height )
    else:
        Height = ylim[1] - ylim[0]
        ytext = ylim[0] + 0.9 * Height

    if UseSemiLogXScale:
        ax.set_xscale( 'log' )
        Width = np.log10( x1lim[1] ) - np.log10( x1lim[0] )
        xtext = 10**( np.log10( x1lim[0] ) + 0.1 * Width )
    else:
        Width = x1lim[1] - x1lim[0]
        xtext = x1lim[0] + 0.1 * Width

    if not PlotTwoVariables:
        line1, = ax.plot([],[],'k-')
    else:
        line1, = ax.plot([],[], label = 'BaryonicMass')
        line2, = ax.plot([],[],label = 'GravitationalMass')

    time_text = plt.text( xtext, ytext, '' )
    #plt.legend()

    ax.grid()

    X1_C = np.array( Names['X1_C'][1] )
    X1   = np.array( Names['X1'  ][1] )

    nNodes = X1.shape[0] // X1_C.shape[0]
    nX1    = X1.shape[0] // nNodes

    def computeCellAverage( Names, nNodes, Field, nSS ):

      x1_C = np.empty( nX1, np.float64 )

      uh = Names[Field][1]

      if   nNodes == 1:
        wq = np.array( [ 1.0 ], np.float64 )
      elif   nNodes == 2:
        wq = np.array( [ 0.5, 0.5 ], np.float64 )
      elif nNodes == 3:
        wq = np.array( [ 5.0, 8.0, 5.0 ], np.float64 ) / 18.0
      else:
        exit( 'Not available for nNodes = {:}'.format( nNodes ) )

      SqrtGm = Names['GF_Sg'][1]

      uK = np.empty( (nSS,nX1), np.float64 )

      for iSS in range( nSS ):
        for iX1 in range( nX1 ):

          iLo = nNodes * iX1
          iHi = iLo + nNodes

          vK = np.sum( wq * SqrtGm[iSS,0,0,iLo:iHi] )

          uK[iSS,iX1] \
            = np.sum( wq * uh[iSS,0,0,iLo:iHi] * SqrtGm[iSS,0,0,iLo:iHi] ) / vK

          if iSS == 0:
            x1_C[iX1] = np.sum( wq * x1[iLo:iHi] )

      return uK, x1_C

    uK, X1_C \
      = computeCellAverage( Names, nNodes, Field, nFiles )

    header  = 'data[0,1:]  = X1_C [km]\n'
    header += 'data[1:,0]  = Time [ms]\n'
    header += 'data[1:,1:] = uK'

    uKK = np.empty( (nFiles+1,nX1+1), np.float64 )
    uKK[0,0 ] = np.nan
    uKK[0,1:] = X1_C
    uKK[1:,0] = Time
    uKK[1:,1:] = uK

    np.savetxt( '{:}_native_{:}.dat'.format( Problem, Field ), \
                uKK, header = header )

    os.system( 'rm -rf __pycache__' )
    exit()

    # Intialize each new frame
    if not PlotTwoVariables:
        def InitializeFrame():
            line1.set_data([],[])
            time_text.set_text('')
            return line1, time_text
    else:
        def InitializeFrame():
            line1.set_data([],[])
            line2.set_data([],[])
            time_text.set_text('')
            return line1, line2,

    # Animation function
    if not PlotTwoVariables:
        def UpdateFrame(t):
            print( '  {:}/{:}'.format( t, nFiles ) )
            y1 = Y1[t][0,0,:]
            line1.set_data( x1, y1 )
            time_text.set_text('time = {:.3e} {:}'.format( Time[t], TimeUnit ) )
            return line1, time_text
    else:
        def UpdateFrame(t):
            y1 = Y1[t][0,0,:]
            y2 = Y2[t][0,0,:]
            line1.set_data( x1, y1 )
            line2.set_data( x1, y2 )
            time_text.set_text('time = {:.3e} {:}'.format( Time[t], TimeUnit ) )
            return line1, line2, time_text

    # Call the animator
    anim = animation.FuncAnimation \
             ( fig, \
               UpdateFrame, \
               init_func = InitializeFrame, \
               frames    = nFiles, \
               interval  = 1000, \
               blit      = True )

    anim.save( SaveFileAs, fps = max( 1, int( nFiles / RunTime ) ), dpi = 300 )

else:

    # 2D animation program adapted from
    # https://matplotlib.org/examples/animation/dynamic_image.html
    fig = plt.figure()

    fig.suptitle( FigTitle, fontsize = 20 )
    def f(t):
        return Y1[t][0,:,:]

    if( PlotTranspose ):
        fig.suptitle( '|{:}-{:}.T|'.format( Field, Field ), \
          fontsize = 20 )
        def f(t):
            return np.abs( Y1[t][0,:,:] - Y1[t][0,:,:].T )

    vmin = +np.inf
    vmax = -np.inf
    for t in range( nFiles ):
        vmin = min( vmin, np.min( Y1[t][0,:,:] ) )
        vmax = max( vmax, np.max( Y1[t][0,:,:] ) )

    print( vmin )
    print( vmax )

    if( PlotTranspose ):
        vmin = np.min( np.abs( Y1[:][0,:,:] - Y1[:][0,:,:].T ) )
        vmax = np.max( np.abs( Y1[:][0,:,:] - Y1[:][0,:,:].T ) )
        if( vmin == 0.0 ): vmin = 1.0e-17

    Norm = None
    if( UseSemiLogYScale ):
        from matplotlib.colors import LogNorm, SymLogNorm

        if( vmin < 0.0 ):

            Norm = SymLogNorm( linthresh = 1.0e2, base = 10 )

        else:

            Norm = LogNorm()

    im = plt.imshow( f(0), \
                     cmap     = cmap, \
                     animated = True, \
                     vmin     = vmin, \
                     vmax     = vmax, \
                     extent   = [ x1.min(), x1.max(), x2.min(), x2.max() ], \
                     aspect   = 'equal', \
                     origin   = 'lower', \
                     norm     = Norm )

    plt.colorbar(im)

    Width  = x1.max() - x1.min()
    Height = x2.max() - x2.min()

    time_text = plt.text( x1.min() + 0.5 * Width, x2.min() + 0.9 * Height, '' )

    def UpdateFrame(t):
        im.set_array( f(t) )
        time_text.set_text('time = {:.3e} {:}'.format( Time[t], TimeUnit ) )
        return im,

    # Call the animator
    anim = animation.FuncAnimation( fig, UpdateFrame, frames = nFiles, \
                                    interval = 100, blit=True)

    #plt.show()
    anim.save( SaveFileAs, fps = 5 )
    plt.close()

os.system( 'rm -f *.pyc' )
os.system( 'rm -rf __pycache__' )
