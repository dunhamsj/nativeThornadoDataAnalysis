#!/usr/bin/env python3

import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
plt.style.use( './publication.sty' )

from UtilitiesModuleHDF import ReadFieldsHDF, GetNorm

print( '' )
print( 'Running PlotFieldsHDF_3D.py...' )
print( '------------------------------' )

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True)
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

############################ User Input ############################

RootPath = THORNADO_DIR + 'SandBox/dgExperiments_Euler_Relativistic_IDEAL/'
suffix = 'Output/'

RootPath = '/home/kkadoogan/'
suffix = 'KHI3D_GPU/'

Problem = 'KelvinHelmholtzInstability'

xL = np.array( [ -0.5, -1.0, -0.5 ], np.float64 )
xU = np.array( [ +0.5, +1.0, +0.5 ], np.float64 )

Field            = 'PF_D'
UseLogScale      = False
UsePhysicalUnits = False

Snapshots         = [0,43]
cmap              = 'jet'
PlotTranspose     = False
PlotContours      = False
UseGeometryFields = True
yScale            = 1.0

UseCustomLimits = True
vmin = 0.0
vmax = 2.0

FigTitle = Problem

############################ End of User Input #####################

Snapshots = np.array( Snapshots )

Polar            = False
CoordinateSystem = 'CARTESIAN'

if Problem == 'StandingAccretionShock':

    Polar            = True
    CoordinateSystem = 'SPHERICAL'

if UsePhysicalUnits:

    TimeUnit   = 'ms'
    LengthUnit = 'km'

else:

    TimeUnit   = ''
    LengthUnit = ''

PathToData = RootPath + suffix + Problem

nFiles = Snapshots.shape[0]

names = ReadFieldsHDF( PathToData, Snapshots, CoordinateSystem, \
                       UsePhysicalUnits, UseGeometryFields )

T  = names['Time'][1]
X1 = names['X1'][1]
X2 = names['X2'][1]
X3 = names['X3'][1]

TimeUnit = names['Time'][0]
X1Unit   = names['X1'][0]
X2Unit   = names['X2'][0]
DataUnit = names[Field][0]

Data = np.empty( (Snapshots.shape[0],X3.shape[0],X2.shape[0],X1.shape[0]), \
                 np.float64 )

Data = names[Field][1]

##### Plotting #####

fig = plt.figure( figsize = (16,9) )
ax  = fig.add_subplot( 111, polar = Polar )
fig.suptitle( FigTitle )

if UseCustomLimits:

    vmin = 0.0
    vmax = 2.0

else:

    vmin = Data.min()
    vmax = Data.max()

Norm = GetNorm( UseLogScale, Data, vmin = vmin, vmax = vmax )

extent = [ xL[0], xU[0], xL[1], xU[1] ]
aspect = ( extent[1] - extent[0] ) / ( extent[3] - extent[2] )

if CoordinateSystem == 'SPHERICAL':

    """
    Taken from:
    https://brushingupscience.com/2016/06/21/
    matplotlib-animations-the-easy-way/
    """

    X2v, X1v = np.meshgrid( X2, X1 )

    im = ax.pcolormesh( X2v, X1v, Data[0,0], \
                        cmap    = cmap, \
                        norm    = Norm, \
                        shading = 'nearest' )

    ax.set_thetamin( 180.0/np.pi * xL[1] )
    ax.set_thetamax( 180.0/np.pi * xU[1] )

    ax.set_theta_zero_location( 'W' ) # z-axis horizontal
    ax.set_theta_direction( -1 )

    ax.text( 0.4, 0.9, 'Time = {:.2e} {:}'.format \
             ( T[iSS], TimeUnit ), \
             transform = ax.transAxes )

else:

    iSS = -1
    iX3 = -1

    X1v, X2v = np.meshgrid( X1, X2 )

    Data -= Data[0]

    if PlotContours:

        levels = np.linspace( vmin, vmax, 30 )

        im = ax.contour( Data[iSS,iX3], \
                         levels = levels, \
                         extent = extent, \
                         cmap   = cmap )

    else:

        im = ax.pcolormesh( X1v, X2v, Data[iSS,iX3], \
                            cmap    = cmap, \
                            norm    = Norm, \
                            shading = 'nearest' )

    ax.text( 0.4, 0.9, 'Time = {:.2e} {:}'.format \
             ( T[0], TimeUnit ), \
             transform = ax.transAxes )

    ax.set_xlabel( 'X1' + ' ' + X1Unit )
    ax.set_ylabel( 'X2' + ' ' + X2Unit )

cbar = fig.colorbar( im )
cbar.set_label( Field + ' ' + DataUnit )

plt.show()
plt.close()

os.system( 'rm -rf __pycache__' )
