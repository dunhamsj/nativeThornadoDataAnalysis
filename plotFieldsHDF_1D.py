#!/usr/bin/env python3

import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
plt.style.use( './publication.sty' )

from UtilitiesModuleHDF import ReadFieldsHDF, GetNorm

print( '' )
print( '  Running PlotFieldsHDF_1D.py...' )
print( '  ------------------------------\n' )

# --- Get user's HOME directory ---

HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---

THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True)
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

############################ User Input ############################

Verbose = True

RootPath = THORNADO_DIR + 'SandBox/StandingAccretionShock_CFA_Perturbed/'
suffix = 'Output/'

#Problem = 'Advection'
#Problem = 'RiemannProblem'
#Problem = 'RiemannProblemSpherical'
Problem = 'StandingAccretionShock'
#Problem = 'StaticTOV'
#Problem = 'DynamicTOV'
#Problem = 'YahilCollapse'
#Problem = 'GravitationalCollapse'

Field            = 'PF_D'
UseLogScale      = True
UsePhysicalUnits = True

UseCustomXLabel = True
xLabel = r'$x^{1}\,\left[\mathrm{km}\right]$'

UseCustomYLabel = True
yLabel = r'$\rho\,\left[\mathrm{cgs}\right]$'

Snapshots         = [0,8]
UseGeometryFields = True
yScale            = 1.0

UseCustomLimits = False
ymin = 0.0
ymax = 1.0

FigTitle = Problem

############################ End of User Input #####################

if Verbose:
    print( '    PathToData: {:}'.format( RootPath + suffix ) )
    print( '       Problem: {:}'.format( Problem ) )
    print( '         Field: {:}'.format( Field ) )

Snapshots = np.array( Snapshots )

CoordinateSystem = 'CARTESIAN'

if    Problem == 'RiemannProblemSpherical' \
   or Problem == 'StandingAccretionShock' \
   or Problem == 'StaticTOV' \
   or Problem == 'YahilCollapse' \
   or Problem == 'GravitationalCollapse':

    CoordinateSystem = 'SPHERICAL'

if UsePhysicalUnits:

    TimeUnit   = 'ms'
    LengthUnit = 'km'

else:

    TimeUnit   = ''
    LengthUnit = ''

PathToData = RootPath + suffix + Problem

nFiles = Snapshots.shape[0]

names \
  = ReadFieldsHDF \
      ( PathToData, Snapshots, \
        CoordinateSystem = CoordinateSystem, \
        UsePhysicalUnits = UsePhysicalUnits, \
        UseGeometryFields = UseGeometryFields )

T  = names['Time'][1]
X1 = names['X1'][1]
X2 = names['X2'][1]
X3 = names['X3'][1]

TimeUnit = names['Time'][0]
X1Unit   = names['X1'][0]
X2Unit   = names['X2'][0]
DataUnit = names[Field][0]

assert ( X1.shape[0] > 1 and X2.shape[0] == 1 and X3.shape[0] == 1 ), \
       'must have nDims = 1'

Data = np.empty( (Snapshots.shape[0],X3.shape[0],X2.shape[0],X1.shape[0]), \
                 np.float64 )

for i in range( nFiles ):

    Data[i] = names[Field][1][i,0,0,:]

##### Plotting #####

fig = plt.figure()# figsize = (16,9) )
ax  = fig.add_subplot( 111 )
ax.set_title( r'$\texttt{{{:}}}$'.format( FigTitle ) )

if UseCustomLimits:

    ymin = 0.0
    ymax = 2.0

else:

    ymin = Data.min()
    ymax = Data.max()

Norm = GetNorm( UseLogScale, Data, vmin = ymin, vmax = ymax )

for i in range( nFiles ):

    ax.plot( X1, Data[i,0,0], '.',label = 'Snapshot {:}'.format( Snapshots[i] ) )

ax.legend()


#xRef = [ 8.0e1, 7.0e1, 6.0e1, 5.5e1, 5.0e1, 4.5e1 ]
#for xx in xRef:
#    ax.axvline( xx, c = 'b' )

if UseLogScale:
    ax.set_yscale( 'log' )
ax.grid()
if UseCustomXLabel:
    ax.set_xlabel( xLabel )
else:
    ax.set_xlabel( 'X1' + ' ' + X1Unit )
if UseCustomYLabel:
    ax.set_ylabel( yLabel )
else:
    ax.set_ylabel( Field + ' ' + DataUnit )

plt.show()

#figName = 'fig.{:}_NativeThornado_{:}.png'.format( Problem, Field )
#plt.savefig( figName, dpi = 300 )
#print( '\n  Saved {:}'.format( figName ) )

plt.close()

os.system( 'rm -rf __pycache__' )
