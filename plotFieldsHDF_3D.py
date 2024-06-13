#!/usr/bin/env python3

import numpy as np
import os

from UtilitiesModuleHDF import ReadFieldsHDF

RootPath  = '/home/kkadoogan/'
Problem   = 'KelvinHelmholtzInstability'
Field     = 'PF_D'
Snapshots = np.array( [0] )

suffix = 'OutputGPU/'
PathToData = RootPath + suffix + Problem
names = ReadFieldsHDF( PathToData, Snapshots, 'CARTESIAN', False, True )
DataG = np.copy( names[Field][1] )

suffix = 'OutputCPU/'
PathToData = RootPath + suffix + Problem
names = ReadFieldsHDF( PathToData, Snapshots, 'CARTESIAN', False, True )
DataC = np.copy( names[Field][1] )

print( DataG.min(), DataG.max() )
print( DataC.min(), DataC.max() )
diff = np.abs( DataG - DataC ) / ( 0.5 * np.abs( DataG + DataC ) )

print( diff.min(), diff.max() )

os.system( 'rm -rf __pycache__' )
