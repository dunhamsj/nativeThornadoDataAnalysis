#!/usr/bin/env python3

# --- Import libraries ---

import warnings
warnings.simplefilter( action = 'ignore' , category = FutureWarning )

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

import numpy as np
import h5py as h5

def GetNorm( UseLogScale, Data, vmin = +1.0e100, vmax = -1.0e100 ):

    if vmin > +1.0e99: vmin = Data.min()
    if vmax < -1.0e99: vmax = Data.max()

    if( UseLogScale ):

        if( np.any( Data <= 0.0 ) ):

            Norm = SymLogNorm( vmin = vmin, vmax = vmax, \
                               linthresh = 1.0e-2, base = 10 )

        else:

            Norm = LogNorm   ( vmin = vmin, vmax = vmax )

    else:

        Norm = plt.Normalize ( vmin = vmin, vmax = vmax )

    return Norm

# End of GetNorm

def ReadFieldsHDF( PathToData, Snapshots, CoordinateSystem, \
                   UsePhysicalUnits, UseGeometryFields = True ):

    FF_root = PathToData + '_FluidFields_'
    GF_root = PathToData + '_GeometryFields_'

    nFiles = len( Snapshots )

    TimeSteps              = np.empty( nFiles, object )
    DataFileNames_Fluid    = np.empty( nFiles, object )
    DataFileNames_Geometry = np.empty( nFiles, object )

    # Arrays to hold data
    Data_FF = np.empty( nFiles, object )
    Data_GF = np.empty( nFiles, object )

    for i in range( 1 ):

        TimeSteps[i] = str( Snapshots[i] ).zfill( 6 )

        DataFileNames_Fluid[i] = FF_root + str( TimeSteps[i] ) + '.h5'
        Data_FF[i]             = h5.File( DataFileNames_Fluid[i], 'r' )

        if( UseGeometryFields ):
            DataFileNames_Geometry[i] \
              = GF_root + str( TimeSteps[i] ) + '.h5'
            Data_GF[i] \
              = h5.File( DataFileNames_Geometry[i], 'r' )

    # Get the spatial grid
    X    = Data_FF[0][ 'Spatial Grid' ]
    X1   = np.array( X[ 'X1' ] )
    X2   = np.array( X[ 'X2' ] )
    X3   = np.array( X[ 'X3' ] )
    X1_C = np.array( X[ 'X1_C' ] )
    X2_C = np.array( X[ 'X2_C' ] )
    X3_C = np.array( X[ 'X3_C' ] )

    T  = np.empty( nFiles, float  )
    FF = np.empty( nFiles, object )
    PF = np.empty( nFiles, object )
    CF = np.empty( nFiles, object )
    AF = np.empty( nFiles, object )
    DF = np.empty( nFiles, object )
    GF = np.empty( nFiles, object )

    # Conserved fields

    CF_D  = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    CF_S1 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    CF_S2 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    CF_S3 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    CF_E  = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    CF_Ne = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    # Primitive fields

    PF_D  = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    PF_V1 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    PF_V2 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    PF_V3 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    PF_E  = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    PF_Ne = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    # Auxiliary fields

    AF_P  = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    AF_Cs = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    AF_Gm = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    AF_Ye = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    # Diagnostic fields

    DF_TCI   \
      = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    DF_Sh_X1 \
      = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    DF_Sh_X2 \
      = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    DF_Sh_X3 \
      = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    # Geometry fields

    GF_CF = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_al = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_NP = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_b1 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_b2 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_b3 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_g1 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_g2 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_g3 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_h1 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_h2 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_h3 = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )
    GF_Sg = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    # Derived fields

    LorentzFactor \
      = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    SpecificEnthalpy \
      = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    PolytropicConstant \
      = np.empty( (nFiles,X3.shape[0],X2.shape[0],X1.shape[0]), np.float64 )

    # Set units

    if UsePhysicalUnits:

        TimeUnit               = 'ms'
        MassDensityUnit        = 'g/cm**3'
        EnergyDensityUnit      = 'erg/cm**3'
        NumberDensityUnit      = '1/cm**3'
        PressureUnit           = 'erg/cm**3'
        PolytropicConstantUnit = 'erg/cm**3/(g/cm**3)**Gamma'

        if CoordinateSystem == 'CARTESIAN':

            X1Unit                = 'km'
            X2Unit                = 'km'
            X3Unit                = 'km'
            MomentumDensityX1Unit = 'g/cm**2/s'
            MomentumDensityX2Unit = 'g/cm**2/s'
            MomentumDensityX3Unit = 'g/cm**2/s'
            VelocityX1Unit        = 'km/s'
            VelocityX2Unit        = 'km/s'
            VelocityX3Unit        = 'km/s'
            Gm11Unit              = ''
            Gm22Unit              = ''
            Gm33Unit              = ''
            h1Unit                = ''
            h2Unit                = ''
            h3Unit                = ''
            SqrtGmUnit            = ''

        elif CoordinateSystem == 'CYLINDRICAL':

            X1Unit                = 'km'
            X2Unit                = 'km'
            X3Unit                = 'rad'
            MomentumDensityX1Unit = 'g/cm**2/s'
            MomentumDensityX2Unit = 'g/cm**2/s'
            MomentumDensityX3Unit = 'g/cm/s'
            VelocityX1Unit        = 'km/s'
            VelocityX2Unit        = 'km/s'
            VelocityX3Unit        = 'rad/s'
            Gm11Unit              = ''
            Gm22Unit              = ''
            Gm33Unit              = 'km**2'
            h1Unit                = ''
            h2Unit                = ''
            h3Unit                = 'km'
            SqrtGmUnit            = 'km'

        elif CoordinateSystem == 'SPHERICAL':

            X1Unit                = 'km'
            X2Unit                = 'rad'
            X3Unit                = 'rad'
            MomentumDensityX1Unit = 'g/cm**2/s'
            MomentumDensityX2Unit = 'g/cm/s'
            MomentumDensityX3Unit = 'g/cm/s'
            VelocityX1Unit        = 'km/s'
            VelocityX2Unit        = 'rad/s'
            VelocityX3Unit        = 'rad/s'
            Gm11Unit              = ''
            Gm22Unit              = 'km**2'
            Gm33Unit              = 'km**2'
            h1Unit                = ''
            h2Unit                = 'km'
            h3Unit                = 'km'
            SqrtGmUnit            = 'km**2'

        else:

          print( 'Invalid coordinate system: {:s}'.format( CoordinateSytem ) )
          exit()

    else:

        TimeUnit               = ''
        MassDensityUnit        = ''
        EnergyDensityUnit      = ''
        NumberDensityUnit      = ''
        PressureUnit           = ''
        PolytropicConstantUnit = ''
        X1Unit                 = ''
        X2Unit                 = ''
        X3Unit                 = ''
        MomentumDensityX1Unit  = ''
        MomentumDensityX2Unit  = ''
        MomentumDensityX3Unit  = ''
        VelocityX1Unit         = ''
        VelocityX2Unit         = ''
        VelocityX3Unit         = ''
        Gm11Unit               = ''
        Gm22Unit               = ''
        Gm33Unit               = ''
        h1Unit                 = ''
        h2Unit                 = ''
        h3Unit                 = ''
        SqrtGmUnit             = ''

    for i in range( nFiles ):

        TimeSteps[i] = str( Snapshots[i] ).zfill( 6 )

        DataFileNames_Fluid[i] = FF_root + str( TimeSteps[i] ) + '.h5'

        if( UseGeometryFields ):

            DataFileNames_Geometry[i] = GF_root + str( TimeSteps[i] ) + '.h5'

            with h5.File( DataFileNames_Geometry[i], 'r' ) as f:

                Data_GF[i] = f

                GF[i] = Data_GF[i]['Geometry Fields']

                GF_CF[i] = GF[i][ 'Conformal Factor'                 ][:][:][:]
                GF_al[i] = GF[i][ 'Lapse Function'                   ][:][:][:]
                GF_NP[i] = GF[i][ 'Newtonian Potential'              ][:][:][:]
                GF_b1[i] = GF[i][ 'Shift Vector (1)'                 ][:][:][:]
                GF_b2[i] = GF[i][ 'Shift Vector (2)'                 ][:][:][:]
                GF_b3[i] = GF[i][ 'Shift Vector (3)'                 ][:][:][:]
                GF_g1[i] = GF[i][ 'Spatial Metric Component (11)'    ][:][:][:]
                GF_g2[i] = GF[i][ 'Spatial Metric Component (22)'    ][:][:][:]
                GF_g3[i] = GF[i][ 'Spatial Metric Component (33)'    ][:][:][:]
                GF_h1[i] = GF[i][ 'Spatial Scale Factor (1)'         ][:][:][:]
                GF_h2[i] = GF[i][ 'Spatial Scale Factor (2)'         ][:][:][:]
                GF_h3[i] = GF[i][ 'Spatial Scale Factor (3)'         ][:][:][:]
                GF_Sg[i] = GF[i][ 'Sqrt Spatial Metric Determinant'  ][:][:][:]

        with h5.File( DataFileNames_Fluid[i], 'r' ) as f:

            Data_FF[i] = f

            # First level groups

            FF[i] = Data_FF[i]['Fluid Fields']

            # Second level groups

            PF[i] = FF[i][ 'Primitive' ]
            CF[i] = FF[i][ 'Conserved' ]
            AF[i] = FF[i][ 'Auxiliary' ]
            DF[i] = FF[i][ 'Diagnostic' ]

            # Third level groups

            CF_D [i] = CF[i][ 'Conserved Baryon Density'       ][:][:][:]
            CF_S1[i] = CF[i][ 'Conserved Momentum Density (1)' ][:][:][:]
            CF_S2[i] = CF[i][ 'Conserved Momentum Density (2)' ][:][:][:]
            CF_S3[i] = CF[i][ 'Conserved Momentum Density (3)' ][:][:][:]
            CF_E [i] = CF[i][ 'Conserved Energy Density'       ][:][:][:]
            CF_Ne[i] = CF[i][ 'Conserved Electron Density'     ][:][:][:]

            PF_D [i] = PF[i][ 'Comoving Baryon Density'   ][:][:][:]
            PF_V1[i] = PF[i][ 'Three-Velocity (1)'        ][:][:][:]
            PF_V2[i] = PF[i][ 'Three-Velocity (2)'        ][:][:][:]
            PF_V3[i] = PF[i][ 'Three-Velocity (3)'        ][:][:][:]
            PF_E [i] = PF[i][ 'Internal Energy Density'   ][:][:][:]
            PF_Ne[i] = PF[i][ 'Comoving Electron Density' ][:][:][:]

            AF_P [i] = AF[i][ 'Pressure'                        ][:][:][:]
            AF_Cs[i] = AF[i][ 'Sound Speed'                     ][:][:][:]
            AF_Gm[i] = AF[i][ 'Ratio of Specific Heats (Gamma)' ][:][:][:]
            AF_Ye[i] = AF[i][ 'Electron Fraction'               ][:][:][:]

            DF_TCI  [i] = DF[i][ 'TCI'        ][:][:][:]
            DF_Sh_X1[i] = DF[i][ 'Shock (X1)' ][:][:][:]
            DF_Sh_X2[i] = DF[i][ 'Shock (X2)' ][:][:][:]
            DF_Sh_X3[i] = DF[i][ 'Shock (X3)' ][:][:][:]

            T[i] = Data_FF[i][ 'Time' ][0]

            # Derived fields

            c_cm = 2.99792458e10
            c_km = 2.99792458e5

            BetaSq = (   GF_g1[i] * PF_V1[i]**2 \
                       + GF_g2[i] * PF_V2[i]**2 \
                       + GF_g3[i] * PF_V3[i]**2 ) / c_km**2

            LorentzFactor[i] = 1.0 / np.sqrt( 1.0 - BetaSq )

            SpecificEnthalpy[i] \
              = ( c_cm**2 + ( PF_E[i] + AF_P[i] ) / PF_D[i] ) / c_cm**2

            PolytropicConstant[i] = AF_P[i] / PF_D[i]**( AF_Gm[i] )

    names =\
    { \
      'Time'               : [ TimeUnit              , T        ], \
      'X1'                 : [ X1Unit                , X1       ], \
      'X2'                 : [ X2Unit                , X2       ], \
      'X3'                 : [ X3Unit                , X3       ], \
      'X1_C'               : [ X1Unit                , X1_C     ], \
      'X2_C'               : [ X2Unit                , X2_C     ], \
      'X3_C'               : [ X3Unit                , X3_C     ], \
      'CF_D'               : [ MassDensityUnit       , CF_D     ]  , \
      'CF_S1'              : [ MomentumDensityX1Unit , CF_S1    ], \
      'CF_S2'              : [ MomentumDensityX2Unit , CF_S2    ], \
      'CF_S3'              : [ MomentumDensityX3Unit , CF_S3    ], \
      'CF_E'               : [ EnergyDensityUnit     , CF_E     ], \
      'CF_Ne'              : [ NumberDensityUnit     , CF_Ne    ], \
      'PF_D'               : [ MassDensityUnit       , PF_D     ], \
      'PF_V1'              : [ VelocityX1Unit        , PF_V1    ], \
      'PF_V2'              : [ VelocityX2Unit        , PF_V2    ], \
      'PF_V3'              : [ VelocityX3Unit        , PF_V3    ], \
      'PF_E'               : [ EnergyDensityUnit     , PF_E     ], \
      'PF_Ne'              : [ NumberDensityUnit     , PF_Ne    ], \
      'AF_P'               : [ PressureUnit          , AF_P     ], \
      'AF_Cs'              : [ VelocityX1Unit        , AF_Cs    ], \
      'AF_Gm'              : [ ''                    , AF_Gm    ], \
      'AF_Ye'              : [ ''                    , AF_Ye    ], \
      'GF_CF'              : [ ''                    , GF_CF    ], \
      'GF_al'              : [ ''                    , GF_al    ], \
      'GF_NP'              : [ ''                    , GF_NP    ], \
      'GF_b1'              : [ VelocityX1Unit        , GF_b1    ], \
      'GF_b2'              : [ VelocityX2Unit        , GF_b2    ], \
      'GF_b2'              : [ VelocityX3Unit        , GF_b3    ], \
      'GF_g1'              : [ Gm11Unit              , GF_g1    ], \
      'GF_g2'              : [ Gm22Unit              , GF_g2    ], \
      'GF_g3'              : [ Gm33Unit              , GF_g3    ], \
      'GF_h1'              : [ h1Unit                , GF_h1    ], \
      'GF_h2'              : [ h2Unit                , GF_h2    ], \
      'GF_h3'              : [ h3Unit                , GF_h3    ], \
      'GF_Sg'              : [ SqrtGmUnit            , GF_Sg    ], \
      'DF_TCI'             : [ ''                    , DF_TCI   ], \
      'DF_Sh_X1'           : [ ''                    , DF_Sh_X1 ], \
      'DF_Sh_X2'           : [ ''                    , DF_Sh_X2 ], \
      'DF_Sh_X3'           : [ ''                    , DF_Sh_X3 ], \
      'LorentzFactor'      : [ ''                    , LorentzFactor ], \
      'SpecificEnthalpy'   : [ ''                    , SpecificEnthalpy ], \
      'PolytropicConstant' : [ PolytropicConstantUnit, PolytropicConstant ]
    }

    return names

# End of ReadFieldsHDF
