#! /usr/bin/python3 -B
from helperFuncs import makeOperators

# Set list of elements and pairs for which to compile forms
#
# 0: marker
# 1: name of element for test space
# 2: order of element for test space
# 3: name of element for trial space
# 4: order of element for trial space
#
elemList = list()
elemList.append( [ "p1", "Lagrange", 1, "Lagrange", 1 ] )
elemList.append( [ "p2", "Lagrange", 2, "Lagrange", 2 ] )
elemList.append( [ "p1_to_p2", "Lagrange", 2, "Lagrange", 1 ] )
elemList.append( [ "p2_to_p1", "Lagrange", 1, "Lagrange", 2 ] )
elemList.append( [ "n1Curl", "N1curl", 1, "N1curl", 1 ] )

# For each element and pair set the forms to compile
forms = {}
forms[ "p1" ] = [ "diffusion", "div_tet", "divt_tet", "mass", "pspg_tet", "stokes_epsilon_tet", "stokes_full_tet" ]
forms[ "p2" ] = [ "diffusion", "div_tet", "divt_tet", "mass", "pspg_tet", "stokes_epsilon_tet", "stokes_full_tet" ]
forms[ "p1_to_p2" ] = [ "divt_tet" ]
forms[ "p2_to_p1" ] = [ "div_tet" ]
forms[ "n1Curl" ] = [ "mass" ]

# For each element and pair set the directory into which we place the
# generated C++ header files (relative to hytegsrc below)
outDir = {}
outDir[ "p1" ] = [ "forms/form_fenics_generated" ]
outDir[ "p2" ] = [ "forms/form_fenics_generated" ]
outDir[ "p1_to_p2" ] = [ "forms/form_fenics_generated" ]
outDir[ "p2_to_p1" ] = [ "forms/form_fenics_generated" ]
outDir[ "n1Curl" ] = [ "forms/form_fenics_generated" ]

# Set output directory
hytegsrc = "../../src/hyteg"

# Decide whether script should be talkative
beVerbose = True

# Set options for invocation of FFC
ffcOpts = ""

# Set to false to not remove temporary files
cleanTemporaries = True

# ----------------------------------------------------------------------------

makeOperators( elemList, forms, outDir, hytegsrc, beVerbose, ffcOpts,
               cleanTemporaries, "3D" )

# ----------------------------------------------------------------------------
