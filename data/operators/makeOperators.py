#! /usr/bin/python3
import subprocess as sproc

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
elemList.append( [ "bubble", "Bubble", 3, "Bubble", 3 ] )
elemList.append( [ "bubble_to_p1", "Bubble", 3, "Lagrange", 1 ] )
elemList.append( [ "p1_to_bubble", "Lagrange", 1, "Bubble", 3 ] )
elemList.append( [ "p1_to_p2", "Lagrange", 1, "Lagrange", 2 ] )
elemList.append( [ "p2_to_p1", "Lagrange", 2, "Lagrange", 1 ] )

# For each element and pair set the forms to compile
forms = {}
forms[ "p1" ] = [ "diffusion", "div", "divt", "mass", "pspg", "stokes_epsilon" ]
forms[ "p2" ] = [ "diffusion", "div", "divt" ]
forms[ "bubble" ] = [ "diffusion" ]
forms[ "bubble_to_p1" ] = [ "divt" ]
forms[ "p1_to_bubble" ] = [ "div" ]
forms[ "p1_to_p2" ] = [ "divt" ]
forms[ "p2_to_p1" ] = [ "div" ]

# For each element and pair set the directory into which we place the
# generated C++ header files (relative to tinyHHGsrc below)
outDir = {}
outDir[ "p1" ] = [ "p1functionspace/generated" ]
outDir[ "p2" ] = [ "p2functionspace/generated" ]
outDir[ "bubble" ] = [ "bubblefunctionspace/generated" ]
outDir[ "bubble_to_p1" ] = [ "mixedoperators/BubbleToP1/generated" ]
outDir[ "p1_to_bubble" ] = [ "mixedoperators/P1ToBubble/generated" ]
outDir[ "p1_to_p2" ] = [ "mixedoperators/generated" ]
outDir[ "p2_to_p1" ] = [ "mixedoperators/generated" ]

# Set output directory
tinyHHGsrc = "../../src/tinyhhg_core"

# Decide whether script should be talkative
beVerbose = True

# Set options for invocation of FFC
ffcOpts = ""

# Set to false to not remove temporary files
cleanTemporaries = True

# ----------------------------------------------------------------------------

def getFileHandle( fname, operation ):
    try:
        fHandle = open( fname, operation )
    except IOError:
        print( "Sorry could not open file '" + fname + "'" )
    return fHandle

def closeFile( handle, fname ):
    try:
        handle.close()
    except IOError:
        print( "Sorry problems closing file '" + fname + "'" )
        
# ----------------------------------------------------------------------------

def getFormFileName( form ):
    return "form_" + form + ".ufl"

def getFFCFileName( elem, form ):
    return elem + "_" + form + ".ufl"

def getHeaderFileName( elem, form ):
    return elem + "_" + form + ".h"

def readFormFile( form ):
    fname = getFormFileName( form )
    fh = getFileHandle( fname, "r" )
    formDescr = fh.read()
    return formDescr

def genFFCFile( eDescr, form ):

    # Disassemble element description
    eName       = eDescr[0]
    eTypeTest   = eDescr[1]
    eOrderTest  = eDescr[2]
    eTypeTrial  = eDescr[3]
    eOrderTrial = eDescr[4]
  
    # Read (template) file with form description
    formDescr = readFormFile( form )

    # Open output file
    ofName = getFFCFileName( eName, form )
    of = getFileHandle( ofName, "w" )

    # Write comment line
    comment = "# generate code with ffc " + ofName + "\n"
    of.write( comment )

    # Assemble UFL element description (test space)
    eLine = "testElement = FiniteElement( \"" + eTypeTest
    eLine += "\", \"triangle\", " + str(eOrderTest) + " )\n"
    of.write( eLine )

    # Assemble UFL element description (trial space)
    eLine = "trialElement = FiniteElement( \"" + eTypeTrial
    eLine += "\", \"triangle\", " + str(eOrderTrial) + " )\n"
    of.write( eLine )

    # Finish script by adding form description
    of.write( formDescr )
    closeFile( of, ofName )

    # We are done
    return ofName
    
def fixIncludeAndWriteOperatorFile( eDescr, form ):

    # Disassemble element description
    eName  = eDescr[0]

    # Determine name of operator file (input)
    # and destination for operator file (output)
    inFileName  = getHeaderFileName( eName, form )
    # outFileName = tinyHHGsrc + "/" + eName + "functionspace/generated/"
    # outFileName += inFileName
    outFileName = tinyHHGsrc + "/" + outDir[ eName ][0] + "/" + inFileName

    # Read in header file
    with open( inFileName, "r" ) as inFile :
        sourceCode = inFile.read()

    # Adapt include directive
    sourceCode = sourceCode.replace( "#include <ufc.h>",
                                     "#include \"tinyhhg_core/fenics/ufc.h\"" )

    # Write header file
    with open( outFileName, "w" ) as outFile :
        outFile.write( sourceCode )

# ----------------------------------------------------------------------------

def logMsg( mesg ):
    if beVerbose:
        print( mesg )

# ----------------------------------------------------------------------------

for elem in elemList:

    # Be verbose
    eName = elem[0]
    if elem[1] == elem[3] and elem[2] == elem[4]:
        logMsg( " Generating forms for element type '" + eName + "'" )
    else:
        logMsg( " Generating forms for mixed element pair '" + eName + "'" )
        
    for form in forms[ eName ]:
        logMsg( " -> current form: " + form )
        ffcFileName = genFFCFile( elem, form )
        command = [ "ffc" ] + ffcOpts.split() + [ ffcFileName ]
        sproc.run( command )
        fixIncludeAndWriteOperatorFile( elem, form )

        if cleanTemporaries:
            sproc.run( [ "rm", getFFCFileName( elem[0], form ) ] )
            sproc.run( [ "rm", getHeaderFileName( elem[0], form ) ] )

# ----------------------------------------------------------------------------
