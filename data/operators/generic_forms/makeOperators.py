#! /usr/bin/python3
import subprocess as sproc

# Set list of elements for which to compile forms
elemList = list()
elemList.append( [ "p1", "Lagrange", 1 ] )
elemList.append( [ "p2", "Lagrange", 2 ] )
elemList.append( [ "bubble", "Bubble", 3 ] )

# For each element set the forms to compile
forms = {}
forms[ "p1" ] = [ "diffusion", "div", "divt", "mass", "pspg", "stokes_epsilon" ]
forms[ "p2" ] = [ "diffusion", "div", "divt" ]
forms[ "bubble" ] = [ "diffusion" ]

# Decide whether script should be talkative
beVerbose = True

# Set options for invocation of FFC
ffcOpts = ""

# Set output directory
tinyHHGsrc = "../../../src/tinyhhg_core"

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

def readFormFile( elem, form ):
    fname = getFormFileName( form )
    fh = getFileHandle( fname, "r" )
    formDescr = fh.read()
    return formDescr

def genFFCFile( eDescr, form ):

    # Disassemble element description
    eName  = eDescr[0]
    eType  = eDescr[1]
    eOrder = eDescr[2]
  
    # Read (template) file with form description
    formDescr = readFormFile( eName, form )

    # Open output file
    ofName = getFFCFileName( eName, form )
    of = getFileHandle( ofName, "w" )

    # Write comment line
    comment = "# generate code with ffc " + ofName + "\n"
    of.write( comment )

    # Assemble UFL element description
    eLine = "element = FiniteElement( \"" + eType + "\", \"triangle\", "
    eLine += str(eOrder) + " )\n"
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
    outFileName = tinyHHGsrc + "/" + eName + "functionspace/generated/"
    outFileName += inFileName

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
    logMsg( " Generating forms for element type '" + eName + "'" )

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
