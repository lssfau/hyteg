import subprocess as sproc

areWeVerbose = False

# ----------------------------------------------------------------------------

def makeOperators( elemList, forms, outDir, hytegsrc, beVerbose, ffcOpts,
                   cleanTemporaries, dimension ):

    global areWeVerbose
    areWeVerbose = beVerbose

    # Check dimension
    if dimension == "2D":
        geometricObject = "triangle"
        dimPrefix = ""
    elif dimension == "3D":
        geometricObject = "tetrahedron"
        dimPrefix = "tet_"
    else:
        print( "ERROR: I don't understand dimension '" + dimension )
        return

    for elem in elemList:

        # Be verbose
        eName = elem[0]
        if elem[1] == elem[3] and elem[2] == elem[4]:
            logMsg( " Generating forms for element type '" + eName + "'" )
        else:
            logMsg( " Generating forms for mixed element pair '" + eName + "'" )
        
        for form in forms[ eName ]:
            logMsg( " -> current form: " + form )
            ffcFileName = genFFCFile( elem, form, geometricObject, dimPrefix )
            command = [ "ffc" ] + ffcOpts.split() + [ ffcFileName ]
            sproc.run( command )
            fixIncludeAndWriteOperatorFile( elem, form, hytegsrc, outDir,
                                            dimPrefix )

            if cleanTemporaries:
                sproc.run([ "rm", getFFCFileName( elem[0], form, dimPrefix )])
                sproc.run([ "rm", getHeaderFileName( elem[0], form, dimPrefix)])

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

def getFFCFileName( elem, form, dimPrefix ):
    return elem + "_" + dimPrefix + form + ".ufl"

def getHeaderFileName( elem, form, dimPrefix ):
    return elem + "_" + dimPrefix + form + ".h"

def readFormFile( form ):
    fname = getFormFileName( form )
    fh = getFileHandle( fname, "r" )
    formDescr = fh.read()
    return formDescr

def genFFCFile( eDescr, form, geometricObject, dimPrefix ):

    # Disassemble element description
    eName       = eDescr[0]
    eTypeTest   = eDescr[1]
    eOrderTest  = eDescr[2]
    eTypeTrial  = eDescr[3]
    eOrderTrial = eDescr[4]
  
    # Read (template) file with form description
    formDescr = readFormFile( form )

    # Open output file
    ofName = getFFCFileName( eName, form, dimPrefix )
    of = getFileHandle( ofName, "w" )

    # Write comment line
    comment = "# generate code with ffc " + ofName + "\n"
    of.write( comment )

    # Assemble UFL element description (test space)
    eLine = "testElement = FiniteElement( \"" + eTypeTest
    eLine += "\", \"" + geometricObject + "\", " + str(eOrderTest) + " )\n"
    of.write( eLine )

    # Assemble UFL element description (trial space)
    eLine = "trialElement = FiniteElement( \"" + eTypeTrial
    eLine += "\", \"" + geometricObject + "\", " + str(eOrderTrial) + " )\n"
    of.write( eLine )

    # Finish script by adding form description
    of.write( formDescr )
    closeFile( of, ofName )

    # We are done
    return ofName

# ----------------------------------------------------------------------------
    
def fixIncludeAndWriteOperatorFile( eDescr, form, hytegsrc, outDir,
                                    dimPrefix ):

    # Disassemble element description
    eName  = eDescr[0]

    # Determine name of operator file (input)
    # and destination for operator file (output)
    inFileName  = getHeaderFileName( eName, form, dimPrefix )
    # outFileName = hytegsrc + "/" + eName + "functionspace/generated/"
    # outFileName += inFileName
    outFileName = hytegsrc + "/" + outDir[ eName ][0] + "/" + inFileName

    # Read in header file
    with open( inFileName, "r" ) as inFile :
        sourceCode = inFile.read()

    # Adapt include directive
    sourceCode = sourceCode.replace( "#include <ufc.h>",
                                     "#include \"hyteg/fenics/ufc.h\"" )

    # Write header file
    with open( outFileName, "w" ) as outFile :
        outFile.write( sourceCode )

# ----------------------------------------------------------------------------

def logMsg( mesg ):
    global areWeVerbose
    if areWeVerbose:
        print( mesg )

# ----------------------------------------------------------------------------
