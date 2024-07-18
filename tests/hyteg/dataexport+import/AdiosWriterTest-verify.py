# -*- coding: utf-8 -*-

import sys
import subprocess as sproc

# Where to find output files to compare with each other
dataDir = "AdiosWriterTest-Output"
refsDir = "AdiosWriterTest-References"

# Base names of output files
testCases = [
    "AdiosWriterTest_2D-P1_level3",
    "AdiosWriterTest_2D-P2_level3",
    "AdiosWriterTest_3D-P1_level2",
    "AdiosWriterTest_3D-P2_level2",
    "AdiosWriterTest_3DSurface-P1_level3",
    "AdiosWriterTest_3DSurface-P2_level3" ]


def filterOutVolatileAttribute( bplsOutput ):
    """Remove attributes whose values are naturally changing between jobs

    BP files contain attributes which are 'volatile' in the sense that they
    depend on job specifics such as e.g. compiler, build type or MPI library.
    This function removes the corresponding lines from bplsOutput."""

    volatileAttrsKeys = [ "git branch", "SHA1 of last commit", "build type", "compiler", "compiler flags", "mpi version" ]

    filtered = []

    for line in bplsOutput:
        if not line.startswith( tuple( volatileAttrsKeys ) ):
            filtered.append( line )

    return filtered


def importReferenceData( case ):
    """Read and filter bpls output from reference file for specific case"""
    refFile = open( refsDir + "/" + case + ".txt", "r" )
    refOutput = refFile.read().splitlines()
    return filterOutVolatileAttribute( refOutput )


def compareData( dataTest, dataRef ):
    """Compare job created bpls output line-by-line with reference output"""

    everythingIsFine = True

    for index, elemTest in enumerate( dataTest ):
        elemRef = dataRef[ index ]
        if elemTest != elemRef:
            print( " ERROR: difference to reference detected:" )
            print( " * reference value --> " + elemRef )
            print( " * value from test --> " + elemTest )
            everythingIsFine = False

    return everythingIsFine

def runBPLSonBPfilesFromJob( case ):
    """Extract info from BP files created by test job using bpls and filter it"""

    tool = "bpls"   # must be in path (but building with ADIOS2 should take care of this)
    options ="-la"  # list summary info and include attributes

    command = tool + " " + options + " " + dataDir + "/" + case + ".bp"
    bplsOutput = sproc.getoutput( command )
    bplsOutput = bplsOutput.splitlines()

    return filterOutVolatileAttribute( bplsOutput )

def executeSingleTest( case ):
    """Compare single BP file to reference"""

    print( "---------------------------------------------------------------" )
    print( f" Comparing {case}.bp to reference" )
    print( "---------------------------------------------------------------" )

    dataJob = runBPLSonBPfilesFromJob( case )
    dataRef = importReferenceData( case )
    checkPassed = compareData( dataJob, dataRef )
    return checkPassed

# -----------------------
#  check all test cases
# -----------------------

# be pessimistic
allChecksPassed = False

for case in testCases:
    allChecksPassed = True if executeSingleTest( case ) else False

if allChecksPassed:
    sys.exit(0)
else:
    sys.exit(1)
