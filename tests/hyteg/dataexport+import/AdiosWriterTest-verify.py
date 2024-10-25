# -*- coding: utf-8 -*-

import sys
import subprocess as sproc
import argparse

# Where to find output files to compare with each other
dataDir = "AdiosWriterTest-Output"
refsDir = "AdiosWriterTest-References"

# Base names of output files (TAG will be replaced by adaptTestCaseNames())
testCases = [
    "AdiosWriterTest_2D-TAG-P1_level3",
    "AdiosWriterTest_2D-TAG-P2_level3",
    "AdiosWriterTest_3D-TAG-P1_level2",
    "AdiosWriterTest_3D-TAG-P2_level2",
    "AdiosWriterTest_3DSurface-TAG-P1_level3",
    "AdiosWriterTest_3DSurface-TAG-P2_level3" ]


def initArgparse() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(
        usage="%(prog)s -nw NUM_WRITERS -fp FP_BITS",
        description="Compare BP-files created with AdiosWriterTest to reference files."
    )

    parser.add_argument( "-n", "--num-writers", required=True, help="number of writers used to create files" )
    parser.add_argument( "-f", "--fp-bits", required=True, help="number of bits in FP-type of files" )

    return parser


def adaptTestCaseNames( testCases, args ):
    """Generate tags from CLI input and insert it into test-case names"""

    tag = f"nw={args.num_writers}-fp={args.fp_bits}"

    realTestCases = []
    for case in testCases:
        realTestCases.append( case.replace( "TAG", tag ) )

    return realTestCases


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


def importConnectivityData( case ):
    """Read connectivity info from reference file for specific case"""

    # remove fp part from case name
    filename = case + "-connectivity.dump"
    filename = filename.replace( "-fp=32", "" )
    filename = filename.replace( "-fp=64", "" )

    dumpFile = open( refsDir + "/" + filename, "r" )
    connectivity = dumpFile.read().splitlines()
    return connectivity


def compareData( dataTest, dataRef ):
    """Compare job created bpls output line-by-line with reference output"""

    everythingIsFine = True

    if len( dataTest ) > len( dataRef ):
        print( " ERROR: test data has more lines than reference!" )
        print( "        test data follows:" )
        for elem in dataTest:
            print( elem )
        everythingIsFine = False
    elif len( dataTest ) < len( dataRef ):
        print( " ERROR: test data has fewer lines than reference!" )
        print( "        test data follows:" )
        for elem in dataTest:
            print( elem )
        everythingIsFine = False
    else:
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

def dumpConnectivity( case ):
    """Dump connectivity from BP files created by test job using bpls and filter it"""

    tool = "bpls"          # must be in path (but building with ADIOS2 should take care of this)
    options ="-d"          # dump a field
    field = "connectivity"

    command = tool + " " + options + " " + dataDir + "/" + case + ".bp" + " " + field
    bplsOutput = sproc.getoutput( command )
    bplsOutput = bplsOutput.splitlines()

    return bplsOutput

def executeSingleTest( case ):
    """Compare data of single BP file to reference(s)"""
    print( "===============================================================" )
    print( f" Comparing {case}.bp to references" )

    status = executeMetaDataTest( case )
    status = status and executeConnectivityDataTest( case )

    return status

def executeMetaDataTest( case ):
    """Compare single BP file to reference"""

    print( " * checking meta-data", end="" )

    dataJob = runBPLSonBPfilesFromJob( case )
    dataRef = importReferenceData( case )
    checkPassed = compareData( dataJob, dataRef )
    if checkPassed:
        print( " -> OKAY" )
    else:
        print( " -> FAILED" )

    return checkPassed

def executeConnectivityDataTest( case ):
    """Compare connectivity info from single BP file to reference"""

    print( " * checking connectivity", end="" )

    dataJob = dumpConnectivity( case )
    dataRef = importConnectivityData( case )
    checkPassed = compareData( dataJob, dataRef )
    if checkPassed:
        print( " -> OKAY" )
    else:
        print( " -> FAILED" )

    return checkPassed

def executeFailingTest( case ):
    """Modify test-data so that check fails"""

    print( "===============================================================" )
    print( f" Comparing {case}.bp to reference" )
    print( " ---------" )

    # dataJob = runBPLSonBPfilesFromJob( case )
    dataJob = importReferenceData( case )

    replacementLine = "  int64_t   NumberOfElements          {1} = <error here!>"
    for idx, line in enumerate( dataJob ):
        if "int64_t   NumberOfElements" in line:
            dataJob[idx] = replacementLine

    dataRef = importReferenceData( case )
    checkPassed = compareData( dataJob, dataRef )
    if checkPassed:
        print( " CHECK: okay --> ERROR not detected!!!!" )
    else:
        print( " CHECK: failed --> GOOD: error was detected!" )
    return not checkPassed

# -----------------------
#  check all test cases
# -----------------------

# be pessimistic
allChecksPassed = False

# we need the number of processes and the fp-type
# to be handed to the script as input
parser = initArgparse()
args = parser.parse_args()
print( f"num_writers = {args.num_writers},  fp_bits = {args.fp_bits}" );

testCases = adaptTestCaseNames( testCases, args )

for case in testCases:
    allChecksPassed = True if executeSingleTest( case ) else False

executeFailingTest( testCases[0] )
allChecksPassed = True if executeSingleTest( case ) else False

print( "===============================================================" )

if allChecksPassed:
    sys.exit(0)
else:
    sys.exit(1)
