#include <sstream>
#include <string>
#include <core/timing/Timer.h>
#include <core/Environment.h>
#include <core/math/Constants.h>

// Primitive management
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

// Function spaces
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"

// P1 Operators
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

// P2 Operators
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"

// Mixed Operators
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"

// PETSc interface
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScSparseMatrix.hpp"
#include "tinyhhg_core/petsc/PETScExportOperatorMatrix.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;
using walberla::math::pi;

using namespace hyteg;


typedef enum { P1MASS, P1LAPLACE, P2MASS, P2LAPLACE, P2P1STOKES } operatorTag;

typedef struct {
  operatorTag oprEnum;
  std::string matName;
  bool elimDirichletBC;
} oprInfo;

std::map<std::string, oprInfo> oprMap = {
  {  "P1Mass",     { P1MASS,     "MassOpP1",     false } },
  {  "P1Diff",     { P1LAPLACE,  "DiffOpP1",     true  } },
  {  "P2Mass",     { P2MASS,     "MassOpP2",     false } },
  {  "P2Diff",     { P2LAPLACE,  "DiffOpP2",     true  } },
  {  "P2P1Stokes", { P2P1STOKES, "StokesOpP2P1", true  } }
};


void showUsage() {

  std::stringstream mesg;

  mesg << "Please specify the following two parameters in the given order:\n\n"
       << "  <level>      on which level do you want the operator to be set up?\n"
       << "  <operator>   which operator do you want to export?\n\n"
       << " Choices available for <operator> are\n";

  for( auto it = oprMap.begin(); it != oprMap.end(); ++it ) {
    mesg << " - " << it->first << "\n";
  }

  WALBERLA_LOG_INFO_ON_ROOT( "" << mesg.str() );
}


int main( int argc, char* argv[] ) {

  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

  WALBERLA_LOG_DEVEL_ON_ROOT( "" );
  WALBERLA_LOG_DEVEL_ON_ROOT( "----------------------------------------------" );
  WALBERLA_LOG_DEVEL_ON_ROOT( "--------- Matrix Export of Operators ---------" );
  WALBERLA_LOG_DEVEL_ON_ROOT( "----------------------------------------------" );
  WALBERLA_LOG_DEVEL_ON_ROOT( "" );

  // Process command-line
  if( argc < 3 || argc > 3 ) {
    showUsage();
    WALBERLA_ABORT( "\n" );
  }

  uint_t level = static_cast<uint_t>( std::stoul( argv[1] ) );
  std::string oprName  = ( argv[2] );
  operatorTag oprTag   = oprMap[ oprName ].oprEnum;
  std::string matName  = oprMap[ oprName ].matName;
  std::string fileName = oprName + ".m";
  bool elim = oprMap[ oprName ].elimDirichletBC;

  // Mesh generation
  WALBERLA_LOG_INFO_ON_ROOT( "Generating criss mesh on unit square" );
  MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
  Point2D cornerLL( { 0.0, 0.0 } );
  Point2D cornerUR( { 1.0, 1.0 } );
  meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CRISS, 1, 1 );

  SetupPrimitiveStorage setupStorage( meshInfo,
                                      uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hyteg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // Operator creation and export
  switch( oprTag ) {

    // --------------
    //  P1 operators
    // --------------
  case P1LAPLACE:
    {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Laplace operator for P1 elements" );
      hyteg::P1ConstantLaplaceOperator opr( storage, level, level );
      exportOperator< P1ConstantLaplaceOperator, P1Function, P1FunctionTag >( opr, fileName, matName, storage, level, elim );
    }
    break;

  case P1MASS:
    {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Mass operator for P1 elements" );
      hyteg::P1ConstantMassOperator opr( storage, level, level );
      exportOperator< P1ConstantMassOperator, P1Function, P1FunctionTag >( opr, fileName, matName, storage, level, elim );
    }
    break;

    // --------------
    //  P2 operators
    // --------------
  case P2LAPLACE:
    {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Laplace operator for P2 elements" );
      hyteg::P2ConstantLaplaceOperator opr( storage, level, level );
      exportOperator< P2ConstantLaplaceOperator, P2Function, P2FunctionTag >( opr, fileName, matName, storage, level, elim );
    }
    break;

  case P2MASS:
    {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Mass operator for P2 elements" );
      hyteg::P2ConstantMassOperator opr( storage, level, level );
      exportOperator< P2ConstantMassOperator, P2Function, P2FunctionTag >( opr, fileName, matName, storage, level, elim );
    }
    break;

    // -----------------
    //  Mixed operators
    // -----------------
  case P2P1STOKES:
    {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Stokes Operator for P2-P1 element" );
      hyteg::P2P1TaylorHoodStokesOperator opr( storage, level, level );
      exportOperator< P2P1TaylorHoodStokesOperator, P2P1TaylorHoodFunction, P2P1TaylorHoodFunctionTag >( opr, fileName, matName, storage, level, elim );
    }
    break;
  }

  WALBERLA_LOG_DEVEL_ON_ROOT( "----------------------------------------------" );

  return 0;
}
