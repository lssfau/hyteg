
#include <hyteg/mesh/MeshInfo.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/p1functionspace/P1Function.hpp>
#include <hyteg/p1functionspace/P1ConstantOperator.hpp>
#include <hyteg/p2functionspace/P2Function.hpp>
#include <hyteg/p2functionspace/P2ConstantOperator.hpp>
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "core/Environment.h"

namespace hyteg {

void TauExtrapolation3D( const MeshInfo & meshInfo, InnerEdgeType innerEdgeType, uint_t lowerLevel )
{
  // Writes the system matrices for an inner edge orientation type for the specified mesh.
  // Steps:
  //  1. Write the lower level P1 operator
  //     -> this is straightforward
  //  2. Write the lower level P2 operator
  //     -> this is also straightforward
  //  3. Write the higher level P1 operator
  //     -> here we need to be careful regarding the enumeration
  //        the lower level P1 entries must be enumerated exactly like the
  //        overlapping higher level P1 entries -> we therefore use an injection
  //        from the P2 enumeration at the lower level

  const uint_t higherLevel = lowerLevel + 1;

  std::array< uint_t, 4 > cellEnumerationDirection;
  std::string innerEdgeTypeString;
  switch ( innerEdgeType )
  {
    case InnerEdgeType::ALWAYS_1_n1_1:
      cellEnumerationDirection = {0, 1, 2, 3 };
      innerEdgeTypeString = "ALWAYS_1_n1_1";
      break;
    case InnerEdgeType::ALWAYS_n1_n1_1:
      cellEnumerationDirection = {1, 0, 2, 3 };
      innerEdgeTypeString = "ALWAYS_n1_n1_1";
      break;
    case InnerEdgeType::ALWAYS_n1_1_1:
      cellEnumerationDirection = {0, 2, 1, 3 };
      innerEdgeTypeString = "ALWAYS_n1_1_1";
      break;
    default:
      WALBERLA_ABORT("Invalid inner edge type." );
      break;
  }

  SetupPrimitiveStorage setupPrimitiveStorage( meshInfo, 1, innerEdgeType );
  setupPrimitiveStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  auto storage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage );

  P1Function< PetscInt > numeratorP1( "numeratorP1", storage, lowerLevel, higherLevel );
  P2Function< PetscInt > numeratorP2( "numeratorP2", storage, lowerLevel, lowerLevel );

  numeratorP1.enumerate( lowerLevel, cellEnumerationDirection );
  numeratorP2.enumerate( lowerLevel, cellEnumerationDirection );

  numeratorP1.assign( numeratorP2, higherLevel, All );
  communication::syncFunctionBetweenPrimitives( numeratorP1, higherLevel );

  P1ConstantLaplaceOperator laplaceP1( storage, lowerLevel, higherLevel );
  P2ConstantLaplaceOperator laplaceP2( storage, lowerLevel, lowerLevel );

  exportOperator< P1ConstantLaplaceOperator, P1Function, P1FunctionTag >( laplaceP1,
                                                                          numeratorP1,
                                                                          "laplace_P1_level_" + std::to_string( lowerLevel ) + "_" + innerEdgeTypeString +
                                                                              ".m",
                                                                          "P1_level_" + std::to_string( lowerLevel ) + "_" + innerEdgeTypeString,
                                                                          storage,
                                                                          lowerLevel,
                                                                          true,
                                                                          true );
  exportOperator< P1ConstantLaplaceOperator, P1Function, P1FunctionTag >( laplaceP1,
                                                                          numeratorP1,
                                                                          "laplace_P1_level_" + std::to_string( higherLevel )  + "_" + innerEdgeTypeString +
                                                                              ".m",
                                                                          "P1_level_" + std::to_string( higherLevel ) + "_" + innerEdgeTypeString,
                                                                          storage,
                                                                          higherLevel,
                                                                          true,
                                                                          true );
  exportOperator< P2ConstantLaplaceOperator, P2Function, P2FunctionTag >( laplaceP2,
                                                                          numeratorP2,
                                                                          "laplace_P2_level_" + std::to_string( lowerLevel )  + "_" + innerEdgeTypeString +
                                                                              ".m",
                                                                          "P2_level_" + std::to_string( lowerLevel ) + "_" + innerEdgeTypeString,
                                                                          storage,
                                                                          lowerLevel,
                                                                          true,
                                                                          true );
}

} // namespace hyteg

int main( int argc, char ** argv )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::MPIManager::instance()->useWorldComm();

  PETScManager petscManager( &argc, &argv );

  hyteg::TauExtrapolation3D( hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), hyteg::InnerEdgeType::ALWAYS_1_n1_1, walberla::uint_c(2) );
  hyteg::TauExtrapolation3D( hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), hyteg::InnerEdgeType::ALWAYS_n1_n1_1, walberla::uint_c(2) );
  hyteg::TauExtrapolation3D( hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), hyteg::InnerEdgeType::ALWAYS_n1_1_1, walberla::uint_c(2) );

  return EXIT_SUCCESS;
}