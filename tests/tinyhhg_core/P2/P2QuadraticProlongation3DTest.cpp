
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/FunctionIterator.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

using hhg::indexing::Index;
using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

void testWeightsInCellVertexDoF()
{
   typedef edgedof::EdgeDoFOrientation eo;

   const uint_t lowerLevel = 2;

   const auto            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P2Function< real_t > u( "u", storage, lowerLevel, lowerLevel + 1 );

   for ( auto it : FunctionIterator< P1Function< real_t > >( u.getVertexDoFFunction(), lowerLevel ) )
   {
      if ( it.isOnMacroCell() && it.isVertexDoF() && it.index() == Index( {1, 1, 1} ) )
      {
         it.value() = 1.0;
      }
   }

   P2toP2QuadraticProlongation prolongationOperator;
   prolongationOperator.prolongate( u, lowerLevel, Inner | NeumannBoundary );

   std::map< eo, uint_t > numNeighborElements = {
       {eo::X, 6},
       {eo::Y, 4},
       {eo::Z, 6},
       {eo::XY, 6},
       {eo::XZ, 4},
       {eo::YZ, 6},
       {eo::XYZ, 4},
   };

   uint_t numModifiedDoFs = 0;
   for ( auto it : FunctionIterator< P1Function< real_t > >( u.getVertexDoFFunction(), lowerLevel + 1 ) )
   {
      if ( it.isOnMacroCell() && it.isVertexDoF() && std::abs( it.value() ) > 1e-8 )
      {
         numModifiedDoFs++;
         WALBERLA_LOG_INFO_ON_ROOT( it );
      }
   }

   for ( auto it : FunctionIterator< EdgeDoFFunction< real_t > >( u.getEdgeDoFFunction(), lowerLevel + 1 ) )
   {
      if ( it.isOnMacroCell() && it.isEdgeDoF() && std::abs( it.value() ) > 1e-8 )
      {
         numModifiedDoFs++;
         WALBERLA_LOG_INFO_ON_ROOT( it );
      }
   }

   WALBERLA_CHECK_EQUAL( numModifiedDoFs, 1 + 2 * 14 + 3 * 24 + 24 )
}


void testWeightsInCellEdgeDoF()
{
  typedef edgedof::EdgeDoFOrientation eo;

  const uint_t lowerLevel = 2;

  const auto            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

  P2Function< real_t > u( "u", storage, lowerLevel, lowerLevel + 1 );

  for ( auto it : FunctionIterator< EdgeDoFFunction< real_t > >( u.getEdgeDoFFunction(), lowerLevel ) )
  {
    if ( it.isOnMacroCell() && it.isEdgeDoF() && it.edgeDoFOrientation() == eo::Z && it.index() == Index( {1, 1, 0} ) )
    {
      WALBERLA_LOG_INFO_ON_ROOT( it );
      it.value() = 1.0;
    }
  }

  P2toP2QuadraticProlongation prolongationOperator;
  prolongationOperator.prolongate( u, lowerLevel, Inner | NeumannBoundary );

  std::map< eo, uint_t > numNeighborElements = {
  {eo::X, 6},
  {eo::Y, 4},
  {eo::Z, 6},
  {eo::XY, 6},
  {eo::XZ, 4},
  {eo::YZ, 6},
  {eo::XYZ, 4},
  };

  uint_t numModifiedDoFs = 0;
  for ( auto it : FunctionIterator< P1Function< real_t > >( u.getVertexDoFFunction(), lowerLevel + 1 ) )
  {
    if ( it.isOnMacroCell() && it.isVertexDoF() && std::abs( it.value() ) > 1e-8 )
    {
      numModifiedDoFs++;
      WALBERLA_LOG_INFO_ON_ROOT( it );
    }
  }

  for ( auto it : FunctionIterator< EdgeDoFFunction< real_t > >( u.getEdgeDoFFunction(), lowerLevel + 1 ) )
  {
    if ( it.isOnMacroCell() && it.isEdgeDoF() && std::abs( it.value() ) > 1e-8 )
    {
      numModifiedDoFs++;
      WALBERLA_LOG_INFO_ON_ROOT( it );
    }
  }

  WALBERLA_CHECK_EQUAL( numModifiedDoFs, 1 + 2 + 6 * 2 + 6 * 1 + 6 * 1 )
}


void testGridTransfer3D( const std::string& meshFile, const uint_t& lowerLevel )
{
   const bool   writeVTK   = true;
   const real_t errorLimit = 1e-15;

   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // WALBERLA_CHECK( storage->hasGlobalCells() );

   if ( writeVTK )
      writeDomainPartitioningVTK( storage, "../../output", "P1LaplaceOperatorTest3D_partitioning" );

   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) -> real_t { return 0.0; };

   std::function< real_t( const hhg::Point3D& ) > one = []( const hhg::Point3D& ) -> real_t { return 1.0; };

   std::function< real_t( const hhg::Point3D& ) > constant = []( const hhg::Point3D& ) -> real_t { return 42.0; };

   std::function< real_t( const hhg::Point3D& ) > linearInX = []( const hhg::Point3D& p ) -> real_t {
      return real_c( 42 ) * p[0];
   };

   std::function< real_t( const hhg::Point3D& ) > linearInXYZ = []( const hhg::Point3D& p ) -> real_t {
      return real_c( 42 ) * p[0] + p[1] + real_c( 1337 ) * p[2];
   };

   std::function< real_t( const hhg::Point3D& ) > quadraticInXYZ = []( const hhg::Point3D& p ) -> real_t {
      return 2. * p[0] * p[0] + 3. * p[0] + 13. + 4. * p[1] + 5. * p[1] * p[1] + p[2] * p[2] + 6.;
   };

   P2Function< real_t > u( "u", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > resultExact( "u_exact", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > err( "err", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > oneFunction( "oneFunction", storage, lowerLevel, lowerLevel + 1 );

   VTKOutput vtkOutput( "../../output", "P2QuadraticProlongationTest3D", storage );
   vtkOutput.add( u );
   vtkOutput.add( resultExact );
   vtkOutput.add( err );

   auto testProlongationResult = [&]( std::function< real_t( const hhg::Point3D& ) > uFunction ) -> real_t {
      u.interpolate( uFunction, lowerLevel, All );
      resultExact.interpolate( uFunction, lowerLevel + 1, Inner | NeumannBoundary );

      P2toP2QuadraticProlongation prolongationOperator;
      prolongationOperator.prolongate( u, lowerLevel, Inner | NeumannBoundary );

      err.assign( {1.0, -1.0}, {u, resultExact}, lowerLevel + 1, Inner | NeumannBoundary );
      const real_t discrErr = err.dotGlobal( err, lowerLevel + 1, Inner | NeumannBoundary );
      return discrErr;
   };

   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 0 );

   // 1. u = const
   // ------------
   //   a) u = 0
   const real_t errorUZero = testProlongationResult( zero );
   WALBERLA_LOG_INFO_ON_ROOT( "u = 0: L2 error: " << errorUZero );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 1 );
   WALBERLA_CHECK_LESS( errorUZero, errorLimit );

   //   b) u = 1
   const real_t errorUOne = testProlongationResult( one );
   WALBERLA_LOG_INFO_ON_ROOT( "u = 1: L2 error: " << errorUOne );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 2 );
   WALBERLA_CHECK_LESS( errorUOne, errorLimit );

   //   c) u = some other constant
   const real_t errorUConstant = testProlongationResult( constant );
   WALBERLA_LOG_INFO_ON_ROOT( "u = const: L2 error: " << errorUConstant );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 3 );
   WALBERLA_CHECK_LESS( errorUConstant, errorLimit );

   // 2. u linear
   // -----------
   //   a) u linear in x
   const real_t errorULinearInX = testProlongationResult( linearInX );
   WALBERLA_LOG_INFO_ON_ROOT( "u linear in x: L2 error: " << errorULinearInX );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 4 );
   WALBERLA_CHECK_LESS( errorULinearInX, errorLimit );

   //   b) u linear in x, y and z
   const real_t errorULinearInXYZ = testProlongationResult( linearInXYZ );
   WALBERLA_LOG_INFO_ON_ROOT( "u linear in x, y and z: L2 error: " << errorULinearInXYZ );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 5 );
   WALBERLA_CHECK_LESS( errorULinearInXYZ, errorLimit );

   // 3. u quadratic
   // --------------
   //   a) u quadratic in x, y and z
   const real_t errorUQuadraticInXYZ = testProlongationResult( quadraticInXYZ );
   WALBERLA_LOG_INFO_ON_ROOT( "u quadratic in x, y and z: L2 error: " << errorUQuadraticInXYZ );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 6 );
   WALBERLA_CHECK_LESS( errorUQuadraticInXYZ, errorLimit );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   testWeightsInCellVertexDoF();
   testWeightsInCellEdgeDoF();

   testGridTransfer3D( "../../data/meshes/quad_8el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/tet_1el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_2el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_4el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 3 );

   return 0;
}
