/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionIterator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"

using hyteg::indexing::Index;
using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

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
         WALBERLA_LOG_DETAIL_ON_ROOT( it );
      }
   }

   for ( auto it : FunctionIterator< EdgeDoFFunction< real_t > >( u.getEdgeDoFFunction(), lowerLevel + 1 ) )
   {
      if ( it.isOnMacroCell() && it.isEdgeDoF() && std::abs( it.value() ) > 1e-8 )
      {
         numModifiedDoFs++;
         WALBERLA_LOG_DETAIL_ON_ROOT( it );
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
   const real_t errorLimit = real_c( std::is_same< real_t, double >() ? 1e-15 : 4e-5 );

   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // WALBERLA_CHECK( storage->hasGlobalCells() );

   if ( writeVTK )
      writeDomainPartitioningVTK( storage, "../../output", "P1LaplaceOperatorTest3D_partitioning" );

   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) -> real_t { return 0.0; };

   std::function< real_t( const hyteg::Point3D& ) > one = []( const hyteg::Point3D& ) -> real_t { return 1.0; };

   std::function< real_t( const hyteg::Point3D& ) > constant = []( const hyteg::Point3D& ) -> real_t { return 42.0; };

   std::function< real_t( const hyteg::Point3D& ) > linearInX = []( const hyteg::Point3D& p ) -> real_t {
      return real_c( 42 ) * p[0];
   };

   std::function< real_t( const hyteg::Point3D& ) > linearInXYZ = []( const hyteg::Point3D& p ) -> real_t {
      return real_c( 42 ) * p[0] + p[1] + real_c( 1337 ) * p[2];
   };

   std::function< real_t( const hyteg::Point3D& ) > quadraticInXYZ = []( const hyteg::Point3D& p ) -> real_t {
      return real_c( 2. * p[0] * p[0] + 3. * p[0] + 13. + 4. * p[1] + 5. * p[1] * p[1] + p[2] * p[2] + 6. );
   };

   P2Function< real_t > u( "u", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > resultExact( "u_exact", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > err( "err", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > oneFunction( "oneFunction", storage, lowerLevel, lowerLevel + 1 );

   VTKOutput vtkOutput( "../../output", "P2QuadraticProlongationTest3D", storage );
   vtkOutput.add( u );
   vtkOutput.add( resultExact );
   vtkOutput.add( err );

   auto testProlongationResult = [&]( std::function< real_t( const hyteg::Point3D& ) > uFunction ) -> real_t {
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
   WALBERLA_LOG_DETAIL_ON_ROOT( "u = 0: L2 error: " << errorUZero );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 1 );
   WALBERLA_CHECK_LESS( errorUZero, errorLimit );

   //   b) u = 1
   const real_t errorUOne = testProlongationResult( one );
   WALBERLA_LOG_DETAIL_ON_ROOT( "u = 1: L2 error: " << errorUOne );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 2 );
   WALBERLA_CHECK_LESS( errorUOne, errorLimit );

   //   c) u = some other constant
   const real_t errorUConstant = testProlongationResult( constant );
   WALBERLA_LOG_DETAIL_ON_ROOT( "u = const: L2 error: " << errorUConstant );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 3 );
   WALBERLA_CHECK_LESS( errorUConstant, errorLimit );

   // 2. u linear
   // -----------
   //   a) u linear in x
   const real_t errorULinearInX = testProlongationResult( linearInX );
   WALBERLA_LOG_DETAIL_ON_ROOT( "u linear in x: L2 error: " << errorULinearInX );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 4 );
   WALBERLA_CHECK_LESS( errorULinearInX, errorLimit );

   //   b) u linear in x, y and z
   const real_t errorULinearInXYZ = testProlongationResult( linearInXYZ );
   WALBERLA_LOG_DETAIL_ON_ROOT( "u linear in x, y and z: L2 error: " << errorULinearInXYZ );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 5 );
   WALBERLA_CHECK_LESS( errorULinearInXYZ, errorLimit );

   // 3. u quadratic
   // --------------
   //   a) u quadratic in x, y and z
   const real_t errorUQuadraticInXYZ = testProlongationResult( quadraticInXYZ );
   WALBERLA_LOG_DETAIL_ON_ROOT( "u quadratic in x, y and z: L2 error: " << errorUQuadraticInXYZ );
   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 6 );
   WALBERLA_CHECK_LESS( errorUQuadraticInXYZ, errorLimit );
}

void testProlongateAndAdd3D( const std::string& meshFile, const uint_t& lowerLevel )
{
   const bool   writeVTK   = true;
   const real_t errorLimit = real_c( std::is_same< real_t, double >() ? 1e-15 : 1e-9 );

   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   if ( writeVTK )
      writeDomainPartitioningVTK( storage, "../../output", "P1LaplaceOperatorTest3D_partitioning" );

   std::function< real_t( const hyteg::Point3D& ) > termProlongate = []( const hyteg::Point3D& x ) {
      return std::sin( x[0] ) + std::sinh( x[1] ) + x[2] * x[2];
   };
   std::function< real_t( const hyteg::Point3D& ) > termAdd = []( const hyteg::Point3D& x ) {
      return std::cos( x[0] ) + std::cosh( x[1] ) + x[2] * x[2] * x[2];
   };

   P2Function< real_t > uProlongate( "uProlongate", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > uAdd( "uAdd", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > uProlongateAndAdd( "uProlongateAndAdd", storage, lowerLevel, lowerLevel + 1 );
   P2Function< real_t > err( "err", storage, lowerLevel, lowerLevel + 1 );

   VTKOutput vtkOutput( "../../output", "P2QuadraticProlongationTest3D_testProlongateAndAdd", storage );
   vtkOutput.add( uProlongate );
   vtkOutput.add( uAdd );
   vtkOutput.add( uProlongateAndAdd );
   vtkOutput.add( err );

   uProlongate.interpolate( termProlongate, lowerLevel, All );
   uProlongate.interpolate( termAdd, lowerLevel + 1, All );

   uAdd.interpolate( termAdd, lowerLevel + 1, All );

   uProlongateAndAdd.interpolate( termProlongate, lowerLevel, All );
   uProlongateAndAdd.interpolate( termAdd, lowerLevel + 1, All );

   P2toP2QuadraticProlongation prolongationOperator;

   prolongationOperator.prolongate( uProlongate, lowerLevel, Inner | NeumannBoundary );
   prolongationOperator.prolongateAndAdd( uProlongateAndAdd, lowerLevel, Inner | NeumannBoundary );

   uProlongate.add( {1.0}, {uAdd}, lowerLevel + 1, Inner | NeumannBoundary );

   err.assign( {1.0, -1.0}, {uProlongate, uProlongateAndAdd}, lowerLevel + 1, Inner | NeumannBoundary );
   const real_t discrErr = err.dotGlobal( err, lowerLevel + 1, Inner | NeumannBoundary );

   if ( writeVTK )
      vtkOutput.write( lowerLevel + 1, 1 );

   WALBERLA_LOG_INFO_ON_ROOT( "L2 error: " << discrErr );
   WALBERLA_CHECK_LESS( discrErr,
                        errorLimit,
                        "Prolongate and add test failed. Level (fine): " + std::to_string( lowerLevel + 1 ) +
                            ", Mesh: " + meshFile );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   testWeightsInCellVertexDoF();
   testWeightsInCellEdgeDoF();

   testGridTransfer3D( "../../data/meshes/quad_8el.msh", 0 );
   testGridTransfer3D( "../../data/meshes/bfs_126el.msh", 0 );
   testGridTransfer3D( "../../data/meshes/3D/tet_1el.msh", 0 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_2el.msh", 0 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_4el.msh", 0 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 0 );
   testGridTransfer3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 0 );

   testGridTransfer3D( "../../data/meshes/quad_8el.msh", 1 );
   testGridTransfer3D( "../../data/meshes/3D/tet_1el.msh", 1 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_2el.msh", 1 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_4el.msh", 1 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 1 );
   testGridTransfer3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 1 );

   testGridTransfer3D( "../../data/meshes/quad_8el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/tet_1el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_2el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_4el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 3 );
   testGridTransfer3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 3 );

   testProlongateAndAdd3D( "../../data/meshes/quad_8el.msh", 0 );
   testProlongateAndAdd3D( "../../data/meshes/bfs_126el.msh", 0 );
   testProlongateAndAdd3D( "../../data/meshes/3D/tet_1el.msh", 0 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_2el.msh", 0 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_4el.msh", 0 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 0 );
   testProlongateAndAdd3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 0 );

   testProlongateAndAdd3D( "../../data/meshes/quad_8el.msh", 1 );
   testProlongateAndAdd3D( "../../data/meshes/3D/tet_1el.msh", 1 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_2el.msh", 1 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_4el.msh", 1 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 1 );
   testProlongateAndAdd3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 1 );

   testProlongateAndAdd3D( "../../data/meshes/quad_8el.msh", 3 );
   testProlongateAndAdd3D( "../../data/meshes/3D/tet_1el.msh", 3 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_2el.msh", 3 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_4el.msh", 3 );
   testProlongateAndAdd3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 3 );
   testProlongateAndAdd3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 3 );

   return 0;
}
