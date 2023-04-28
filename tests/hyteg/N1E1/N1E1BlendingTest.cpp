/*
 * Copyright (c) 2023 Daniel Bauer.
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

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using namespace hyteg;

std::shared_ptr< PrimitiveStorage > storageWithAffineBlending()
{
   const auto            meshInfo = MeshInfo::meshSymmetricCuboid( { 0, 0, 0 }, { 1, 1, 1 }, 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   Matrix3r B;
   // clang-format off
   B <<  0.1, 0.2, 0.3,
        -0.8, 0.5, 0.0,
         1.0, 1.0, 0.5;
   // clang-format on
   const Point3D b{ 2.0, 3.0, 4.0 };
   AffineMap3D::setMap( setupStorage, B, b );

   return std::make_shared< PrimitiveStorage >( setupStorage );
}

std::shared_ptr< PrimitiveStorage > storageWithNonAffineBlending()
{
   // the resoultion is as low as possible to make the error that would result
   // from incorrect blending as large as possible
   const auto            meshInfo = MeshInfo::meshSphericalShell( 2, 2, 0.5, 1.0 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   IcosahedralShellMap::setMap( setupStorage );

   return std::make_shared< PrimitiveStorage >( setupStorage );
}

template < typename F >
std::pair< real_t, real_t > test( std::shared_ptr< PrimitiveStorage > storage,
                                  const F                             exact,
                                  const uint_t                        level,
                                  const bool                          check    = true,
                                  const bool                          writeVTK = false )
{
   using namespace n1e1;

   const uint_t minLevel                    = level;
   const uint_t maxLevel                    = level;
   const uint_t numRandomEvaluationsPerCell = 10;
   walberla::math::seedRandomGenerator( 42 );

   if ( writeVTK )
   {
      writeDomainPartitioningVTK( storage, "../../output", "N1E1BlendingTest_partitioning" );
   }

   // interpolate
   N1E1VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
   u.interpolate( exact, maxLevel, DoFType::All );

   // by communicating, we test the implementation on faces and edges, too
   u.communicate< Face, Cell >( maxLevel );
   u.communicate< Edge, Cell >( maxLevel );

   real_t sum = 0.0;
   real_t max = 0.0;
   int    n   = 0;

   for ( const auto& it : storage->getCells() )
   {
      Cell& cell = *it.second;

      for ( uint_t i = 0; i < numRandomEvaluationsPerCell; ++i )
      {
         const std::array< Point3D, 4 >& verts = cell.getCoordinates();

         const real_t r0 = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         const real_t r1 = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         const real_t r2 = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         const real_t r3 = real_c( walberla::math::realRandom( 0.0, 1.0 ) );

         const Point3D compCoords = ( r0 * verts[0] + r1 * verts[1] + r2 * verts[2] + r3 * verts[3] ) / ( r0 + r1 + r2 + r3 );

         Point3D physCoords;
         cell.getGeometryMap()->evalF( compCoords, physCoords );

         Point3D eval;
         auto    success = u.evaluate( physCoords, level, eval );
         WALBERLA_CHECK( success );

         Point3D exactAtPhys;
         if constexpr ( std::is_same< F, Point3D >::value )
         {
            exactAtPhys = exact;
         }
         else
         {
            exactAtPhys = exact( physCoords );
         }

         const real_t err = ( eval - exactAtPhys ).norm();
         sum += err;
         max = std::max( max, err );
         n += 1;

         if ( check )
         {
            WALBERLA_CHECK_FLOAT_EQUAL( eval[0], exactAtPhys[0], "Test3D: wrong X-coordinate at " << physCoords << "." );
            WALBERLA_CHECK_FLOAT_EQUAL( eval[1], exactAtPhys[1], "Test3D: wrong Y-coordinate at " << physCoords << "." );
            WALBERLA_CHECK_FLOAT_EQUAL( eval[2], exactAtPhys[2], "Test3D: wrong Z-coordinate at " << physCoords << "." );
         }
      }
   }

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "N1E1BlendingTest", storage );
      vtk.add( u );
      vtk.add( *u.getDoFs() );
      vtk.write( level );
   }

   return { sum / real_c( n ), max };
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   // analytical function in N1E1 FE space
   const Point3D                                    a0{ 1, 2, 3 };
   const Point3D                                    a1{ 4, 5, 6 };
   const std::function< Point3D( const Point3D& ) > exact = [&]( const Point3D& x ) { return ( a0 + a1.cross( x ) ).eval(); };

   WALBERLA_LOG_INFO_ON_ROOT( "Affine non-constant" )
   test( storageWithAffineBlending(), exact, 3 );

   WALBERLA_LOG_INFO_ON_ROOT( "Affine constant" )
   test( storageWithAffineBlending(), a0, 3 );

   // With non-affine blending, functions can not be represented exactly in computational space.
   // Therefore, we test linear convergence (by some margin).
   const real_t convFactor = 0.51;

   {
      WALBERLA_LOG_INFO_ON_ROOT( "Non-affine non-constant" )
      auto [avgCoarse, maxCoarse] = test( storageWithNonAffineBlending(), a0, 3, false );
      auto [avgFine, maxFine]     = test( storageWithNonAffineBlending(), a0, 4, false );

      WALBERLA_CHECK_LESS( avgFine, 0.1 )
      WALBERLA_CHECK_LESS( avgFine, convFactor * avgCoarse )
      WALBERLA_CHECK_LESS( maxFine, convFactor * maxCoarse )
   }

   {
      WALBERLA_LOG_INFO_ON_ROOT( "Non-affine constant" )
      auto [avgCoarse, maxCoarse] = test( storageWithNonAffineBlending(), a0, 3, false );
      auto [avgFine, maxFine]     = test( storageWithNonAffineBlending(), a0, 4, false );

      WALBERLA_CHECK_LESS( avgFine, 0.1 )
      WALBERLA_CHECK_LESS( avgFine, convFactor * avgCoarse )
      WALBERLA_CHECK_LESS( maxFine, convFactor * maxCoarse )
   }
   return EXIT_SUCCESS;
}
