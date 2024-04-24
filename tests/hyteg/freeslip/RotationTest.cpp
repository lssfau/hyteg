/*
 * Copyright (c) 2024 Ponsuganth Ilangovan
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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1RotationOperator.hpp"
#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_c;
using walberla::real_t;

using namespace hyteg;

/**
 * @brief This test mainly checks if the vector field is reconstructed back after applying R^T R
 *        Intermediately it also checks if the radial component of the rotated field is approximately correct
 */

template < typename StokesFunctionType, typename ProjectNormalOperatorType, bool use3D >
static void testRotation()
{
   const int level = 2;

   real_t rMin = 0.5;
   real_t rMax = 1.0;

   // real_t rMean = ( rMin + rMax ) / 2;

   auto meshInfo = MeshInfo::emptyMeshInfo();
   if ( use3D )
   {
      meshInfo = MeshInfo::meshSphericalShell( 2, 2, rMin, rMax );
   }
   else
   {
      meshInfo = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CRISS, 6, 6 );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( use3D )
   {
      IcosahedralShellMap::setMap( setupStorage );
   }
   else
   {
      AnnulusMap::setMap( setupStorage );
   }

   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto normalInterpolant = [=]( const Point3D& p ) {
      real_t norm = p.norm();
      return p / norm;
   };

   auto normalFunction = [=]( const Point3D& p, Point3D& n ) -> void { n = normalInterpolant( p ); };

   real_t ux = 1.23;
   real_t uy = 2.31;
   real_t uz = 3.12;

   std::function< real_t( const Point3D& ) > uX = [&]( const Point3D& ) { return ux; };
   std::function< real_t( const Point3D& ) > uY = [&]( const Point3D& ) { return uy; };
   std::function< real_t( const Point3D& ) > uZ = [&]( const Point3D& ) { return uz; };

   ProjectNormalOperatorType rotationOperator( storage, level, level, normalFunction );

   BoundaryCondition bcVelocity;

   bcVelocity.createAllInnerBC();
   bcVelocity.createFreeslipBC( "FSInnerAndOuter",
                                { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );

   StokesFunctionType u( "u", storage, level, level, bcVelocity );

   u.uvw()[0].interpolate( uX, level );
   u.uvw()[1].interpolate( uY, level );
   if ( storage->hasGlobalCells() )
      u.uvw()[2].interpolate( uZ, level );

   // Apply rotation
   rotationOperator.rotate( u, level, FreeslipBoundary );

   // we check if the rotation makes the ndim component ~ equal to the radial component
   if ( use3D )
   {
      real_t urTol = 1e-1;

      real_t urMax = u.uvw().component( 2 ).getMaxValue( level, FreeslipBoundary );

      real_t urVal = std::sqrt( ux * ux + uy * uy + uz * uz );

      // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "urVal = %4.7e, urMax = %4.7e", urVal, urMax ) );

      WALBERLA_CHECK_LESS( std::abs( urMax - urVal ), urTol );
   }
   else
   {
      real_t urTol = 1e-2;

      real_t urMax = u.uvw().component( 1 ).getMaxValue( level, FreeslipBoundary );

      real_t urVal = std::sqrt( ux * ux + uy * uy );

      // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "urVal = %4.7e, urMax = %4.7e", urVal, urMax ) );

      WALBERLA_CHECK_LESS( std::abs( urMax - urVal ), urTol );
   }

   // Apply rotation transpose
   rotationOperator.rotate( u, level, FreeslipBoundary, true );

   auto dp = std::is_same< real_t, double >();

   real_t uxMax = u.uvw().component( 0 ).getMaxValue( level, FreeslipBoundary );
   real_t uxMin = u.uvw().component( 0 ).getMinValue( level, FreeslipBoundary );

   real_t uyMax = u.uvw().component( 1 ).getMaxValue( level, FreeslipBoundary );
   real_t uyMin = u.uvw().component( 1 ).getMinValue( level, FreeslipBoundary );

   // WALBERLA_LOG_INFO_ON_ROOT(
   //     walberla::format( "uxMax = %4.7e, uxMin = %4.7e, uyMax = %4.7e, uyMin = %4.7e", uxMax, uxMin, uyMax, uyMin ) );

   real_t dpTol = dp ? 1e-14 : 2e-6;

   // Check if we get back the original values after applying R^T R u
   WALBERLA_CHECK_LESS( std::abs( uxMax - uxMin ), dpTol );
   WALBERLA_CHECK_LESS( std::abs( uyMax - uyMin ), dpTol );

   WALBERLA_CHECK_LESS( std::abs( uxMax - ux ), dpTol );
   WALBERLA_CHECK_LESS( std::abs( uyMax - uy ), dpTol );

   if ( use3D )
   {
      real_t uzMax = u.uvw().component( 2 ).getMaxValue( level, FreeslipBoundary );
      real_t uzMin = u.uvw().component( 2 ).getMinValue( level, FreeslipBoundary );

      // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "uzMax = %4.7e, uzMin = %4.7e", uzMax, uzMin ) );

      WALBERLA_CHECK_LESS( std::abs( uzMax - uzMin ), dpTol );
      WALBERLA_CHECK_LESS( std::abs( uzMax - uz ), dpTol );
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "normal projection P1-P1 in 2D" );
   testRotation< P1StokesFunction< real_t >, P1RotationOperator, false >();

   WALBERLA_LOG_INFO_ON_ROOT( "normal projection P2-P1-TH in 2D" );
   testRotation< P2P1TaylorHoodFunction< real_t >, P2RotationOperator, false >();

   WALBERLA_LOG_INFO_ON_ROOT( "normal projection P1-P1 in 3D" );
   testRotation< P1StokesFunction< real_t >, P1RotationOperator, true >();

   WALBERLA_LOG_INFO_ON_ROOT( "normal projection P2-P1-TH in 3D" );
   testRotation< P2P1TaylorHoodFunction< real_t >, P2RotationOperator, true >();
}
