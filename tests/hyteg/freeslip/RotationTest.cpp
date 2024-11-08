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

template < typename StokesFunctionType, typename RotationOperatorType, bool use3D >
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

   real_t ur     = 1.23;
   real_t uTheta = 1.1;
   real_t uPhi   = 0.85;

   std::function< std::array< real_t, 3u >( const Point3D& ) > getRThetaPhi = [&]( const Point3D& x ) {
      real_t r   = x.norm();
      real_t phi = std::atan2( x[1], x[0] );
      if ( storage->hasGlobalCells() )
      {
         real_t theta = std::acos( x[2] / r );
         return std::array< real_t, 3u >( { r, theta, phi } );
      }
      else
      {
         return std::array< real_t, 3u >( { r, phi, 0.0 } );
      }
   };

   std::function< std::array< real_t, 5u >( const Point3D& ) > getRSinThetaCosThetaSinPhiCosPhi = [&]( const Point3D& x ) {
      std::array< real_t, 3u > rThetaPhi = getRThetaPhi( x );
      if ( storage->hasGlobalCells() )
      {
         real_t sinTheta = std::sin( rThetaPhi[1] );
         real_t cosTheta = std::cos( rThetaPhi[1] );
         real_t sinPhi   = std::sin( rThetaPhi[2] );
         real_t cosPhi   = std::cos( rThetaPhi[2] );

         return std::array< real_t, 5u >( { rThetaPhi[0], sinTheta, cosTheta, sinPhi, cosPhi } );
      }
      else
      {
         real_t sinPhi = std::sin( rThetaPhi[1] );
         real_t cosPhi = std::cos( rThetaPhi[1] );

         return std::array< real_t, 5u >( { rThetaPhi[0], sinPhi, cosPhi, 0.0, 0.0 } );
      }
   };

   std::function< real_t( const Point3D& ) > uX = [&]( const Point3D& x ) {
      std::array< real_t, 5u > vals = getRSinThetaCosThetaSinPhiCosPhi( x );

      // real_t r = vals[0];

      if ( storage->hasGlobalCells() )
      {
         real_t sTheta = vals[1];
         real_t cTheta = vals[2];

         real_t sPhi = vals[3];
         real_t cPhi = vals[4];

         return sTheta * cPhi * ur + cTheta * cPhi * uTheta - sPhi * uPhi;
      }
      else
      {
         real_t sPhi = vals[1];
         real_t cPhi = vals[2];

         return ur * cPhi - uPhi * sPhi;
      }
   };

   std::function< real_t( const Point3D& ) > uY = [&]( const Point3D& x ) {
      std::array< real_t, 5u > vals = getRSinThetaCosThetaSinPhiCosPhi( x );

      // real_t r = vals[0];

      if ( storage->hasGlobalCells() )
      {
         real_t sTheta = vals[1];
         real_t cTheta = vals[2];

         real_t sPhi = vals[3];
         real_t cPhi = vals[4];

         return sTheta * sPhi * ur + cTheta * sPhi * uTheta + cPhi * uPhi;
      }
      else
      {
         real_t sPhi = vals[1];
         real_t cPhi = vals[2];

         return ur * sPhi + uPhi * cPhi;
      }
   };

   std::function< real_t( const Point3D& ) > uZ = [&]( const Point3D& x ) {
      std::array< real_t, 5u > vals = getRSinThetaCosThetaSinPhiCosPhi( x );

      // real_t r = vals[0];

      if ( storage->hasGlobalCells() )
      {
         real_t sTheta = vals[1];
         real_t cTheta = vals[2];

         // real_t sPhi = vals[3];
         // real_t cPhi = vals[4];

         return cTheta * ur - sTheta * uTheta;
      }
      else
      {
         WALBERLA_ABORT( "Shouldn't be here!" );
      }
   };

   RotationOperatorType rotationOperator( storage, level, level, normalFunction );

   BoundaryCondition bcVelocity;

   bcVelocity.createAllInnerBC();
   bcVelocity.createFreeslipBC( "FSInnerAndOuter",
                                { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );

   StokesFunctionType uRef( "uRef", storage, level, level, bcVelocity );
   StokesFunctionType u( "u", storage, level, level, bcVelocity );
   StokesFunctionType error( "error", storage, level, level, bcVelocity );

   if ( storage->hasGlobalCells() )
   {
      u.uvw().interpolate( { uX, uY, uZ }, level );
      uRef.uvw().interpolate( { uX, uY, uZ }, level );
   }
   else
   {
      u.uvw().interpolate( { uX, uY }, level );
      uRef.uvw().interpolate( { uX, uY }, level );
   }

   // VTKOutput vtkOutput( "../../output", "RotationTest", storage );

   // vtkOutput.add( u );
   // vtkOutput.add( uRef );
   // vtkOutput.write( level, 0u );

   // Apply rotation
   rotationOperator.rotate( u, level, FreeslipBoundary );
   // vtkOutput.write(level, 1u);

   auto dp = std::is_same< real_t, double >();

   real_t dpTol = dp ? 1e-14 : 2e-6;

   // we check if the rotation makes the ndim component ~ equal to the radial component
   uint_t radialComponentIdx = use3D ? 2u : 1u;

   real_t urMax = u.uvw().component( radialComponentIdx ).getMaxDoFValue( level, FreeslipBoundary );
   // real_t urMin = u.uvw().component( radialComponentIdx ).getMinDoFValue( level, FreeslipBoundary );

   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "urVal = %4.7e, urMax = %4.7e, urMin = %4.7e", ur, urMax, urMin ) );

   // WALBERLA_CHECK_LESS( std::abs( urMax - urMin ), dpTol );
   WALBERLA_CHECK_LESS( std::abs( urMax - ur ), dpTol );

   // Apply rotation transpose
   rotationOperator.rotate( u, level, FreeslipBoundary, true );
   // vtkOutput.write(level, 2u);

   // Check if we get back the original values after applying R^T R u
   error.assign( { 1.0, -1.0 }, { u, uRef }, level, All );
   real_t errorNorm = error.dotGlobal( error, level );

   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "errorNorm = %4.7e", errorNorm ) );

   WALBERLA_CHECK_LESS( errorNorm, dpTol );
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
