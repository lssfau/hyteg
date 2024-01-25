/*
 * Copyright (c) 2017-2023 Daniel Drzisga, Dominik Thoennes, Nils Kohl,
 * Marcus Mohr.
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

#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/composites/P0P1UpwindOperator.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

/// Setting from Kuzmin transport tutorial section 4.4.6.1

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo meshInfo = hyteg::MeshInfo::meshRectangle( Point2D(  0, 0  ), Point2D(  1, 1  ), MeshInfo::CRISS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   const uint_t minLevel = 2;
   const uint_t maxLevel = 8;
   // Comment in for "production" run
   // const uint_t timesteps = 6000;
   const uint_t timesteps = 40;
   real_t       dt        = 2 * walberla::math::pi / 6000.;
   WALBERLA_LOG_DEVEL( "dt = " << dt )

   auto r = []( const hyteg::Point3D& x, const hyteg::Point3D& x0, const real_t& r0 ) -> real_t {
      return ( 1 / r0 ) * std::sqrt( std::pow( x[0] - x0[0], 2 ) + std::pow( x[1] - x0[1], 2 ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > conicalBody = [&]( const hyteg::Point3D& x ) -> real_t {
      const Point3D x0(  0.5, 0.25, 0.0  );
      const real_t  r0 = 0.15;
      if ( r( x, x0, r0 ) <= 1. )
         return 1 - r( x, x0, r0 );
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > gaussianCone = [&]( const hyteg::Point3D& x ) -> real_t {
      const Point3D x0(  0.25, 0.5, 0.0  );
      const real_t  r0 = 0.15;
      if ( r( x, x0, r0 ) <= 1. )
         return ( 1 + std::cos( walberla::math::pi * r( x, x0, r0 ) ) ) * 0.25;
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > slottedCylinder = [&]( const hyteg::Point3D& x ) -> real_t {
      const Point3D x0(  0.5, 0.75, 0.0  );
      const real_t  r0 = 0.15;
      if ( ( r( x, x0, r0 ) <= 1. ) && ( std::abs( x[0] - x0[0] ) >= 0.025 || x[1] >= 0.85 ) )
         return 1;
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > initialBodies = [&]( const hyteg::Point3D& x ) -> real_t {
      return conicalBody( x ) + gaussianCone( x ) + slottedCylinder( x );
   };

   auto vel_x = []( const hyteg::Point3D& x ) -> real_t { return 0.5 - x[1]; };

   auto vel_y = []( const hyteg::Point3D& x ) -> real_t { return x[0] - 0.5; };

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1u );

   std::shared_ptr< hyteg::P0Function< real_t > > c_old =
       std::make_shared< hyteg::P0Function< real_t > >( "c_old", storage, minLevel, maxLevel );
   std::shared_ptr< hyteg::P0Function< real_t > > c =
       std::make_shared< hyteg::P0Function< real_t > >( "c", storage, minLevel, maxLevel );
   std::shared_ptr< hyteg::P0Function< real_t > > c_final =
       std::make_shared< hyteg::P0Function< real_t > >( "c_final", storage, minLevel, maxLevel );
   std::shared_ptr< hyteg::P0Function< real_t > > c_error =
       std::make_shared< hyteg::P0Function< real_t > >( "c_error", storage, minLevel, maxLevel );
   std::shared_ptr< hyteg::P1VectorFunction< real_t > > velocityPtr =
       std::make_shared< hyteg::P1VectorFunction< real_t > >( "velocity", storage, minLevel, maxLevel );

   // hyteg::DG0P1UpwindOperator transport( storage, *velocityPtr, minLevel, maxLevel );
   hyteg::P0P1UpwindOperator transport( storage, *velocityPtr, minLevel, maxLevel );

   velocityPtr->interpolate( { vel_x, vel_y }, maxLevel );
   c_old->interpolate( initialBodies, maxLevel );
   c->interpolate( initialBodies, maxLevel );
   c_final->interpolate( initialBodies, maxLevel );

   bool             storeOutput = false;
   hyteg::VTKOutput vtkOutput( "../../output", "Upwind2DCircularConvectionTest", storage, 10 );
   vtkOutput.add( *velocityPtr );
   vtkOutput.add( *c_old );
   vtkOutput.add( *c );
   vtkOutput.add( *c_final );
   if ( storeOutput )
   {
      vtkOutput.write( maxLevel );
   }

   real_t max_error;

   for ( uint_t i = 1; i <= timesteps; i++ )
   {
      if ( i % 10 == 0 )
      {
         c_error->assign( { 1.0, -1.0 }, { *c, *c_final }, maxLevel, All );
         // const auto max_error = c_error->getMaxMagnitude( maxLevel, All );
         max_error = c_error->getMaxMagnitude( maxLevel );
         WALBERLA_LOG_INFO_ON_ROOT( "Timestep: " << i << ", max error magnitude = " << max_error );
      }

      transport.apply( *c_old, *c, maxLevel, hyteg::Inner, Replace );
      c->assign( { 1.0, -dt }, { *c_old, *c }, maxLevel, hyteg::Inner );

      if ( storeOutput )
      {
         vtkOutput.write( maxLevel, i );
      }

      c_old.swap( c );
   }

   // Check values here have no physical meaning, I only inserted them
   // to check that algorithmic behaviour did not change [we test after
   // 40 iterations]
   WALBERLA_CHECK_LESS( max_error, real_c(1.000) )
   WALBERLA_CHECK_GREATER( max_error, real_c(0.999) )

   return EXIT_SUCCESS;
}
