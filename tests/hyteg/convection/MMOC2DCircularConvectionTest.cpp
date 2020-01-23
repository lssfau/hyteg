/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <core/timing/Timer.h>

#include "hyteg/composites/P1MMOCTransport.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
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

   MeshInfo meshInfo = hyteg::MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, 3, 3 );
   // MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   const uint_t minLevel   = 2;
   const uint_t maxLevel   = 5;
   const uint_t outerSteps = 100;
   const uint_t innerSteps = 100;
   real_t       dt         = 1.0 * std::pow( 2.0, -walberla::real_c( maxLevel + 1 ) );
   WALBERLA_LOG_DEVEL( "dt = " << dt )

   auto r = []( const hyteg::Point3D& x, const hyteg::Point3D& x0, const real_t& r0 ) -> real_t {
      return ( 1 / r0 ) * std::sqrt( std::pow( x[0] - x0[0], 2 ) + std::pow( x[1] - x0[1], 2 ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > conicalBody = [&]( const hyteg::Point3D& x ) -> real_t {
      const Point3D x0( {0.5, 0.25, 0.0} );
      const real_t  r0 = 0.15;
      if ( r( x, x0, r0 ) <= 1. )
         return 1 - r( x, Point3D( {0.5, 0.25, 0.0} ), 0.15 );
      else
         return 0.0;
   };

   auto vel_x = []( const hyteg::Point3D& x ) -> real_t {
      if ( (x - Point3D({0.5, 0.5, 0})).norm() < 0.45 )
         return 0.5 - x[1];
      else
         return 0;
   };

   auto vel_y = []( const hyteg::Point3D& x ) -> real_t {
      if ( (x - Point3D({0.5, 0.5, 0})).norm() < 0.45 )
         return x[0] - 0.5;
      else
         return 0;
   };

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "MMOC2DCircularConvectionTest_Domain" );

   P1Function< real_t > c( "c", storage, minLevel, maxLevel );
   P1Function< real_t > cInitial( "cInitial", storage, minLevel, maxLevel );
   P1Function< real_t > cError( "cError", storage, minLevel, maxLevel );
   P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   P1Function< real_t > v( "v", storage, minLevel, maxLevel );
   P1Function< real_t > w( "w", storage, minLevel, maxLevel );

   hyteg::P1MMOCTransport transport( storage, minLevel, maxLevel, TimeSteppingScheme::ExplicitEuler, false );

   u.interpolate( vel_x, maxLevel );
   v.interpolate( vel_y, maxLevel );
   c.interpolate( conicalBody, maxLevel );
   cInitial.interpolate( conicalBody, maxLevel );

   hyteg::VTKOutput vtkOutput( "../../output", "MMOC2DCircularConvectionTest", storage, 100 );

   vtkOutput.add( u );
   vtkOutput.add( v );
   vtkOutput.add( c );
   vtkOutput.add( cInitial );

   cError.assign( {1.0, -1.0}, {c, cInitial}, maxLevel, All );
   auto max_error = cError.getMaxMagnitude( maxLevel, All );
   auto max_temp  = c.getMaxMagnitude( maxLevel, All );
   auto l2_temp =
       std::sqrt( c.dotGlobal( c, maxLevel, All ) / real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );

   vtkOutput.write( maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Timestep: " << 0 << ", max error magnitude = " << max_error << ", max temp = " << max_temp
                                           << ", l2-discr temp: " << l2_temp );

   for ( uint_t i = 1; i <= outerSteps; i++ )
   {
      transport.step( c, u, v, w, maxLevel, Inner, dt, innerSteps );

      cError.assign( {1.0, -1.0}, {c, cInitial}, maxLevel, All );
      max_error = cError.getMaxMagnitude( maxLevel, All );
      max_temp  = c.getMaxMagnitude( maxLevel, All );
      l2_temp =
          std::sqrt( c.dotGlobal( c, maxLevel, All ) / real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );
      WALBERLA_LOG_INFO_ON_ROOT( "Timestep: " << i * innerSteps << ", max error magnitude = " << max_error
                                              << ", max temp = " << max_temp << ", l2-discr temp: " << l2_temp );

      vtkOutput.write( maxLevel, i * innerSteps );
   }

   return EXIT_SUCCESS;
}
