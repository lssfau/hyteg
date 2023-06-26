/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/MeshQuality.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;


real_t errorE1( const uint_t&                   level,
                const P2Function< real_t >&     c,
                const P2Function< real_t >&     solution,
                const P2Function< real_t >&     tmp0,
                const P2Function< real_t >&     tmp1,
                const P2ConstantRowSumOperator& lumpedMass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E1 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::abs( values[0] - values[1] ); };

   tmp0.interpolate( E1, {solution, c}, level, All );
   lumpedMass.apply( tmp0, tmp1, level, All );
   return tmp1.sumGlobal( level, All );
}

real_t errorE2( const uint_t&                   level,
                const P2Function< real_t >&     c,
                const P2Function< real_t >&     solution,
                const P2Function< real_t >&     tmp0,
                const P2Function< real_t >&     tmp1,
                const P2ConstantRowSumOperator& lumpedMass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E2 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::pow( std::abs( values[0] - values[1] ), 2 ); };

   tmp0.interpolate( E2, {solution, c}, level, All );
   lumpedMass.apply( tmp0, tmp1, level, All );
   return std::sqrt( tmp1.sumGlobal( level, All ) );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo meshInfo = hyteg::MeshInfo::meshCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   storage->getTimingTree()->start( "Total" );

   const uint_t minLevel   = 2;
   const uint_t maxLevel   = 4;

   const real_t hMin       = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t hMax       = MeshQuality::getMaximalEdgeLength( storage, maxLevel );
   const real_t tEnd       = 2 * walberla::math::pi;
   const real_t dt         = tEnd / 60.;
   const uint_t stepsTotal = 2; // uint_c( tEnd / dt );

   const uint_t outerSteps = 2; // 60;
   const uint_t innerSteps = stepsTotal / outerSteps;

   const bool vtk = false;

   const auto unknowns = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( unknowns )

   WALBERLA_LOG_INFO_ON_ROOT( "Circular convection" )
   WALBERLA_LOG_INFO_ON_ROOT( " - dt:                                           " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - h (min, max):                                 " << hMin << ", " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - level:                                        " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( " - inner time steps (no temperature evaluation): " << innerSteps )
   WALBERLA_LOG_INFO_ON_ROOT( " - outer time steps:                             " << outerSteps )
   WALBERLA_LOG_INFO_ON_ROOT( " - steps until circle completed:                 " << stepsTotal )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   auto r = []( const hyteg::Point3D& x, const hyteg::Point3D& x0, const real_t& r0 ) -> real_t {
     return ( 1 / r0 ) * std::sqrt( std::pow( x[0] - x0[0], 2 ) + std::pow( x[1] - x0[1], 2 ) + std::pow( x[2] - x0[2], 2 ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > conicalBody = [&]( const hyteg::Point3D& x ) -> real_t {
     const Point3D x0( 0.5, 0.25, 0.5 );
     const real_t  r0 = 0.15;
     if ( r( x, x0, r0 ) <= 1. )
        return 1 - r( x, x0, r0 );
     else
        return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > gaussianCone = [&]( const hyteg::Point3D& x ) -> real_t {
     const Point3D x0( 0.25, 0.5, 0.5 );
     const real_t  r0 = 0.15;
     if ( r( x, x0, r0 ) <= 1. )
        return ( 1 + std::cos( walberla::math::pi * r( x, x0, r0 ) ) ) * 0.25;
     else
        return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > slottedCylinder = [&]( const hyteg::Point3D& x ) -> real_t {
     const Point3D x0( 0.5, 0.75, 0.5 );
     const real_t  r0 = 0.15;
     if ( ( r( x, x0, r0 ) <= 1. ) && ( std::abs( x[0] - x0[0] ) >= 0.025 || x[1] >= 0.85 ) )
        return 1;
     else
        return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > initialBodies = [&]( const hyteg::Point3D& x ) -> real_t {
     return conicalBody( x ) + gaussianCone( x ) + slottedCylinder( x );
   };

   auto vel_x = []( const hyteg::Point3D& x ) -> real_t {
        return 0.5 - x[1];
   };

   auto vel_y = []( const hyteg::Point3D& x ) -> real_t {
        return x[0] - 0.5;
   };

   auto zero = []( const hyteg::Point3D& ) -> real_t {
        return real_t(0);
   };

   writeDomainPartitioningVTK( storage, "../../output", "MMOC3DCircularConvectionTest_Domain" );

   typedef P2Function< real_t >   FunctionType;
   typedef P2ConstantMassOperator MassOperator;

   FunctionType c( "c", storage, minLevel, maxLevel );
   FunctionType cInitial( "cInitial", storage, minLevel, maxLevel );
   FunctionType cError( "cError", storage, minLevel, maxLevel );
   FunctionType cMass( "cError", storage, minLevel, maxLevel );
   FunctionType velocityMagnitude( "velocityMagnitude", storage, minLevel, maxLevel );
   typename FunctionTrait< FunctionType >::AssocVectorFunctionType uvw( "uvw", storage, minLevel, maxLevel );

   FunctionType tmp0( "tmp0", storage, minLevel, maxLevel );
   FunctionType tmp1( "tmp1", storage, minLevel, maxLevel );

   MassOperator                  M( storage, minLevel, maxLevel );
   MMOCTransport< FunctionType > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   uvw.interpolate( { vel_x, vel_y, zero }, maxLevel );
   c.interpolate( initialBodies, maxLevel );
   cInitial.interpolate( initialBodies, maxLevel );

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > magnitude = []( const Point3D&,
                                                                                           const std::vector< real_t >& values ) {
     return Point3D( values[0], values[1], values[2] ).norm();
   };

   velocityMagnitude.interpolate( magnitude, {uvw[0], uvw[1], uvw[2]}, maxLevel, All );

   hyteg::VTKOutput vtkOutput( "../../output", "MMOC3DCircularConvectionTest", storage );

   vtkOutput.add( uvw );
   vtkOutput.add( c );
   vtkOutput.add( cError );
   vtkOutput.add( cInitial );

   WALBERLA_LOG_INFO_ON_ROOT( " outer step | timestep | max temperature | total mass | mass lost since last outer step " )
   WALBERLA_LOG_INFO_ON_ROOT( "------------+----------+-----------------+------------+---------------------------------" )

   cError.assign( {1.0, -1.0}, {c, cInitial}, maxLevel, All );
   auto max_temp = c.getMaxMagnitude( maxLevel, All );
   M.apply( c, cMass, maxLevel, All );
   auto total_mass = cMass.sumGlobal( maxLevel );

   if ( vtk )
      vtkOutput.write( maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %10d | %8d | %15.3e | %10.3e | %30.2f%% ", 0, 0, max_temp, total_mass, 0, 0. ) )

   for ( uint_t i = 1; i <= outerSteps; i++ )
   {
      transport.step( c, uvw, uvw, maxLevel, Inner, dt, innerSteps, i == 1 );

      cError.assign( {1.0, -1.0}, {c, cInitial}, maxLevel, All );
      max_temp = c.getMaxMagnitude( maxLevel, All );
      M.apply( c, cMass, maxLevel, All );
      auto total_mass_new  = cMass.sumGlobal( maxLevel );
      auto total_mass_lost = 1.0 - ( total_mass_new / total_mass );
      total_mass           = total_mass_new;

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          " %10d | %8d | %15.3e | %10.3e | %30.2f%% ", i, i * innerSteps, max_temp, total_mass, total_mass_lost * 100. ) )

      if ( vtk )
         vtkOutput.write( maxLevel, i * innerSteps );
   }

   auto p2MassFormFenics =
       std::make_shared< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >();
   P2RowSumForm                       rowSumMass( p2MassFormFenics );
   P2ConstantOperator< P2RowSumForm > MLumped( storage, minLevel, maxLevel, rowSumMass );

   auto error_E1 = errorE1( maxLevel, c, cInitial, tmp0, tmp1, MLumped );
   auto error_E2 = errorE2( maxLevel, c, cInitial, tmp0, tmp1, MLumped );

   WALBERLA_LOG_INFO_ON_ROOT( "E1: " << walberla::format( "%5.4e", error_E1 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "E2: " << walberla::format( "%5.4e", error_E2 ) );

   // WALBERLA_CHECK_LESS( error_E1, 1.2e-03 );
   // WALBERLA_CHECK_LESS( error_E2, 1.2e-02 );

   storage->getTimingTree()->stop( "Total" );
   WALBERLA_LOG_INFO_ON_ROOT( storage->getTimingTree()->getCopyWithRemainder() );

   writeTimingTreeJSON( *storage->getTimingTree(), "/tmp/mmoc_timing.json" );

   return EXIT_SUCCESS;
}
