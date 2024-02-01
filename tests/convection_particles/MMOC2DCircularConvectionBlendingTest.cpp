/*
 * Copyright (c) 2017-2020 Nils Kohl.
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
#include <hyteg/geometry/AnnulusMap.hpp>

#include "hyteg/MeshQuality.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

/// Comparison of the HyTeG MMOC scheme with the solvers in Kuzmin 2010: A Guide to Numerical Methodsfor Transport Equations
/// Setting from transport tutorial section 4.4.6.1: Solid Body Rotation

/// Error definitions in section 4.4.6: Numerical examples

namespace hyteg {

template < typename FunctionType, typename MassOperator >
real_t errorE1( const uint_t&       level,
                const FunctionType& c,
                const FunctionType& solution,
                const FunctionType& tmp0,
                const FunctionType& tmp1,
                const MassOperator& mass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E1 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::abs( values[0] - values[1] ); };

   tmp0.interpolate( E1, {solution, c}, level, All );
   mass.apply( tmp0, tmp1, level, All );
   return tmp1.sumGlobal( level, All );
}

template < typename FunctionType, typename MassOperator >
real_t errorE2( const uint_t&       level,
                const FunctionType& c,
                const FunctionType& solution,
                const FunctionType& tmp0,
                const FunctionType& tmp1,
                const MassOperator& mass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E2 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::pow( std::abs( values[0] - values[1] ), 2 ); };

   tmp0.interpolate( E2, {solution, c}, level, All );
   mass.apply( tmp0, tmp1, level, All );
   return std::sqrt( tmp1.sumGlobal( level, All ) );
}

template < typename FunctionType, typename MassOperator >
void test( uint_t level, real_t errorE1Limit, real_t errorE2Limit )
{
   // MeshInfo meshInfo = MeshInfo::meshAnnulus( 0.5, 1.5, 0.0, 2.0 * walberla::math::pi, MeshInfo::CROSS, 6, 2 );
   MeshInfo meshInfo     = MeshInfo::meshAnnulus( 0.5, 1.5, MeshInfo::CROSS, 6, 2 );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   for ( auto it : setupStorage->getFaces() )
   {
      Face& face = *it.second;
      setupStorage->setGeometryMap( face.getID(), std::make_shared< AnnulusMap >( face ) );
   }
   for ( auto it : setupStorage->getEdges() )
   {
      Edge& edge = *it.second;
      setupStorage->setGeometryMap( edge.getID(), std::make_shared< AnnulusMap >( edge, *setupStorage ) );
   }

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   storage->getTimingTree()->start( "Total" );

   const uint_t minLevel   = 2;
   const uint_t maxLevel   = level;
   const real_t hMin       = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t hMax       = MeshQuality::getMaximalEdgeLength( storage, maxLevel );
   const real_t tEnd       = 2 * walberla::math::pi;
   const real_t dt         = tEnd / 150;
   const uint_t stepsTotal = uint_c( tEnd / dt );

   const bool vtk = false;

   const uint_t innerSteps = 10;
   const uint_t outerSteps = stepsTotal / innerSteps;

   WALBERLA_LOG_INFO_ON_ROOT( "Circular convection" )
   WALBERLA_LOG_INFO_ON_ROOT( " - dt:                                           " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - h (min, max):                                 " << hMin << ", " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - level:                                        " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( " - inner time steps (no temperature evaluation): " << innerSteps )
   WALBERLA_LOG_INFO_ON_ROOT( " - outer time steps:                             " << outerSteps )
   WALBERLA_LOG_INFO_ON_ROOT( " - steps until circle completed:                 " << stepsTotal )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   auto r = []( const hyteg::Point3D& x, const hyteg::Point3D& x0, const real_t& r0 ) -> real_t {
      return ( 1 / r0 ) * std::sqrt( std::pow( x[0] - x0[0], 2 ) + std::pow( x[1] - x0[1], 2 ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > conicalBody = [&]( const hyteg::Point3D& x ) -> real_t {
      const Point3D x0( 0, -0.75, 0.0 );
      const real_t  r0 = 0.15;
      if ( r( x, x0, r0 ) <= 1. )
         return 1 - r( x, x0, r0 );
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > gaussianCone = [&]( const hyteg::Point3D& x ) -> real_t {
      const Point3D x0( -0.75, 0.0, 0.0 );
      const real_t  r0 = 0.15;
      if ( r( x, x0, r0 ) <= 1. )
         return ( 1 + std::cos( walberla::math::pi * r( x, x0, r0 ) ) ) * 0.25;
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > slottedCylinder = [&]( const hyteg::Point3D& x ) -> real_t {
      const Point3D x0( 0.0, 0.75, 0.0 );
      const real_t  r0 = 0.15;
      if ( ( r( x, x0, r0 ) <= 1. ) && ( std::abs( x[0] - x0[0] ) >= 0.025 || x[1] >= 0.85 ) )
         return 1;
      else
         return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > initialBodies = [&]( const hyteg::Point3D& x ) -> real_t {
      return conicalBody( x ) + gaussianCone( x ) + slottedCylinder( x );
   };

   auto vel_x = []( const hyteg::Point3D& x ) -> real_t { return -x[1]; };

   auto vel_y = []( const hyteg::Point3D& x ) -> real_t { return x[0]; };

   writeDomainPartitioningVTK( storage, "../../output", "MMOC2DCircularConvectionBlendingTest_Domain" );

   FunctionType c( "c", storage, minLevel, maxLevel );
   FunctionType cInitial( "cInitial", storage, minLevel, maxLevel );
   FunctionType cError( "cError", storage, minLevel, maxLevel );
   FunctionType cMass( "cError", storage, minLevel, maxLevel );
   FunctionType velocityMagnitude( "velocityMagnitude", storage, minLevel, maxLevel );
   typename FunctionTrait< FunctionType >::AssocVectorFunctionType uv( "uv", storage, minLevel, maxLevel );

   FunctionType tmp0( "tmp0", storage, minLevel, maxLevel );
   FunctionType tmp1( "tmp1", storage, minLevel, maxLevel );

   MassOperator                  M( storage, minLevel, maxLevel );
   MMOCTransport< FunctionType > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   uv.interpolate( {vel_x, vel_y}, maxLevel );
   c.interpolate( initialBodies, maxLevel );
   cInitial.interpolate( initialBodies, maxLevel );

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > magnitude = []( const Point3D&,
                                                                                           const std::vector< real_t >& values ) {
      return Point3D( values[0], values[1], 0 ).norm();
   };

   velocityMagnitude.interpolate( magnitude, {uv[0], uv[1]}, maxLevel, All );

   VTKOutput vtkOutput( "../../output", "MMOC2DCircularConvectionBlendingTest", storage );

   vtkOutput.add( uv );
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
   {
      vtkOutput.write( maxLevel );
   }

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %10d | %8d | %15.3e | %10.3e | %30.2f%% ", 0, 0, max_temp, total_mass, 0, 0. ) )

   for ( uint_t i = 1; i <= outerSteps; i++ )
   {
      transport.step( c, uv, uv, maxLevel, Inner, dt, innerSteps, i == 1 );

      cError.assign( {1.0, -1.0}, {c, cInitial}, maxLevel, All );
      max_temp = c.getMaxMagnitude( maxLevel, All );
      M.apply( c, cMass, maxLevel, All );
      auto total_mass_new  = cMass.sumGlobal( maxLevel );
      auto total_mass_lost = 1.0 - ( total_mass_new / total_mass );
      // WALBERLA_CHECK_LESS( std::abs( total_mass_lost ), 1e-12 );
      total_mass = total_mass_new;

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          " %10d | %8d | %15.3e | %10.3e | %30.2f%% ", i, i * innerSteps, max_temp, total_mass, total_mass_lost * 100. ) )

      if ( vtk )
      {
         vtkOutput.write( maxLevel, i );
      }
   }

   auto error_E1 = errorE1( maxLevel, c, cInitial, tmp0, tmp1, M );
   auto error_E2 = errorE2( maxLevel, c, cInitial, tmp0, tmp1, M );

   WALBERLA_LOG_INFO_ON_ROOT( "E1: " << walberla::format( "%5.4e", error_E1 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "E2: " << walberla::format( "%5.4e", error_E2 ) );

   WALBERLA_CHECK_LESS( error_E1, errorE1Limit );
   WALBERLA_CHECK_LESS( error_E2, errorE2Limit );

   storage->getTimingTree()->stop( "Total" );
   WALBERLA_LOG_INFO_ON_ROOT( storage->getTimingTree()->getCopyWithRemainder() );
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::test< hyteg::P1Function< real_t >, hyteg::P1ElementwiseBlendingMassOperator >( 6, 2.2e-04, 1.2e-03 );
   hyteg::test< hyteg::P2Function< real_t >, hyteg::P2ElementwiseBlendingMassOperator >( 5, 9.0e-07, 6.0e-06 );
}
