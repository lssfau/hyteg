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
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "constantStencilOperator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

class SwirlVelocityX
{
 public:
   SwirlVelocityX( const uint_t& steps )
   : steps_( steps )
   , currentStep_( 0 )
   {}

   real_t operator()( const Point3D& p )
   {
      return std::sin( pi * p[0] ) * std::sin( pi * p[0] ) * std::sin( 2 * pi * p[1] ) *
             std::cos( pi * ( real_c( currentStep_ ) / real_c( steps_ ) ) );
   }

   void step() { currentStep_++; }

 private:
   uint_t steps_;
   uint_t currentStep_;
};

class SwirlVelocityY
{
 public:
   SwirlVelocityY( const uint_t& steps )
   : steps_( steps )
   , currentStep_( 0 )
   {}

   real_t operator()( const Point3D& p )
   {
      return -std::sin( pi * p[1] ) * std::sin( pi * p[1] ) * std::sin( 2 * pi * p[0] ) *
             std::cos( pi * ( real_c( currentStep_ ) / real_c( steps_ ) ) );
   }

   void step() { currentStep_++; }

 private:
   uint_t steps_;
   uint_t currentStep_;
};

/// Error definitions in section 4.4.6: Numerical examples

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

real_t errorE1( const uint_t&                   level,
                const P2Function< real_t >&     c,
                const P2Function< real_t >&     solution,
                const P2Function< real_t >&     tmp0,
                const P2Function< real_t >&     tmp1,
                const P2ConstantMassOperator &  mass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E1 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::abs( values[0] - values[1] ); };

   tmp0.interpolate( E1, {solution, c}, level, All );
   mass.apply( tmp0, tmp1, level, All );
   return tmp1.sumGlobal( level, All );
}

real_t errorE2( const uint_t&                   level,
                const P2Function< real_t >&     c,
                const P2Function< real_t >&     solution,
                const P2Function< real_t >&     tmp0,
                const P2Function< real_t >&     tmp1,
                const P2ConstantMassOperator&  mass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E2 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::pow( std::abs( values[0] - values[1] ), 2 ); };

   tmp0.interpolate( E2, {solution, c}, level, All );
   mass.apply( tmp0, tmp1, level, All );
   return std::sqrt( tmp1.sumGlobal( level, All ) );
}

/// Setting from Kuzmin transport tutorial section 4.4.6.1
int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo meshInfo = hyteg::MeshInfo::meshRectangle( Point2D( 0, 0 ), Point2D( 1, 1 ), MeshInfo::CRISS, 1, 1 );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   const uint_t minLevel = 2;
   const uint_t maxLevel = 5;
   real_t       dt       = 1e-01;
   const real_t tEnd     = 1.5;
   const uint_t steps    = uint_c( tEnd / dt );

   WALBERLA_LOG_INFO_ON_ROOT( "Swirling flow convection" )
   WALBERLA_LOG_INFO_ON_ROOT( " - dt:         " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - level:      " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( " - time steps: " << steps )
   WALBERLA_LOG_INFO_ON_ROOT( " - total time: " << tEnd )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   auto r = []( const hyteg::Point3D& x, const hyteg::Point3D& x0, const real_t& r0 ) -> real_t {
     return ( 1 / r0 ) * std::sqrt( std::pow( x[0] - x0[0], 2 ) + std::pow( x[1] - x0[1], 2 ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > conicalBody = [&]( const hyteg::Point3D& x ) -> real_t {
     const Point3D x0( 0.5, 0.25, 0.0 );
     const real_t  r0 = 0.15;
     if ( r( x, x0, r0 ) <= 1. )
        return 1 - r( x, x0, r0 );
     else
        return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > gaussianCone = [&]( const hyteg::Point3D& x ) -> real_t {
     const Point3D x0( 0.25, 0.5, 0.0 );
     const real_t  r0 = 0.15;
     if ( r( x, x0, r0 ) <= 1. )
        return ( 1 + std::cos( walberla::math::pi * r( x, x0, r0 ) ) ) * 0.25;
     else
        return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > slottedCylinder = [&]( const hyteg::Point3D& x ) -> real_t {
     const Point3D x0( 0.5, 0.75, 0.0 );
     const real_t  r0 = 0.15;
     if ( ( r( x, x0, r0 ) <= 1. ) && ( std::abs( x[0] - x0[0] ) >= 0.025 || x[1] >= 0.85 ) )
        return 1;
     else
        return 0.0;
   };

   std::function< real_t( const hyteg::Point3D& ) > initialBodies = [&]( const hyteg::Point3D& x ) -> real_t {
     return conicalBody( x ) + gaussianCone( x ) + slottedCylinder( x );
   };

   SwirlVelocityX vel_x( steps );
   SwirlVelocityY vel_y( steps );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   writeDomainPartitioningVTK( storage, "../../output", "MMOC2DSwirlingFlowConvectionTest_Domain" );

   typedef P2Function< real_t >   FunctionType;
   typedef P2ConstantMassOperator MassOperator;

   FunctionType c( "c", storage, minLevel, maxLevel );
   FunctionType cInitial( "cInitial", storage, minLevel, maxLevel );
   FunctionType cError( "cError", storage, minLevel, maxLevel );
   FunctionType cMass( "cMass", storage, minLevel, maxLevel );
   FunctionType tmp0( "tmp0", storage, minLevel, maxLevel );
   FunctionType                                                    tmp1( "tmp1", storage, minLevel, maxLevel );
   typename FunctionTrait< FunctionType >::AssocVectorFunctionType uv( "uv", storage, minLevel, maxLevel );
   typename FunctionTrait< FunctionType >::AssocVectorFunctionType uvLastTimeStep( "uvLast", storage, minLevel, maxLevel );

   MassOperator                  M( storage, minLevel, maxLevel );
   MMOCTransport< FunctionType > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   c.interpolate( initialBodies, maxLevel );
   cInitial.interpolate( initialBodies, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( " timestep | max temperature | total mass | mass lost since last outer step " )
   WALBERLA_LOG_INFO_ON_ROOT( "----------+-----------------+------------+---------------------------------" )

   real_t total_mass = 0.0;

   cError.assign( {1.0, -1.0}, {c, cInitial}, maxLevel, All );

//   vtkOutput.write( maxLevel, 0 );

   for ( uint_t i = 1; i <= steps; i++ )
   {
      uvLastTimeStep.interpolate( { vel_x, vel_y }, maxLevel );
      vel_x.step();
      vel_y.step();
      uv.interpolate( { vel_x, vel_y }, maxLevel );

      transport.step( c, uv, uvLastTimeStep, maxLevel, Inner, dt, 1, i == 1 );

      cError.assign( {1.0, -1.0}, {c, cInitial}, maxLevel, All );
      auto max_temp = c.getMaxMagnitude( maxLevel, All );
      M.apply( c, cMass, maxLevel, All );
      auto total_mass_new  = cMass.sumGlobal( maxLevel );
      auto total_mass_lost = 1.0 - ( total_mass_new / total_mass );
      total_mass           = total_mass_new;

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( " %8d | %15.3e | %10.3e | %30.2f%% ", i, max_temp, total_mass, total_mass_lost * 100. ) )

//      vtkOutput.write( maxLevel, i );
   }

   auto p2MassFormFenics =
       std::make_shared< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >();
   P2RowSumForm                       rowSumMass( p2MassFormFenics );
   P2ConstantOperator< P2RowSumForm > MLumped( storage, minLevel, maxLevel, rowSumMass );

   auto error_E1 = errorE1( maxLevel, c, cInitial, tmp0, tmp1, MLumped );
   auto error_E2 = errorE2( maxLevel, c, cInitial, tmp0, tmp1, MLumped );
   auto error_E1_consistent = errorE1( maxLevel, c, cInitial, tmp0, tmp1, M );
   auto error_E2_consistent = errorE2( maxLevel, c, cInitial, tmp0, tmp1, M );

   WALBERLA_LOG_INFO_ON_ROOT( "E1 lumped:     " << walberla::format( "%5.4e", error_E1 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "E1 consistent: " << walberla::format( "%5.4e", error_E1_consistent ) );

   WALBERLA_LOG_INFO_ON_ROOT( "E2 lumped:     " << walberla::format( "%5.4e", error_E2 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "E2 consistent: " << walberla::format( "%5.4e", error_E2_consistent ) );

   WALBERLA_CHECK_LESS( error_E1, 6.1e-05 );
   WALBERLA_CHECK_LESS( error_E2, 4.0e-04 );

   return EXIT_SUCCESS;
}
