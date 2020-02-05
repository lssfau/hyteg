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

#include "hyteg/composites/MMOCTransport.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using walberla::math::pi;

using namespace hyteg;

/// Setting from Kuzmin transport tutorial section 4.4.6.3

class Solution
{
 public:
   Solution( real_t diffusivity, Point3D p0, real_t t0 )
   : diffusivity_( diffusivity )
   , p0_( p0 )
   , t_( t0 )
   {}

   real_t operator()( const Point3D& p )
   {
      auto x_hat    = p0_[0] * std::cos( t_ ) - p0_[1] * std::sin( t_ );
      auto y_hat    = -p0_[0] * std::sin( t_ ) + p0_[1] * std::cos( t_ );
      auto exponent = -std::pow( r( p, x_hat, y_hat ), 2 ) / ( 4.0 * diffusivity_ * t_ );
      return ( 1.0 / ( 4.0 * pi * diffusivity_ * t_ ) ) * std::exp( exponent );
   }

   real_t r( const hyteg::Point3D& p, const real_t& x_hat, const real_t& y_hat )
   {
      return std::sqrt( std::pow( p[0] - x_hat, 2 ) + std::pow( p[1] - y_hat, 2 ) );
   }

   void inc( real_t dt ) { t_ += dt; }

 private:
   real_t  diffusivity_;
   Point3D p0_;
   real_t  t_;
};

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo meshInfo = hyteg::MeshInfo::meshRectangle( Point2D( {-1, -1} ), Point2D( {1, 1} ), MeshInfo::CRISS, 3, 3 );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   storage->getTimingTree()->start( "Entire test" );

   writeDomainPartitioningVTK( storage, "../../output", "P2UnsteadyConvectionDiffusion2DTest_Domain" );

   real_t minLength = std::numeric_limits< real_t >::max();
   real_t maxLength = 0;
   for ( const auto& edge : setupStorage.getEdges() )
   {
      auto edgeLength = edge.second->getLength();
      if ( edgeLength < minLength )
         minLength = edgeLength;
      if ( edgeLength > maxLength )
         maxLength = edgeLength;
   }

   const uint_t  minLevel = 2;
   const uint_t  maxLevel = 5;
   const real_t  hMax     = maxLength / levelinfo::num_microedges_per_edge( maxLevel );
   const real_t  hMin     = minLength / levelinfo::num_microedges_per_edge( maxLevel );
   const real_t  dt       = 1e-3;
   const real_t  t0       = 0.5 * pi;
   const Point3D p0( {0, 0.5, 0} );
   const uint_t  steps           = 20;
   const real_t  stepsRevolution = 2 * pi / dt;
   const bool    vtk             = false;

   const real_t diffusivity = 1e-3;

   WALBERLA_LOG_INFO_ON_ROOT( "Circular convection-diffusion" )
   WALBERLA_LOG_INFO_ON_ROOT( " - dt:                                           " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - h (min, max):                                 " << hMin << ", " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - level:                                        " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( " - diffusivity:                                  " << diffusivity )
   WALBERLA_LOG_INFO_ON_ROOT( " - steps for revolution:                         " << stepsRevolution )

   Solution solution( diffusivity, p0, t0 );

   auto vel_x = []( const hyteg::Point3D& x ) -> real_t { return -x[1]; };

   auto vel_y = []( const hyteg::Point3D& x ) -> real_t { return x[0]; };

   typedef P2Function< real_t >        FunctionType;
   typedef P2ConstantMassOperator      MassOperator;
   typedef P2UnsteadyDiffusionOperator UnsteadyDiffusionOperator;

   FunctionType c( "c", storage, minLevel, maxLevel );
   FunctionType cError( "cError", storage, minLevel, maxLevel );
   FunctionType cSolution( "cSolution", storage, minLevel, maxLevel );
   FunctionType cMass( "cMass", storage, minLevel, maxLevel );
   FunctionType u( "u", storage, minLevel, maxLevel );
   FunctionType v( "v", storage, minLevel, maxLevel );
   FunctionType w( "w", storage, minLevel, maxLevel );
   FunctionType f( "f", storage, minLevel, maxLevel );

   UnsteadyDiffusionOperator     diffusionOperator( storage, minLevel, maxLevel, dt, diffusivity );
   MassOperator                  M( storage, minLevel, maxLevel );
   MMOCTransport< FunctionType > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   auto coarseGridSolver = std::make_shared< CGSolver< P2UnsteadyDiffusionOperator > >( storage, minLevel, minLevel );
   auto smoother         = std::make_shared< GaussSeidelSmoother< P2UnsteadyDiffusionOperator > >();
   auto restriction      = std::make_shared< P2toP2QuadraticRestriction >();
   auto prolongation     = std::make_shared< P2toP2QuadraticProlongation >();
   auto solver           = std::make_shared< GeometricMultigridSolver< P2UnsteadyDiffusionOperator > >(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, 2, 2 );
   auto solverLoop = std::make_shared< SolverLoop< P2UnsteadyDiffusionOperator > >( solver, 1 );

   UnsteadyDiffusion< P2Function< real_t >, P2UnsteadyDiffusionOperator, P2ConstantMassOperator > diffusionSolver(
       storage, minLevel, maxLevel, solverLoop );

   u.interpolate( vel_x, maxLevel );
   v.interpolate( vel_y, maxLevel );
   c.interpolate( solution, maxLevel );

   hyteg::VTKOutput vtkOutput( "../../output", "P2UnsteadyConvectionDiffusion2DTest", storage );

   vtkOutput.add( u );
   vtkOutput.add( v );
   vtkOutput.add( c );
   vtkOutput.add( cSolution );
   vtkOutput.add( cError );

   WALBERLA_LOG_INFO_ON_ROOT( " timestep | discr. L2 error |     E_peak | total mass | mass lost since last outer step " )
   WALBERLA_LOG_INFO_ON_ROOT( "----------+-----------------+------------+------------+---------------------------------" )

   cSolution.interpolate( solution, maxLevel, All );

   // various errors
   cError.assign( {1.0, -1.0}, {c, cSolution}, maxLevel, All );

   auto discrL2            = std::sqrt( cError.dotGlobal( cError, maxLevel, All ) /
                             real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   auto maxTempApproximate = c.getMaxMagnitude( maxLevel, All );
   auto maxTempAnalytical  = cSolution.getMaxMagnitude( maxLevel, All );
   auto E_peak             = maxTempApproximate / maxTempAnalytical - 1;

   M.apply( c, cMass, maxLevel, All );
   auto total_mass = cMass.sumGlobal( maxLevel );

   if ( vtk )
      vtkOutput.write( maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( " %8d | %15.3e | %10.3e | %10.3e | %30.2f%% ", 0, discrL2, E_peak, total_mass, 0. ) )

   for ( uint_t i = 1; i <= steps; i++ )
   {
      solution.inc( dt );

      c.interpolate( solution, maxLevel, DirichletBoundary );
      cSolution.interpolate( solution, maxLevel, All );

      transport.step( c, u, v, w, maxLevel, Inner, dt, 1 );
      diffusionSolver.step( diffusionOperator, M, c, f, maxLevel, Inner );

      // various errors
      cError.assign( {1.0, -1.0}, {c, cSolution}, maxLevel, All );

      discrL2            = std::sqrt( cError.dotGlobal( cError, maxLevel, All ) /
                           real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
      maxTempApproximate = c.getMaxMagnitude( maxLevel, All );
      maxTempAnalytical  = cSolution.getMaxMagnitude( maxLevel, All );
      E_peak             = maxTempApproximate / maxTempAnalytical - 1;

      M.apply( c, cMass, maxLevel, All );
      auto total_mass_new  = cMass.sumGlobal( maxLevel );
      auto total_mass_lost = 1.0 - ( total_mass_new / total_mass );
      total_mass           = total_mass_new;

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          " %8d | %15.3e | %10.3e | %10.3e | %30.2f%% ", i, discrL2, E_peak, total_mass, total_mass_lost * 100. ) )

      WALBERLA_CHECK_LESS( std::abs( E_peak ), 5e-04 )
      WALBERLA_CHECK_LESS( discrL2, 2e-03 )

      if ( vtk )
         vtkOutput.write( maxLevel, i );
   }

   storage->getTimingTree()->stop( "Entire test" );

   // WALBERLA_LOG_INFO_ON_ROOT( storage->getTimingTree()->getCopyWithRemainder() );

   return EXIT_SUCCESS;
}
