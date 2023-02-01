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

#include <utility>

#include "hyteg/MeshQuality.hpp"
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

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

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
   , p0_(std::move( p0 ))
   , t_( t0 )
   {}

   real_t operator()( const Point3D& p )
   {
      // auto x_hat    = p0_[0] * std::cos( t_ ) - p0_[1] * std::sin( t_ );
      // auto y_hat    = -p0_[0] * std::sin( t_ ) + p0_[1] * std::cos( t_ );
      real_t x_hat    = p0_[0] + t_;
      real_t y_hat    = real_c( 0 );
      real_t exponent = -std::pow( r( p, x_hat, y_hat ), real_c( 2 ) ) / ( real_c( 4.0 ) * diffusivity_ * t_ );
      return ( real_c( 1.0 ) / ( real_c( 4.0 ) * pi * diffusivity_ * t_ ) ) * std::exp( exponent );
   }

   static real_t r( const hyteg::Point3D& p, const real_t& x_hat, const real_t& y_hat )
   {
      return std::sqrt( std::pow( p[0] - x_hat, real_c( 2 ) ) + std::pow( p[1] - y_hat, real_c( 2 ) ) );
   }

   void inc( real_t dt ) { t_ += dt; }

 private:
   real_t  diffusivity_;
   Point3D p0_;
   real_t  t_;
};


real_t errorE1( const uint_t&                 level,
                const P2Function< real_t >&   c,
                const P2Function< real_t >&   solution,
                const P2Function< real_t >&   tmp0,
                const P2Function< real_t >&   tmp1,
                const P2ConstantMassOperator& mass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E1 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::abs( values[0] - values[1] ); };

   tmp0.interpolate( E1, {solution, c}, level, All );
   mass.apply( tmp0, tmp1, level, All );
   return tmp1.sumGlobal( level, All );
}

real_t errorE2( const uint_t&                 level,
                const P2Function< real_t >&   c,
                const P2Function< real_t >&   solution,
                const P2Function< real_t >&   tmp0,
                const P2Function< real_t >&   tmp1,
                const P2ConstantMassOperator& mass )
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > E2 =
       []( const Point3D&, const std::vector< real_t >& values ) { return std::pow( std::abs( values[0] - values[1] ), 2 ); };

   tmp0.interpolate( E2, {solution, c}, level, All );
   mass.apply( tmp0, tmp1, level, All );
   return std::sqrt( tmp1.sumGlobal( level, All ) );
}

void runTest( uint_t maxLevel, uint_t steps, uint_t timeSteppingScheme, std::map<uint_t, std::pair<real_t,real_t>> checks )
{
   MeshInfo meshInfo     = hyteg::MeshInfo::meshRectangle( Point2D( {-1, -1} ), Point2D( {5, 1} ), MeshInfo::CRISS, 6, 1 );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   storage->getTimingTree()->start( "Entire test" );

   writeDomainPartitioningVTK( storage, "../../output", "P2UnsteadyConvectionDiffusion2DTest_Domain" );

   const uint_t  minLevel = 2;
   const real_t  hMin     = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t  hMax     = MeshQuality::getMaximalEdgeLength( storage, maxLevel );
   const real_t  t0       = 0.5 * pi;
   const Point3D p0( {-0.5, 0, 0} );
   const bool    vtk   = false;
   const real_t  tFinal = 1.0;
   const real_t dt = tFinal / real_c(steps);

   const real_t diffusivity = 1e-00;

   WALBERLA_LOG_INFO_ON_ROOT( "Convection-diffusion test" )
   WALBERLA_LOG_INFO_ON_ROOT( " - dt:                                           " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - h (min, max):                                 " << hMin << ", " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - level:                                        " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( " - diffusivity:                                  " << diffusivity )
   WALBERLA_LOG_INFO_ON_ROOT( " - time stepping scheme diffusion:               " << timeSteppingScheme )

   DiffusionTimeIntegrator timeIntegrator = timeSteppingScheme == 0 ? DiffusionTimeIntegrator::ImplicitEuler : DiffusionTimeIntegrator::CrankNicolson;

   Solution solution( diffusivity, p0, t0 );

   auto vel_x = []( const hyteg::Point3D& ) -> real_t { return 1; };
   auto vel_y = []( const hyteg::Point3D& ) -> real_t { return 0; };

   typedef P2Function< real_t >                FunctionType;
   typedef P2ConstantLaplaceOperator           LaplaceOperator;
   typedef P2ConstantMassOperator              MassOperator;
   typedef P2ConstantUnsteadyDiffusionOperator UnsteadyDiffusionOperator;

   FunctionType c( "c", storage, minLevel, maxLevel );
   FunctionType cOld( "cOld", storage, minLevel, maxLevel );
   FunctionType cError( "cError", storage, minLevel, maxLevel );
   FunctionType cSolution( "cSolution", storage, minLevel, maxLevel );
   FunctionType cMass( "cMass", storage, minLevel, maxLevel );
   FunctionType tmp0( "tmp0", storage, minLevel, maxLevel );
   FunctionType tmp1( "tmp1", storage, minLevel, maxLevel );
   FunctionType f( "f", storage, minLevel, maxLevel );
   FunctionType fOld( "fOld", storage, minLevel, maxLevel );
   typename FunctionTrait< FunctionType >::AssocVectorFunctionType uv( "uv", storage, minLevel, maxLevel );

   UnsteadyDiffusionOperator     diffusionOperator( storage, minLevel, maxLevel, dt, diffusivity, timeIntegrator );
   LaplaceOperator               L( storage, minLevel, maxLevel );
   MassOperator                  M( storage, minLevel, maxLevel );
   MMOCTransport< FunctionType > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   auto coarseGridSolver = std::make_shared< CGSolver< P2ConstantUnsteadyDiffusionOperator > >( storage, minLevel, maxLevel );
   auto smoother         = std::make_shared< GaussSeidelSmoother< P2ConstantUnsteadyDiffusionOperator > >();
   auto restriction      = std::make_shared< P2toP2QuadraticRestriction >();
   auto prolongation     = std::make_shared< P2toP2QuadraticProlongation >();
   auto solver           = std::make_shared< GeometricMultigridSolver< P2ConstantUnsteadyDiffusionOperator > >(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, 2, 2 );
   auto solverLoop = std::make_shared< SolverLoop< P2ConstantUnsteadyDiffusionOperator > >( solver, 1 );

   UnsteadyDiffusion< P2Function< real_t >, P2ConstantUnsteadyDiffusionOperator, P2ConstantLaplaceOperator, P2ConstantMassOperator > diffusionSolver(
       storage, minLevel, maxLevel, solverLoop );

   uv.interpolate( { vel_x, vel_y }, maxLevel );
   c.interpolate( solution, maxLevel );

   hyteg::VTKOutput vtkOutput( "../../output", "P2UnsteadyConvectionDiffusion2DTest", storage );

   vtkOutput.add( uv );
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

   real_t total_mass_lost = 0;

   auto output = [&](real_t discrL2Eps, real_t peakEps, uint_t steps){
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          " %8d | %15.3e | %10.3e | %10.3e | %30.2f%% ", steps, discrL2, E_peak, total_mass, total_mass_lost * 100. ) )

      auto error_E1_consistent = errorE1( maxLevel, c, cSolution, tmp0, tmp1, M );
      auto error_E2_consistent = errorE2( maxLevel, c, cSolution, tmp0, tmp1, M );

      WALBERLA_LOG_INFO_ON_ROOT( "E1 consistent: " << walberla::format( "%5.4e", error_E1_consistent ) );
      WALBERLA_LOG_INFO_ON_ROOT( "E2 consistent: " << walberla::format( "%5.4e", error_E2_consistent ) );

      WALBERLA_CHECK_LESS( std::abs( E_peak ), peakEps )
      WALBERLA_CHECK_LESS( discrL2, discrL2Eps )
   };

   for ( uint_t i = 1; i <= steps; i++ )
   {
      solution.inc( dt );

      cSolution.interpolate( solution, maxLevel, All );

      transport.step( c, uv, uv, maxLevel, Inner, dt, 1 );

      fOld.assign( {1.0}, {f}, maxLevel, All );
      cOld.assign( {1.0}, {c}, maxLevel, All );

      c.interpolate( solution, maxLevel, DirichletBoundary );

      diffusionSolver.step( diffusionOperator, L, M, c, cOld, f, fOld, maxLevel, Inner );

      // various errors
      cError.assign( { 1.0, -1.0 }, { c, cSolution }, maxLevel, All );

      discrL2            = std::sqrt( cError.dotGlobal( cError, maxLevel, All ) /
                           real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
      maxTempApproximate = c.getMaxMagnitude( maxLevel, All );
      maxTempAnalytical  = cSolution.getMaxMagnitude( maxLevel, All );
      E_peak             = maxTempApproximate / maxTempAnalytical - 1;

      M.apply( c, cMass, maxLevel, All );
      real_t total_mass_new = cMass.sumGlobal( maxLevel );
      total_mass_lost       = real_c( 1.0 ) - ( total_mass_new / total_mass );
      total_mass            = total_mass_new;

      if (checks.find(i) != checks.end()){
         output(checks[i].first,checks[i].second, i);
      }
      if ( vtk )
         vtkOutput.write( maxLevel, i );
   }



   storage->getTimingTree()->stop( "Entire test" );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   uint_t timeSteppingScheme = 0;
   uint_t maxLevel           = 4;
   auto dp = std::is_same<real_t, double>();
   auto checks = std::map<uint_t, std::pair<real_t,real_t>>();
   checks[10] = {real_c( 1.4e-04 ), real_c( 8.7e-03 )};
   checks[20] = {real_c( 9.0e-05 ), real_c( 3.6e-03 )};
   checks[30] = {real_c( dp ? 3e-05: 4.6e-05 ), real_c( dp ? 1.5e-03 : 3e-03 )};
   checks[40] = {real_c( dp ? 2e-05 : 4.5e-05 ), real_c( dp ? 2.0e-03 : 3e-03 )};
   runTest( maxLevel, 40, timeSteppingScheme, checks );

   return EXIT_SUCCESS;
}
