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

#include "Helpers.hpp"

namespace hyteg {
namespace moc_benchmarks {

void solve( const MeshInfo&         meshInfo,
            bool                    setAnnulusMap,
            Solution&               solution,
            Solution&               velocityX,
            Solution&               velocityY,
            real_t                  dt,
            real_t                  diffusivity,
            uint_t                  level,
            DiffusionTimeIntegrator diffusionTimeIntegrator,
            bool                    enableDiffusion,
            bool                    resetParticles,
            uint_t                  numTimeSteps,
            bool                    vtk,
            const std::string&      benchmarkName,
            uint_t                  printInterval,
            uint_t                  vtkInterval )
{
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   if ( setAnnulusMap )
   {
      AnnulusMap::setMap( *setupStorage );
   }
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   const uint_t unknowns = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
   const real_t hMin     = MeshQuality::getMinimalEdgeLength( storage, level );
   const real_t hMax     = MeshQuality::getMaximalEdgeLength( storage, level );

   typedef P2Function< real_t >                   FunctionType;
   typedef P2ElementwiseBlendingLaplaceOperator   LaplaceOperator;
   typedef P2ElementwiseBlendingMassOperator      MassOperator;
   typedef P2ElementwiseUnsteadyDiffusionOperator UnsteadyDiffusionOperator;

   FunctionType c( "c", storage, level, level );
   FunctionType cOld( "cOld", storage, level, level );
   FunctionType cError( "cError", storage, level, level );
   FunctionType cSolution( "cSolution", storage, level, level );
   FunctionType cMass( "cMass", storage, level, level );
   FunctionType tmp( "tmp", storage, level, level );
   FunctionType u( "u", storage, level, level );
   FunctionType v( "v", storage, level, level );
   FunctionType w( "w", storage, level, level );

   UnsteadyDiffusionOperator     diffusionOperator( storage, level, level, dt, diffusivity, diffusionTimeIntegrator );
   LaplaceOperator               L( storage, level, level );
   MassOperator                  M( storage, level, level );
   MMOCTransport< FunctionType > transport( storage, setupStorage, level, level, TimeSteppingScheme::RK4 );

   auto cgSolver = std::make_shared< CGSolver< P2ElementwiseUnsteadyDiffusionOperator > >( storage, level, level );

   UnsteadyDiffusion< FunctionType, UnsteadyDiffusionOperator, LaplaceOperator, MassOperator > diffusionSolver(
       storage, level, level, cgSolver );

   c.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level );
   cSolution.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level );
   u.interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityX ) ), level );
   v.interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityY ) ), level );

   cError.assign( {1.0, -1.0}, {c, cSolution}, level, All );

   auto       discrL2     = normL2( cError, tmp, M, level, Inner );
   auto       maxPeakDiff = maxPeakDifference( c, cSolution, level, All );
   auto       mass        = globalMass( c, tmp, M, level, All );
   const auto initialMass = mass;
   auto       massChange  = ( mass / initialMass ) - 1.0;
   real_t     timeTotal   = 0;

   hyteg::VTKOutput vtkOutput( "./output", benchmarkName, storage, vtkInterval );

   vtkOutput.add( u );
   vtkOutput.add( v );
   vtkOutput.add( c );
   vtkOutput.add( cSolution );
   vtkOutput.add( cError );

   if ( vtk )
      vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << benchmarkName )
   WALBERLA_LOG_INFO_ON_ROOT( " - time discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dt:                                           " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( "   + time steps:                                   " << numTimeSteps )
   WALBERLA_LOG_INFO_ON_ROOT( "   + time final:                                   " << real_c( numTimeSteps ) * dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + level:                                        " << level )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns (== particles), including boundary:  " << unknowns )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( "   + annulus blending:                             " << ( setAnnulusMap ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( " - advection-diffusion settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + diffusivity:                                  "
                              << ( enableDiffusion ? std::to_string( diffusivity ) : "disabled (== 0)" ) )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + diffusion time integrator:                    "
       << ( enableDiffusion ?
                ( diffusionTimeIntegrator == DiffusionTimeIntegrator::ImplicitEuler ? "implicit Euler" : "Crank-Nicolson" ) :
                "disabled" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + particle reset:                               " << ( resetParticles ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK:                                          " << ( vtk ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + print interval:                               " << printInterval )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   WALBERLA_LOG_INFO_ON_ROOT( " timestep | time total | discr. L2 error | max peak diff. | total mass | mass change " )
   WALBERLA_LOG_INFO_ON_ROOT( "----------+------------+-----------------+----------------+------------+-------------" )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8s | %10.5f | %15.3e | %14.3e | %10.3e | %11.2f%% ",
                                                "initial",
                                                timeTotal,
                                                discrL2,
                                                maxPeakDiff,
                                                mass,
                                                massChange * 100 ) )

   for ( uint_t i = 1; i <= numTimeSteps; i++ )
   {
      cOld.assign( {1.0}, {c}, level, All );

      transport.step( c, u, v, w, level, Inner, dt, 1, i == 1 || resetParticles );

      timeTotal += dt;

      solution.incTime( dt );
      velocityX.incTime( dt );
      velocityY.incTime( dt );

      cSolution.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level );
      u.interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityX ) ), level );
      v.interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityY ) ), level );

      c.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level, DirichletBoundary );

      if ( enableDiffusion )
      {
         diffusionSolver.step( diffusionOperator, L, M, c, cOld, level, Inner );
      }

      cError.assign( {1.0, -1.0}, {c, cSolution}, level, All );

      discrL2     = normL2( cError, tmp, M, level, Inner );
      maxPeakDiff = maxPeakDifference( c, cSolution, level, All );
      mass        = globalMass( c, tmp, M, level, All );
      massChange  = ( mass / initialMass ) - 1.0;

      if ( ( printInterval == 0 && i == numTimeSteps ) || ( printInterval > 0 && i % printInterval == 0 ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8d | %10.5f | %15.3e | %14.3e | %10.3e | %11.2f%% ",
                                                      i,
                                                      timeTotal,
                                                      discrL2,
                                                      maxPeakDiff,
                                                      mass,
                                                      massChange * 100 ) )
      }

      if ( vtk )
         vtkOutput.write( level, i );
   }
}

} // namespace moc_benchmarks
} // namespace hyteg