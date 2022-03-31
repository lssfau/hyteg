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

using walberla::int_c;

#if 0
static std::string getDateTimeID()
{
   std::vector< char > cTimeString( 64 );
   WALBERLA_ROOT_SECTION()
   {
      std::time_t t;
      std::time( &t );
      std::strftime( cTimeString.data(), 64, "%F_%H-%M-%S", std::localtime( &t ) );
   }

   walberla::mpi::broadcastObject( cTimeString );

   std::string timeString( cTimeString.data() );
   return timeString;
}
#endif

template < typename FunctionType_T, typename LaplaceOperator_T, typename MassOperator_T, typename UnsteadyDiffusionOperator_T >
void solve( MeshInfo&               meshInfo,
            bool                    setBlendingMap,
            Solution&               solution,
            Solution&               velocityX,
            Solution&               velocityY,
            Solution&               velocityZ,
            real_t                  dt,
            real_t                  diffusivity,
            uint_t                  level,
            DiffusionTimeIntegrator diffusionTimeIntegrator,
            bool                    enableDiffusion,
            bool                    strangSplitting,
            bool                    resetParticles,
            uint_t                  resetParticlesInterval,
            bool                    adjustedAdvection,
            bool                    globalMaxLimiter,
            bool                    setParticlesOutsideDomainToZero,
            uint_t                  numTimeSteps,
            LoadBalancingOptions    lbOptions,
            bool                    vtk,
            bool                    vtkOutputVelocity,
            const std::string&      benchmarkName,
            uint_t                  printInterval,
            uint_t                  vtkInterval,
            bool                    verbose,
            std::string             dbFile )
{
   walberla::WcTimer localTimer;

   const bool outputTimingJSON = true;

   std::shared_ptr< PrimitiveStorage > storage;

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Creating storage ..." );
      printCurrentMemoryUsage();
   }

   {
      auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
          meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      switch ( lbOptions.type )
      {
      case 1:
         loadbalancing::greedy( *setupStorage );
         break;
      case 2:
         loadbalancing::roundRobinVolume( *setupStorage );
         break;
      default:
         break;
      }

      setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      if ( setBlendingMap )
      {
         if ( setupStorage->getNumberOfCells() == 0 )
         {
            AnnulusMap::setMap( *setupStorage );
         }
         else
         {
            IcosahedralShellMap::setMap( *setupStorage );
         }
      }

      storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );
      meshInfo.clear();
   }

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Creating storage ... done." );
      printCurrentMemoryUsage();
   }

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, "vtk/", benchmarkName + "_domain" );
   }

   auto timer = storage->getTimingTree();
   timer->start( "Setup" );

   const uint_t unknowns = numberOfGlobalDoFs< typename FunctionType_T::Tag >( *storage, level );
   const real_t hMin     = MeshQuality::getMinimalEdgeLength( storage, level );
   const real_t hMax     = MeshQuality::getMaximalEdgeLength( storage, level );

   if ( !resetParticles )
   {
      resetParticlesInterval = numTimeSteps;
   }

   const bool forcedParticleReset = adjustedAdvection || enableDiffusion;
   if ( forcedParticleReset )
   {
      resetParticles         = true;
      resetParticlesInterval = 1;
   }

   const std::string elementType = ( std::is_same< FunctionType_T, P1Function< real_t > >::value ? "P1" : "P2" );

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << benchmarkName )
   WALBERLA_LOG_INFO_ON_ROOT( " - time discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dt:                                           " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( "   + time steps:                                   " << numTimeSteps )
   WALBERLA_LOG_INFO_ON_ROOT( "   + time final:                                   " << real_c( numTimeSteps ) * dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + element type:                                 " << elementType )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + level:                                        " << level )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns (== particles), including boundary:  " << unknowns )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( "   + blending:                                     " << ( setBlendingMap ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( " - advection-diffusion settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + diffusivity:                                  "
                              << ( enableDiffusion ? std::to_string( diffusivity ) : "disabled (== 0)" ) )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + diffusion time integrator:                    "
       << ( enableDiffusion ?
                ( diffusionTimeIntegrator == DiffusionTimeIntegrator::ImplicitEuler ? "implicit Euler" : "Crank-Nicolson" ) :
                "disabled" ) )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + Strang splitting:                             " << ( enableDiffusion && strangSplitting ? "yes" : "disabled" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + adjusted advection:                           " << ( adjustedAdvection ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + particle reset:                               "
                              << ( resetParticles ? "yes" : "no" ) << ( forcedParticleReset ? " (forced)" : "" ) )
   if ( resetParticles )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + particle reset interval:                      " << resetParticlesInterval );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + global max limiter:                           " << ( globalMaxLimiter ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK:                                          " << ( vtk ? "yes" : "no" ) )
   if ( vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + VTK interval:                                 " << vtkInterval )
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + print interval:                               " << printInterval )
   WALBERLA_LOG_INFO_ON_ROOT( "   + database file:                                " << dbFile )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   const auto domainInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( domainInfo );

   FixedSizeSQLDB db( dbFile );

   const std::string gitsha = gitSHA1();

   db.setVariableEntry( "ts", uint_c( 0 ) );

   db.setConstantEntry( "git_sha", gitsha );
   db.setConstantEntry( "element_type", elementType );
   db.setConstantEntry( "dt", dt );
   db.setConstantEntry( "num_ts", numTimeSteps );
   db.setConstantEntry( "level", level );
   db.setConstantEntry( "unknowns", unknowns );
   db.setConstantEntry( "h_min", hMin );
   db.setConstantEntry( "h_max", hMax );
   db.setConstantEntry( "num_macro_cells", storage->getNumberOfGlobalCells() );
   db.setConstantEntry( "num_macro_faces", storage->getNumberOfGlobalFaces() );
   db.setConstantEntry( "num_macro_edges", storage->getNumberOfGlobalEdges() );
   db.setConstantEntry( "num_macro_vertices", storage->getNumberOfGlobalVertices() );
   db.setConstantEntry( "num_macro_primitives", storage->getNumberOfGlobalPrimitives() );
   db.setConstantEntry( "diffusivity", diffusivity );
   db.setConstantEntry( "pure_advection", !enableDiffusion );
   db.setConstantEntry( "strang_splitting", strangSplitting );
   db.setConstantEntry( "adjusted_advection", adjustedAdvection );
   db.setConstantEntry( "particle_reset", resetParticles );
   db.setConstantEntry( "particle_reset_interval", resetParticlesInterval );
   db.setConstantEntry( "global_max_limiter", globalMaxLimiter );

   typedef FunctionType_T              FunctionType;
   typedef LaplaceOperator_T           LaplaceOperator;
   typedef MassOperator_T              MassOperator;
   typedef UnsteadyDiffusionOperator_T UnsteadyDiffusionOperator;

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ..." );
      printCurrentMemoryUsage();
   }

   FunctionType c( "c", storage, level, level );
   FunctionType cOld( "cOld", storage, level, level );
   FunctionType cError( "cError", storage, level, level );
   FunctionType cSolution( "cSolution", storage, level, level );
   FunctionType cMass( "cMass", storage, level, level );
   FunctionType tmp( "tmp", storage, level, level );
   FunctionType tmp2( "tmp2", storage, level, level );

   typename FunctionTrait< FunctionType >::AssocVectorFunctionType uvw( "uvw", storage, level, level );
   typename FunctionTrait< FunctionType >::AssocVectorFunctionType uvwLast( "uvwLast", storage, level, level );
   // FunctionType u( "u", storage, level, level );
   // FunctionType v( "v", storage, level, level );
   // FunctionType w( "w", storage, level, level );
   // FunctionType uLast( "uLast", storage, level, level );
   // FunctionType vLast( "vLast", storage, level, level );
   // FunctionType wLast( "wLast", storage, level, level );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ... done." );
      printCurrentMemoryUsage();
   }

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Preparing operators and solvers ..." );
      printCurrentMemoryUsage();
   }

   const real_t                  diffusionDt = strangSplitting ? 0.5 * dt : dt;
   UnsteadyDiffusionOperator     diffusionOperator( storage, level, level, diffusionDt, diffusivity, diffusionTimeIntegrator );
   LaplaceOperator               L( storage, level, level );
   MassOperator                  M( storage, level, level );
   MMOCTransport< FunctionType > transport( storage, level, level, TimeSteppingScheme::RK4 );

   std::shared_ptr< Solver< UnsteadyDiffusionOperator > > solver;

   if ( enableDiffusion )
   {
      solver = std::make_shared< CGSolver< UnsteadyDiffusionOperator > >(
          storage, level, level, std::numeric_limits< uint_t >::max(), 1e-12 );
   }

   UnsteadyDiffusion< FunctionType, UnsteadyDiffusionOperator, LaplaceOperator, MassOperator > diffusionSolver(
       storage, level, level, solver );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Preparing operators and solvers ... done." );
      printCurrentMemoryUsage();
   }

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Interpolating solution ..." );
      printCurrentMemoryUsage();
   }

   c.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level );
   cSolution.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level );
   uvw[0].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityX ) ), level );
   uvw[1].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityY ) ), level );
   if ( storage->hasGlobalCells() )
   {
      uvw[2].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityZ ) ), level );
   }

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Interpolating solution ... done." );
      printCurrentMemoryUsage();
   }

   cError.assign( { 1.0, -1.0 }, { c, cSolution }, level, All );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Calculating relevant quantities ..." );
      printCurrentMemoryUsage();
   }

   auto       discrL2     = normL2( cError, tmp, M, level, Inner );
   auto       maxPeakDiff = maxPeakDifference( c, cSolution, level, All );
   auto       spuriousOsc = spuriousOscillations( c, level, All );
   auto       mass        = globalMass( c, tmp, M, level, All );
   const auto initialMass = mass;
   auto       massChange  = ( mass / initialMass ) - 1.0;
   real_t     timeTotal   = 0;
   real_t     vMax        = 0;

   vMax = velocityMaxMagnitude( uvw, tmp, tmp2, level, All );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Calculating relevant quantities ... done." );
      printCurrentMemoryUsage();
   }

   hyteg::VTKOutput vtkOutput( "./vtk", benchmarkName, storage, vtkInterval );

   if ( vtkOutputVelocity )
   {
      vtkOutput.add( uvw );
   }
   vtkOutput.add( c );
   vtkOutput.add( cSolution );
   vtkOutput.add( cError );

   if ( vtk )
      vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT(
       " timestep | time total | discr. L2 error | max peak diff. | spu. osc. | total mass | mass change | vel max mag " )
   WALBERLA_LOG_INFO_ON_ROOT(
       "----------+------------+-----------------+----------------+-----------+------------+-------------+------------ " )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8s | %10.5f | %15.3e | %14.3e | %9.3e | %10.3e | %10.2f%% | %10.3e ",
                                                "initial",
                                                timeTotal,
                                                discrL2,
                                                maxPeakDiff,
                                                spuriousOsc,
                                                mass,
                                                massChange * 100,
                                                vMax ) )

   db.setVariableEntry( "sim_time", timeTotal );
   db.setVariableEntry( "error_l2", discrL2 );
   db.setVariableEntry( "error_peak", maxPeakDiff );
   db.setVariableEntry( "spurious_osc", spuriousOsc );
   db.setVariableEntry( "mass", mass );
   db.setVariableEntry( "mass_change", massChange );
   db.setVariableEntry( "v_max", vMax );
   db.setVariableEntry( "run_time_advection", real_c( 0 ) );

   db.writeRowOnRoot();

   timer->stop( "Setup" );

   auto startTimeX = velocityX.currentTime();
   auto startTimeY = velocityY.currentTime();
   auto startTimeZ = velocityZ.currentTime();

   for ( uint_t i = 1; i <= numTimeSteps; i++ )
   {
      timer->start( "Simulation" );

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "timestep " << i )

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "interpolating velocity" )

      // To implement re-init intervals > 1 for time dependent velocity fields, we need to
      // reverse the the velocity field over the period in which the solution is transported
      // in the Lagrangian domain.
      auto lagrangeIntervalLength   = resetParticlesInterval;
      auto lagrangeIntervalIdx      = ( i - 1 ) / lagrangeIntervalLength;
      auto lagrangeIntervalInnerIdx = ( i - 1 ) % lagrangeIntervalLength;
      auto tsVelocityCurrent =
          lagrangeIntervalLength * lagrangeIntervalIdx + ( lagrangeIntervalLength - lagrangeIntervalInnerIdx ) - 1;
      auto tsVelocityNext = lagrangeIntervalLength * lagrangeIntervalIdx + ( lagrangeIntervalLength - lagrangeIntervalInnerIdx );

      velocityX.setTime( startTimeX + dt * real_c( tsVelocityCurrent ) );
      velocityY.setTime( startTimeY + dt * real_c( tsVelocityCurrent ) );
      velocityZ.setTime( startTimeZ + dt * real_c( tsVelocityCurrent ) );

      uvwLast[0].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityX ) ), level );
      uvwLast[1].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityY ) ), level );
      if ( storage->hasGlobalCells() )
      {
         uvwLast[2].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityZ ) ), level );
      }

      velocityX.setTime( startTimeX + dt * real_c( tsVelocityNext ) );
      velocityY.setTime( startTimeY + dt * real_c( tsVelocityNext ) );
      velocityZ.setTime( startTimeZ + dt * real_c( tsVelocityNext ) );

      uvw[0].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityX ) ), level );
      uvw[1].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityY ) ), level );
      if ( storage->hasGlobalCells() )
      {
         uvw[2].interpolate( std::function< real_t( const Point3D& ) >( std::ref( velocityZ ) ), level );
      }

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "performing transport step" )

      real_t advectionTimeStepRunTime;

      vMax = velocityMaxMagnitude( uvw, tmp, tmp2, level, All );

      if ( enableDiffusion && strangSplitting )
      {
         solution.incTime( 0.5 * dt );

         cOld.assign( { 1.0 }, { c }, level, All );

         c.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level, DirichletBoundary );

         if ( verbose )
            WALBERLA_LOG_INFO_ON_ROOT( "performing diffusion step" )
         diffusionSolver.step( diffusionOperator, L, M, c, cOld, level, Inner );
      }

      if ( adjustedAdvection )
      {
         const real_t adjustedAdvectionPertubation = 0.1 * ( hMin / vMax );
         localTimer.start();
         transport.step( c,
                         uvw,
                         uvwLast,
                         level,
                         Inner,
                         dt,
                         1,
                         M,
                         0.0,
                         adjustedAdvectionPertubation,
                         globalMaxLimiter,
                         setParticlesOutsideDomainToZero );
         localTimer.end();
         advectionTimeStepRunTime = localTimer.last();
      }
      else
      {
         localTimer.start();
         transport.step( c,
                         uvw,
                         uvwLast,
                         level,
                         Inner,
                         dt,
                         1,
                         i == 1 || ( resetParticles && ( i - 1 ) % resetParticlesInterval == 0 ),
                         globalMaxLimiter,
                         setParticlesOutsideDomainToZero );
         localTimer.end();
         advectionTimeStepRunTime = localTimer.last();
      }

      if ( strangSplitting )
      {
         solution.incTime( 0.5 * dt );
      }
      else
      {
         solution.incTime( dt );
      }

      cOld.assign( { 1.0 }, { c }, level, All );

      c.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level, DirichletBoundary );

      if ( enableDiffusion )
      {
         if ( verbose )
            WALBERLA_LOG_INFO_ON_ROOT( "performing diffusion step" )
         diffusionSolver.step( diffusionOperator, L, M, c, cOld, level, Inner );
      }

      cSolution.interpolate( std::function< real_t( const Point3D& ) >( std::ref( solution ) ), level );
      cError.assign( { 1.0, -1.0 }, { c, cSolution }, level, All );

      timeTotal += dt;

      if ( ( printInterval == 0 && i == numTimeSteps ) || ( printInterval > 0 && i % printInterval == 0 ) )
      {
         discrL2     = normL2( cError, tmp, M, level, Inner );
         maxPeakDiff = maxPeakDifference( c, cSolution, level, All );
         spuriousOsc = spuriousOscillations( c, level, All );
         mass        = globalMass( c, tmp, M, level, All );
         massChange  = ( mass / initialMass ) - 1.0;

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8d | %10.5f | %15.3e | %14.3e | %9.3e | %10.3e | %10.2f%% | %10.3e",
                                                      i,
                                                      timeTotal,
                                                      discrL2,
                                                      maxPeakDiff,
                                                      spuriousOsc,
                                                      mass,
                                                      massChange * 100,
                                                      vMax ) )
      }

      if ( vtk )
         vtkOutput.write( level, i );

      db.setVariableEntry( "ts", i );
      db.setVariableEntry( "sim_time", timeTotal );
      db.setVariableEntry( "error_l2", discrL2 );
      db.setVariableEntry( "error_peak", maxPeakDiff );
      db.setVariableEntry( "spurious_osc", spuriousOsc );
      db.setVariableEntry( "mass", mass );
      db.setVariableEntry( "mass_change", massChange );
      db.setVariableEntry( "run_time_advection", advectionTimeStepRunTime );
      db.setVariableEntry( "v_max", vMax );

      db.writeRowOnRoot();

      timer->stop( "Simulation" );
      if ( outputTimingJSON )
      {
         writeTimingTreeJSON( *timer, dbFile + walberla::format( "_Timing_ts_%04d.json", i ) );
      }
   }
}

template void
    solve< P1Function< real_t >, P1ConstantLaplaceOperator, P1ConstantMassOperator, P1ConstantUnsteadyDiffusionOperator >(
        MeshInfo&               meshInfo,
        bool                    setBlendingMap,
        Solution&               solution,
        Solution&               velocityX,
        Solution&               velocityY,
        Solution&               velocityZ,
        real_t                  dt,
        real_t                  diffusivity,
        uint_t                  level,
        DiffusionTimeIntegrator diffusionTimeIntegrator,
        bool                    enableDiffusion,
        bool                    strangSplitting,
        bool                    resetParticles,
        uint_t                  resetParticlesInterval,
        bool                    adjustedAdvection,
        bool                    globalMaxLimiter,
        bool                    setParticlesOutsideDomainToZero,
        uint_t                  numTimeSteps,
        LoadBalancingOptions    lbOptions,
        bool                    vtk,
        bool                    vtkOutputVelocity,
        const std::string&      benchmarkName,
        uint_t                  printInterval,
        uint_t                  vtkInterval,
        bool                    verbose,
        std::string             dbFile );

template void solve< P2Function< real_t >,
                     P2ElementwiseBlendingLaplaceOperator,
                     P2ElementwiseBlendingMassOperator,
                     P2ElementwiseUnsteadyDiffusionOperator >( MeshInfo&               meshInfo,
                                                               bool                    setBlendingMap,
                                                               Solution&               solution,
                                                               Solution&               velocityX,
                                                               Solution&               velocityY,
                                                               Solution&               velocityZ,
                                                               real_t                  dt,
                                                               real_t                  diffusivity,
                                                               uint_t                  level,
                                                               DiffusionTimeIntegrator diffusionTimeIntegrator,
                                                               bool                    enableDiffusion,
                                                               bool                    strangSplitting,
                                                               bool                    resetParticles,
                                                               uint_t                  resetParticlesInterval,
                                                               bool                    adjustedAdvection,
                                                               bool                    globalMaxLimiter,
                                                               bool                    setParticlesOutsideDomainToZero,
                                                               uint_t                  numTimeSteps,
                                                               LoadBalancingOptions    lbOptions,
                                                               bool                    vtk,
                                                               bool                    vtkOutputVelocity,
                                                               const std::string&      benchmarkName,
                                                               uint_t                  printInterval,
                                                               uint_t                  vtkInterval,
                                                               bool                    verbose,
                                                               std::string             dbFile );

} // namespace moc_benchmarks
} // namespace hyteg
