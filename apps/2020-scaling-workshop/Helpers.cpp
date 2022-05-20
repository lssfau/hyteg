/*
 * Copyright (c) 2017-2020 Nils Kohl, Dominik Thoennes.
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
namespace scaling_workshop {

using walberla::int_c;

void writeDataHeader()
{
   WALBERLA_LOG_INFO_ON_ROOT(
       " iteration | iteration type | FMG level | L2 error velocity | L2 residual velocity | L2 error pressure | L2 residual pressure " );
   WALBERLA_LOG_INFO_ON_ROOT(
       "-----------+----------------+-----------+-------------------+----------------------+-------------------+----------------------" );
}

// iterationType: I: initial
//                V: v-cycle
//                F: FMG
void writeDataRow( uint_t          iteration,
                   std::string     iterationType,
                   uint_t          fmgLevel,
                   real_t          errorL2Velocity,
                   real_t          residualL2Velocity,
                   real_t          errorL2Pressure,
                   real_t          residualL2Pressure,
                   FixedSizeSQLDB& db )
{
   std::string fmgLevelString = "-";
   if ( iterationType == "F" )
   {
      fmgLevelString = walberla::format( "%9d", fmgLevel );
   }

   db.setVariableEntry( "iteration", iteration );
   db.setVariableEntry( "iteration_type", iterationType );
   db.setVariableEntry( "fmg_level", fmgLevel );
   db.setVariableEntry( "error_l2_velocity", errorL2Velocity );
   db.setVariableEntry( "residual_l2_velocity", residualL2Velocity );
   db.setVariableEntry( "error_l2_pressure", errorL2Pressure );
   db.setVariableEntry( "residual_l2_pressure", residualL2Pressure );

   db.writeRowOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9d | %14s | %9s | %17.5e | %20.5e | %17.5e | %20.5e ",
                                                iteration,
                                                iterationType.c_str(),
                                                fmgLevelString.c_str(),
                                                errorL2Velocity,
                                                residualL2Velocity,
                                                errorL2Pressure,
                                                residualL2Pressure ) );
}

template < template < typename > class StokesFunction,
           typename StokesOperator,
           typename VelocityMassOperator,
           typename PressureMassOperator,
           typename Restriction,
           typename Prolongation,
           typename FMGProlongation >
void solveRHS0Implementation( const std::shared_ptr< PrimitiveStorage >&              storage,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionU,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionV,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionW,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionP,
                              const std::function< real_t( const hyteg::Point3D& ) >& initialU,
                              const std::function< real_t( const hyteg::Point3D& ) >& initialV,
                              const std::function< real_t( const hyteg::Point3D& ) >& initialW,
                              const std::function< real_t( const hyteg::Point3D& ) >& initialP,
                              const std::function< real_t( const hyteg::Point3D& ) >& rhsU,
                              const std::function< real_t( const hyteg::Point3D& ) >& rhsV,
                              const std::function< real_t( const hyteg::Point3D& ) >& rhsW,
                              uint_t                                                  minLevel,
                              uint_t                                                  maxLevel,
                              MultigridSettings                                       multigridSettings,
                              SmootherSettings                                        smootherSettings,
                              CoarseGridSettings                                      coarseGridSettings,
                              bool                                                    projectPressure,
                              bool                                                    projectPressurefterRestriction,
                              bool                                                    vtk,
                              const std::string&                                      benchmarkName,
                              FixedSizeSQLDB                                          db,
                              std::string                                             timingFile,
                              bool                                                    RHSisZero )
{
   auto timer = storage->getTimingTree();
   timer->start( "Total" );
   timer->start( "Setup" );

   printGitInfo();

   const DoFType errorFlag = Inner | NeumannBoundary;

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, "vtk/", benchmarkName + "_domain" );
   }

   const uint_t unknowns = numberOfGlobalDoFs< typename StokesFunction< real_t >::Tag >(
       *storage, maxLevel, walberla::mpi::MPIManager::instance()->comm(), true );
   const real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel, true );
   const real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel, true );
   const auto   discretization =
       ( std::is_same< typename StokesFunction< real_t >::Tag, P2P1TaylorHoodFunctionTag >::value ? "P2-P1" : "P1-P1" );

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << benchmarkName )
   WALBERLA_LOG_INFO_ON_ROOT( " - parallelism: " )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + MPI processes:                                " << walberla::mpi::MPIManager::instance()->numProcesses() );
   WALBERLA_OPENMP_SECTION()
   {
      WALBERLA_LOG_INFO_ON_ROOT(
          "   + OpenMP threads per MPI process:               " << hyteg::OpenMPManager::instance()->numThreads() );
   }
   WALBERLA_NON_OPENMP_SECTION() { WALBERLA_LOG_INFO_ON_ROOT( "   + OpenMP disabled" ); }
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + elements:                                     " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + min level:                                    " << minLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + max level:                                    " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns, including boundary:                 " << unknowns )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - multigrid settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + pre smooth:                                   " << multigridSettings.preSmooth )
   WALBERLA_LOG_INFO_ON_ROOT( "   + post smooth:                                  " << multigridSettings.postSmooth )
   WALBERLA_LOG_INFO_ON_ROOT( "   + inc smooth:                                   " << multigridSettings.incSmooth )
   WALBERLA_LOG_INFO_ON_ROOT( "   + num cycles:                                   " << multigridSettings.numCycles )
   WALBERLA_LOG_INFO_ON_ROOT( "   + FMG inner iterations:                         " << multigridSettings.fmgInnerIterations )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + absolute residual tolerance:                  " << multigridSettings.absoluteResidualTolerance )
   WALBERLA_LOG_INFO_ON_ROOT( " - smoother settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + num GS velocity:                              " << smootherSettings.numGSVelocity )
   WALBERLA_LOG_INFO_ON_ROOT( "   + symmetric velocity:                           " << smootherSettings.symmGSVelocity )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + omega estimation:                             " << ( smootherSettings.estimateOmega ? "yes" : "no" ) )
   if ( smootherSettings.estimateOmega )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + omega estimation level:                       " << smootherSettings.omegaEstimationLevel )
      WALBERLA_LOG_INFO_ON_ROOT(
          "   + omega estimation iterations:                  " << smootherSettings.omegaEstimationIterations )
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + fixed omega:                                  " << smootherSettings.omega )
   }
   WALBERLA_LOG_INFO_ON_ROOT( " - coarse grid settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + coarse grid solver:                           " << coarseGridSettings.solverType )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + absolute residual tolerance:                  " << coarseGridSettings.absoluteResidualTolerance )
   WALBERLA_LOG_INFO_ON_ROOT( "   + max iterations:                               " << coarseGridSettings.maxIterations )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK:                                          " << ( vtk ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   const auto storageGlobalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageGlobalInfo );

   uint_t iteration = 0;

   WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ..." )

   StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   StokesFunction< real_t > f =
       RHSisZero ? StokesFunction< real_t >( "f", storage, minLevel, maxLevel > minLevel ? maxLevel - 1 : maxLevel ) :
                   StokesFunction< real_t >( "f", storage, minLevel, maxLevel );
   StokesFunction< real_t > r =
       RHSisZero ? StokesFunction< real_t >( "r", storage, 0, 0 ) : StokesFunction< real_t >( "r", storage, minLevel, maxLevel );
   StokesFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Allocating and assembling operators ..." )

   StokesOperator       A( storage, minLevel, maxLevel );
   VelocityMassOperator velocityMassOperator( storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating solution ..." )

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      u.uvw().interpolate( { initialU, initialV, initialW }, level, All );
      u.uvw().interpolate( { solutionU, solutionV, solutionW }, level, DirichletBoundary );
      u.p().interpolate( initialP, level, All );

      if ( RHSisZero )
      {
         if ( level < maxLevel )
         {
            f.uvw().interpolate( 0, level, All );
            f.p().interpolate( 0, level, All );
         }
      }
      else
      {
         tmp.uvw().interpolate( { rhsU,rhsV, rhsW }, level, All );

         velocityMassOperator.apply( tmp.uvw()[0], f.uvw()[0], level, All );
         velocityMassOperator.apply( tmp.uvw()[1], f.uvw()[1], level, All );
         velocityMassOperator.apply( tmp.uvw()[2], f.uvw()[2], level, All );

         f.p().interpolate( 0, level, All );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up VTK ..." )

   VTKOutput vtkOutput( "vtk", "TME", storage );
   vtkOutput.add( u );
   vtkOutput.add( tmp );
   if ( !RHSisZero )
   {
      vtkOutput.add( f );
      vtkOutput.add( r );
   }

   if ( vtk )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   real_t residualL2Velocity;
   real_t residualL2Pressure;
   real_t errorL2Velocity;
   real_t errorL2Pressure;

   WALBERLA_LOG_INFO_ON_ROOT( "Preparing solvers ..." )

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >( projectPressurefterRestriction );

   std::shared_ptr< Solver< typename StokesOperator::VelocityOperator_T > > forwardVelocitySmoother =
       std::make_shared< SORSmoother< typename StokesOperator::VelocityOperator_T > >( 1.0 );
   std::shared_ptr< Solver< typename StokesOperator::VelocityOperator_T > > symmetricVelocitySmoother =
       std::make_shared< SymmetricSORSmoother< typename StokesOperator::VelocityOperator_T > >( 1.0 );

   auto uzawaForwardVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, forwardVelocitySmoother );
   auto uzawaSymmetricVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, symmetricVelocitySmoother );

   auto uzawaVelocityPreconditioner =
       ( smootherSettings.symmGSVelocity ? uzawaSymmetricVelocityPreconditioner : uzawaForwardVelocityPreconditioner );

   std::vector< uint_t > rhsZeroLevels = { maxLevel };
   auto                  smoother      = std::make_shared< UzawaSmoother< StokesOperator > >( storage,
                                                                        uzawaVelocityPreconditioner,
                                                                        tmp,
                                                                        minLevel,
                                                                        maxLevel,
                                                                        smootherSettings.omega,
                                                                        Inner | NeumannBoundary,
                                                                        smootherSettings.numGSVelocity,
                                                                        false,
                                                                        1,
                                                                        RHSisZero,
                                                                        rhsZeroLevels );

   if ( smootherSettings.estimateOmega )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "Estimating omega (" << smootherSettings.omegaEstimationIterations
                                                      << " power iterations on level " << smootherSettings.omegaEstimationLevel
                                                      << ") ..." );
      smootherSettings.omega = estimateUzawaRelaxationParameter( storage,
                                                                 uzawaSymmetricVelocityPreconditioner,
                                                                 smootherSettings.omegaEstimationLevel,
                                                                 smootherSettings.omegaEstimationIterations,
                                                                 smootherSettings.numGSVelocity );
      if ( std::isnan( smootherSettings.omega ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Omega is NaN, setting to: " << smootherSettings.omega );
         smootherSettings.omega = 0;
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Setting omega to estimate: " << smootherSettings.omega );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      smoother->setRelaxationParameter( smootherSettings.omega );
   }

   // coarse grid solver type:
   // 0: MUMPS                          (PETSc)
   // 1: block preconditioned MINRES    (PETSc)

   std::shared_ptr< Solver< StokesOperator > > coarseGridSolverInternal;

   if ( coarseGridSettings.solverType == 0 )
   {
      auto petscSolverInternalTmp = std::make_shared< PETScLUSolver< StokesOperator > >( storage, minLevel );

      coarseGridSolverInternal = petscSolverInternalTmp;
   }
   else if ( coarseGridSettings.solverType == 1 )
   {
      auto petscSolverInternalTmp = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
          storage, minLevel, coarseGridSettings.absoluteResidualTolerance, coarseGridSettings.maxIterations, 1 );
      petscSolverInternalTmp->setVerbose( true );
      coarseGridSolverInternal = petscSolverInternalTmp;
   }

   auto coarseGridSolver = std::make_shared< TimedSolver< StokesOperator > >( coarseGridSolverInternal );

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                          tmp,
                                                                                          smoother,
                                                                                          coarseGridSolver,
                                                                                          restrictionOperator,
                                                                                          prolongationOperator,
                                                                                          minLevel,
                                                                                          maxLevel,
                                                                                          multigridSettings.preSmooth,
                                                                                          multigridSettings.postSmooth,
                                                                                          multigridSettings.incSmooth,
                                                                                          CycleType::VCYCLE,
                                                                                          RHSisZero,
                                                                                          0 );

   auto fmgProlongation = std::make_shared< FMGProlongation >();

   auto postCycle = [&]( uint_t currentLevel ) {
      smoother->setRHSZeroLevels( { currentLevel + 1 } );

      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), currentLevel );
      }

      errorAndResidual( A,
                        u,
                        f,
                        r,
                        tmp,
                        solutionU,
                        solutionV,
                        solutionW,
                        solutionP,
                        currentLevel,
                        errorFlag,
                        RHSisZero,
                        residualL2Velocity,
                        residualL2Pressure,
                        errorL2Velocity,
                        errorL2Pressure );

      writeDataRow( iteration, "F", currentLevel, errorL2Velocity, residualL2Velocity, errorL2Pressure, residualL2Pressure, db );
      iteration++;
   };

   FullMultigridSolver< StokesOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, multigridSettings.fmgInnerIterations, postCycle );

   WALBERLA_LOG_INFO_ON_ROOT( "Calculating initial error and residual ..." )

   errorAndResidual( A,
                     u,
                     f,
                     r,
                     tmp,
                     solutionU,
                     solutionV,
                     solutionW,
                     solutionP,
                     maxLevel,
                     errorFlag,
                     RHSisZero,
                     residualL2Velocity,
                     residualL2Pressure,
                     errorL2Velocity,
                     errorL2Pressure );

   WALBERLA_LOG_INFO_ON_ROOT( "Obtaining function allocation info ..." )

   printFunctionAllocationInfo( *storage, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "Gathering memory usage info ..." )

   //double sumGBAllocatedRUsage, minGBAllocatedRUsage, maxGBAllocatedRUsage;
   //double sumGBAllocatedPetsc, minGBAllocatedPetsc, maxGBAllocatedPetsc;

   auto [sumGBAllocatedRUsage, minGBAllocatedRUsage, maxGBAllocatedRUsage] =
       printAndGetCurrentMemoryUsage( MemoryUsageDeterminationType::C_RUSAGE );
   auto [sumGBAllocatedPetsc, minGBAllocatedPetsc, maxGBAllocatedPetsc] =
       printAndGetCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );

   writeDataHeader();

   db.setConstantEntry( "discretization", discretization );

   db.setConstantEntry( "min_level", minLevel );
   db.setConstantEntry( "max_level", maxLevel );
   db.setConstantEntry( "unknowns", unknowns );
   db.setConstantEntry( "h_min", hMin );
   db.setConstantEntry( "h_max", hMax );
   db.setConstantEntry( "num_macro_cells", storage->getNumberOfGlobalCells() );
   db.setConstantEntry( "num_macro_faces", storage->getNumberOfGlobalFaces() );
   db.setConstantEntry( "num_macro_edges", storage->getNumberOfGlobalEdges() );
   db.setConstantEntry( "num_macro_vertices", storage->getNumberOfGlobalVertices() );
   db.setConstantEntry( "num_macro_primitives", storage->getNumberOfGlobalPrimitives() );

   db.setConstantEntry( "num_cycles", multigridSettings.numCycles );
   db.setConstantEntry( "pre_smooth", multigridSettings.preSmooth );
   db.setConstantEntry( "post_smooth", multigridSettings.postSmooth );
   db.setConstantEntry( "inc_smooth", multigridSettings.incSmooth );
   db.setConstantEntry( "fmg_inner_iterations", multigridSettings.fmgInnerIterations );

   db.setConstantEntry( "omega", smootherSettings.omega );
   db.setConstantEntry( "omega_estimation_iterations", smootherSettings.omegaEstimationIterations );
   db.setConstantEntry( "omega_estimation_level", smootherSettings.omegaEstimationLevel );
   db.setConstantEntry( "estimate_omega", smootherSettings.estimateOmega );
   db.setConstantEntry( "symm_gs_velocity", smootherSettings.symmGSVelocity );
   db.setConstantEntry( "num_gs_velocity", smootherSettings.numGSVelocity );

   db.setConstantEntry( "coarse_grid_solver_type", coarseGridSettings.solverType );
   db.setConstantEntry( "coarse_grid_max_iterations", coarseGridSettings.maxIterations );
   db.setConstantEntry( "coarse_grid_absolute_residual_tolerance", coarseGridSettings.absoluteResidualTolerance );

   db.setConstantEntry( "allocated_mem_rusage_gb_sum", sumGBAllocatedRUsage );
   db.setConstantEntry( "allocated_mem_rusage_gb_min", minGBAllocatedRUsage );
   db.setConstantEntry( "allocated_mem_rusage_gb_max", maxGBAllocatedRUsage );

   db.setConstantEntry( "allocated_mem_petsc_gb_sum", sumGBAllocatedPetsc );
   db.setConstantEntry( "allocated_mem_petsc_gb_min", minGBAllocatedPetsc );
   db.setConstantEntry( "allocated_mem_petsc_gb_max", maxGBAllocatedPetsc );

   writeDataRow( iteration, "I", 0, errorL2Velocity, residualL2Velocity, errorL2Pressure, residualL2Pressure, db );
   iteration++;

   timer->stop( "Setup" );
   timer->start( "Solve" );

   if ( multigridSettings.fmgInnerIterations > 0 )
   {
      smoother->setRHSZeroLevels( { minLevel } );
      fullMultigridSolver.solve( A, u, f, maxLevel );
   }
   else
   {
      multigridSolver->solve( A, u, f, maxLevel );

      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), maxLevel );
      }

      errorAndResidual( A,
                        u,
                        f,
                        r,
                        tmp,
                        solutionU,
                        solutionV,
                        solutionW,
                        solutionP,
                        maxLevel,
                        errorFlag,
                        RHSisZero,
                        residualL2Velocity,
                        residualL2Pressure,
                        errorL2Velocity,
                        errorL2Pressure );

      writeDataRow( iteration, "V", 0, errorL2Velocity, residualL2Velocity, errorL2Pressure, residualL2Pressure, db );
      iteration++;
   }

   smoother->setRHSZeroLevels( { maxLevel } );

   for ( uint_t i = 1; i < multigridSettings.numCycles && ( multigridSettings.absoluteResidualTolerance < residualL2Velocity ||
                                                            multigridSettings.absoluteResidualTolerance < residualL2Pressure );
         i++ )
   {
      multigridSolver->solve( A, u, f, maxLevel );

      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), maxLevel );
      }

      errorAndResidual( A,
                        u,
                        f,
                        r,
                        tmp,
                        solutionU,
                        solutionV,
                        solutionW,
                        solutionP,
                        maxLevel,
                        errorFlag,
                        RHSisZero,
                        residualL2Velocity,
                        residualL2Pressure,
                        errorL2Velocity,
                        errorL2Pressure );

      writeDataRow( iteration, "V", 0, errorL2Velocity, residualL2Velocity, errorL2Pressure, residualL2Pressure, db );
      iteration++;
   }

   timer->stop( "Solve" );

   if ( vtk )
   {
      vtkOutput.write( maxLevel, 1 );
   }

   timer->stop( "Total" );
   writeTimingTreeJSON( *timer, timingFile );
}

void solve( const std::shared_ptr< PrimitiveStorage >&              storage,
            Discretization                                          discretization,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionU,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionV,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionW,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionP,
            const std::function< real_t( const hyteg::Point3D& ) >& initialU,
            const std::function< real_t( const hyteg::Point3D& ) >& initialV,
            const std::function< real_t( const hyteg::Point3D& ) >& initialW,
            const std::function< real_t( const hyteg::Point3D& ) >& initialP,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsU,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsV,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsW,
            uint_t                                                  minLevel,
            uint_t                                                  maxLevel,
            MultigridSettings                                       multigridSettings,
            SmootherSettings                                        smootherSettings,
            CoarseGridSettings                                      coarseGridSettings,
            bool                                                    projectPressure,
            bool                                                    projectPressurefterRestriction,
            bool                                                    vtk,
            std::string                                             benchmarkName,
            FixedSizeSQLDB                                          db,
            std::string                                             timingFile,
            bool                                                    RHSisZero )
{
   if ( discretization == Discretization::P2_P1 )
   {
      solveRHS0Implementation< P2P1TaylorHoodFunction,
                               P2P1TaylorHoodStokesOperator,
                               P2ConstantMassOperator,
                               P1ConstantMassOperator,
                               P2P1StokesToP2P1StokesRestriction,
                               P2P1StokesToP2P1StokesProlongation,
                               P2P1StokesToP2P1StokesProlongation >( storage,
                                                                     solutionU,
                                                                     solutionV,
                                                                     solutionW,
                                                                     solutionP,
                                                                     initialU,
                                                                     initialV,
                                                                     initialW,
                                                                     initialP,
                                                                     rhsU,
                                                                     rhsV,
                                                                     rhsW,
                                                                     minLevel,
                                                                     maxLevel,
                                                                     multigridSettings,
                                                                     smootherSettings,
                                                                     coarseGridSettings,
                                                                     projectPressure,
                                                                     projectPressurefterRestriction,
                                                                     vtk,
                                                                     benchmarkName,
                                                                     db,
                                                                     timingFile,
                                                                     RHSisZero );
   }
   else
   {
      solveRHS0Implementation< P1StokesFunction,
                               P1P1StokesOperator,
                               P1ConstantMassOperator,
                               P1ConstantMassOperator,
                               P1P1StokesToP1P1StokesRestriction,
                               P1P1StokesToP1P1StokesProlongation,
                               P1P1StokesToP1P1StokesProlongation >( storage,
                                                                     solutionU,
                                                                     solutionV,
                                                                     solutionW,
                                                                     solutionP,
                                                                     initialU,
                                                                     initialV,
                                                                     initialW,
                                                                     initialP,
                                                                     rhsU,
                                                                     rhsV,
                                                                     rhsW,
                                                                     minLevel,
                                                                     maxLevel,
                                                                     multigridSettings,
                                                                     smootherSettings,
                                                                     coarseGridSettings,
                                                                     projectPressure,
                                                                     projectPressurefterRestriction,
                                                                     vtk,
                                                                     benchmarkName,
                                                                     db,
                                                                     timingFile,
                                                                     RHSisZero );
   }
}

void domainInfo( const std::shared_ptr< PrimitiveStorage >& storage,
                 Discretization                             discretization,
                 uint_t                                     minLevel,
                 uint_t                                     maxLevel,
                 bool                                       vtk,
                 std::string                                benchmarkName )
{
   printGitInfo();
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, "vtk/", benchmarkName + "_domain" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Unknowns (incl. boundary) and total allocated floating point values per Stokes function:" )
   WALBERLA_LOG_INFO_ON_ROOT( "" )
   unsigned long long unknownsSum        = 0;
   unsigned long long allocatedMemorySum = 0;
   WALBERLA_LOG_INFO_ON_ROOT( " level |       unknowns |          total " )
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      unsigned long long unknowns             = 0;
      unsigned long long allocatedMemoryLevel = 0;
      if ( discretization == Discretization::P2_P1 )
      {
         unknowns             = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
         allocatedMemoryLevel = p2p1globalFunctionMemorySize( level, storage );
      }
      else
      {
         unknowns             = numberOfGlobalDoFs< P1StokesFunctionTag >( *storage, level );
         allocatedMemoryLevel = p1p1globalFunctionMemorySize( level, storage );
      }
      unknownsSum += unknowns;
      allocatedMemorySum += allocatedMemoryLevel;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %5d | %14llu | %14llu ", level, unknowns, allocatedMemoryLevel ) )
   }
   WALBERLA_LOG_INFO_ON_ROOT( " ------+----------------+--------------- " )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "   sum | %14llu | %14llu ", unknownsSum, allocatedMemorySum ) )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   const real_t hMin                 = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t hMax                 = MeshQuality::getMaximalEdgeLength( storage, maxLevel );
   const auto   discretizationString = ( discretization == Discretization::P2_P1 ? "P2-P1" : "P1-P1" );

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << benchmarkName )
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + elements:                                     " << discretizationString );
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + min level:                                    " << minLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + max level:                                    " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK:                                          " << ( vtk ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   const auto storageGlobalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageGlobalInfo );

   WALBERLA_LOG_INFO_ON_ROOT( "" )
}

} // namespace scaling_workshop
} // namespace hyteg
