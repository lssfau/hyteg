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
namespace tme_benchmarks {

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

void writeDataHeader()
{
   WALBERLA_LOG_INFO_ON_ROOT(
       " iteration | iteration type | FMG level | L2 error velocity | L2 alg err. rel. vel. | L2 err. grid inc. vel. | L2 residual velocity | L2 error pressure | L2 alg err. rel. pre. | L2 err. grid inc. pre. | L2 residual pressure " );
   WALBERLA_LOG_INFO_ON_ROOT(
       "-----------+----------------+-----------+-------------------+-----------------------+----------------------+-------------------+-----------------------+----------------------" );
}

// iterationType: I: initial
//                V: v-cycle
//                F: FMG
void writeDataRow( uint_t          iteration,
                   std::string     iterationType,
                   uint_t          fmgLevel,
                   real_t          errorL2Velocity,
                   real_t          errorRelL2Velocity,
                   real_t          residualL2Velocity,
                   real_t          errorL2Pressure,
                   real_t          errorRelL2Pressure,
                   real_t          residualL2Pressure,
                   real_t          errorL2GridIncrementVelocity,
                   real_t          errorL2GridIncrementPressure,
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
   db.setVariableEntry( "error_rel_l2_velocity", errorRelL2Velocity );
   db.setVariableEntry( "residual_l2_velocity", residualL2Velocity );
   db.setVariableEntry( "error_l2_pressure", errorL2Pressure );
   db.setVariableEntry( "error_rel_l2_pressure", errorRelL2Pressure );
   db.setVariableEntry( "residual_l2_pressure", residualL2Pressure );
   db.setVariableEntry( "error_l2_grid_increment_velocity", errorL2GridIncrementVelocity );
   db.setVariableEntry( "error_l2_grid_increment_pressure", errorL2GridIncrementPressure );

   db.writeRowOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( " %9d | %14s | %9s | %17.5e | %21.5e | %22.5e | %20.5e | %17.5e | %21.5e | %22.5e | %20.5e ",
                         iteration,
                         iterationType.c_str(),
                         fmgLevelString.c_str(),
                         errorL2Velocity,
                         errorRelL2Velocity,
                         errorL2GridIncrementVelocity,
                         residualL2Velocity,
                         errorL2Pressure,
                         errorRelL2Pressure,
                         errorL2GridIncrementPressure,
                         residualL2Pressure ) );
}

template < template < typename > class StokesFunction,
           typename StokesOperator,
           typename VelocityMassOperator,
           typename PressureMassOperator,
           typename Restriction,
           typename Prolongation,
           typename FMGProlongation >
void solveImplementation( const std::shared_ptr< PrimitiveStorage >&              storage,
                          StokesFunction< real_t >&                               uSolution,
                          const std::function< real_t( const hyteg::Point3D& ) >& solutionU,
                          const std::function< real_t( const hyteg::Point3D& ) >& solutionV,
                          const std::function< real_t( const hyteg::Point3D& ) >& solutionW,
                          const std::function< real_t( const hyteg::Point3D& ) >& solutionP,
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
                          bool                                                    calculateDiscretizationError,
                          bool                                                    discretizationErrorWasCalculated,
                          uint_t                                                  errorCalculationLevelIncrement,
                          bool                                                    solveWithCoarseGridSolverOnEachFMGLevel,
                          bool                                                    vtk,
                          const std::string&                                      benchmarkName,
                          bool                                                    verbose,
                          std::string                                             dbFile )
{
   printGitInfo();

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, "vtk/", benchmarkName + "_domain" );
   }

   auto timer = storage->getTimingTree();
   timer->start( "Total" );
   timer->start( "Setup" );

   const uint_t unknowns   = numberOfGlobalDoFs< typename StokesFunction< real_t >::Tag >( *storage, maxLevel );
   const real_t hMin       = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t hMax       = MeshQuality::getMaximalEdgeLength( storage, maxLevel );
   const auto   errorLevel = maxLevel + errorCalculationLevelIncrement;
   const auto   discretization =
       ( std::is_same< typename StokesFunction< real_t >::Tag, P2P1TaylorHoodFunctionTag >::value ? "P2-P1" : "P1-P1" );

   if ( calculateDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### CALCULATING DISCRETIZATION ERROR ###" );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << benchmarkName )
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
   WALBERLA_LOG_INFO_ON_ROOT( "   + error interpolation level:                    " << errorLevel )
   WALBERLA_LOG_INFO_ON_ROOT(
       "   + solve on each FMG level:                      " << ( solveWithCoarseGridSolverOnEachFMGLevel ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK:                                          " << ( vtk ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + verbose:                                      " << verbose )
   WALBERLA_LOG_INFO_ON_ROOT( "   + database file:                                " << dbFile )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   std::map< uint_t, uint_t > unknownsLevel;
   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      unknownsLevel[l] = numberOfGlobalDoFs< typename StokesFunction< real_t >::Tag >( *storage, l );
      WALBERLA_LOG_INFO_ON_ROOT( "unknowns on level " << l << ": " << unknownsLevel[l] );
   }

   FixedSizeSQLDB db( dbFile );

   uint_t iteration = 0;

   const auto errorFlag = Inner | NeumannBoundary | FreeslipBoundary;

   StokesFunction< real_t > u( "u", storage, minLevel, errorLevel );
   StokesFunction< real_t > f( "f", storage, minLevel, errorLevel );
   StokesFunction< real_t > r( "r", storage, minLevel, errorLevel );
   StokesFunction< real_t > err( "error", storage, minLevel, errorLevel );
   StokesFunction< real_t > errDiscr( "error_discr", storage, minLevel, errorLevel );
   StokesFunction< real_t > tmp( "tmp", storage, minLevel, errorLevel );
   StokesFunction< real_t > tmpFMG( "tmpFMG", storage, minLevel, errorLevel );
   StokesFunction< real_t > exact( "exact", storage, minLevel, errorLevel );
   // compute on level l the quantity
   //   || u_l - I_{l-1}^l u_{l-1} ||_2
   // which should converge as the discretization error
   StokesFunction< real_t > uFMGSolution( "u_fmg_solution", storage, minLevel, errorLevel );
   StokesFunction< real_t > errGridIncrement( "error_grid_increment", storage, minLevel, errorLevel );

   StokesOperator       A( storage, minLevel, errorLevel );
   VelocityMassOperator MVelocity( storage, minLevel, errorLevel );
   PressureMassOperator MPressure( storage, minLevel, errorLevel );

   for ( uint_t level = minLevel; level <= errorLevel; level++ )
   {
      exact.uvw().interpolate( { solutionU, solutionV, solutionW },level, All );
      exact.p().interpolate( solutionP, level, All );

      if ( projectPressure )
      {
         vertexdof::projectMean( exact.p(), level );
      }

      u.uvw().interpolate( { solutionU, solutionV, solutionW },level, DirichletBoundary );
      tmp.uvw().interpolate( { rhsU, rhsV, rhsW }, level, All );

      MVelocity.apply( tmp.uvw()[0], f.uvw()[0], level, All );
      MVelocity.apply( tmp.uvw()[1], f.uvw()[1], level, All );
      MVelocity.apply( tmp.uvw()[2], f.uvw()[2], level, All );
   }

   const auto solutionNormVelocity = normL2Velocity( uSolution, MVelocity, tmp, maxLevel, errorFlag );
   const auto solutionNormPressure = normL2Scalar( uSolution.p(), MPressure, tmp.p(), maxLevel, errorFlag );

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
   auto smoother = std::make_shared< UzawaSmoother< StokesOperator > >( storage,
                                                                        uzawaVelocityPreconditioner,
                                                                        minLevel,
                                                                        maxLevel,
                                                                        smootherSettings.omega,
                                                                        Inner | NeumannBoundary,
                                                                        smootherSettings.numGSVelocity );

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

   // iterative or direct solver to solve on each level - obtaining discretization accuracy
   // used to indicate incremental error reduction rate
   std::vector< std::shared_ptr< Solver< StokesOperator > > > fmgLevelWiseSolver;
   if ( solveWithCoarseGridSolverOnEachFMGLevel )
   {
      for ( uint_t l = 0; l <= maxLevel; l++ )
      {
         if ( coarseGridSettings.solverType == 0 )
         {
            auto slvr = std::make_shared< PETScLUSolver< StokesOperator > >( storage, l );
            fmgLevelWiseSolver.push_back( slvr );
         }
         else
         {
            auto slvr = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
                storage, l, coarseGridSettings.absoluteResidualTolerance, coarseGridSettings.maxIterations, 1 );
            slvr->setVerbose( true );
            fmgLevelWiseSolver.push_back( slvr );
         }
      }
   }

   auto coarseGridSolver = std::make_shared< TimedSolver< StokesOperator > >( coarseGridSolverInternal );

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                          smoother,
                                                                                          coarseGridSolver,
                                                                                          restrictionOperator,
                                                                                          prolongationOperator,
                                                                                          minLevel,
                                                                                          maxLevel,
                                                                                          multigridSettings.preSmooth,
                                                                                          multigridSettings.postSmooth,
                                                                                          multigridSettings.incSmooth,
                                                                                          CycleType::VCYCLE );

   auto fmgProlongation = std::make_shared< FMGProlongation >();

   auto postCycle = [&]( uint_t currentLevel ) {

      if (solveWithCoarseGridSolverOnEachFMGLevel)
      {
         fmgLevelWiseSolver[currentLevel]->solve( A, u, f, currentLevel );
      }


      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), currentLevel );
      }

      for ( uint_t l = minLevel; l <= errorLevel; l++ )
      {
         tmpFMG.interpolate( 0, l, All );
      }

      tmpFMG.assign( { 1.0 }, { u }, currentLevel, All );
      uFMGSolution.assign( { 1.0 }, { u }, currentLevel, All );

      const auto fmgErrorLevel = currentLevel + errorCalculationLevelIncrement;
      for ( uint_t prolongationSourceLevel = currentLevel; prolongationSourceLevel < fmgErrorLevel; prolongationSourceLevel++ )
      {
         prolongationOperator->prolongate( tmpFMG, prolongationSourceLevel, All );
      }

      residual( tmpFMG, f, A, tmp, currentLevel, errorFlag, r );
      error( tmpFMG, exact, fmgErrorLevel, errorFlag, err );
      if ( currentLevel == maxLevel )
      {
         error( u, uSolution, maxLevel, errorFlag, errDiscr );
      }

      // grid increment error
      if ( currentLevel > 0 )
      {
         errGridIncrement.assign( { 1.0 }, { uFMGSolution }, currentLevel - 1, All );
         prolongationOperator->prolongate( errGridIncrement, currentLevel - 1, All );
         errGridIncrement.assign( { 1.0, -1.0 }, { uFMGSolution, errGridIncrement }, currentLevel, All );
      }

      auto errorL2VelocityFMG    = normL2Velocity( err, MVelocity, tmp, fmgErrorLevel, errorFlag );
      auto residualL2VelocityFMG = normL2Velocity( r, MVelocity, tmp, currentLevel, errorFlag );

      auto errorL2PressureFMG    = normL2Scalar( err.p(), MPressure, tmp.p(), fmgErrorLevel, errorFlag );
      auto residualL2PressureFMG = normL2Scalar( r.p(), MPressure, tmp.p(), currentLevel, errorFlag );

      auto errorRelL2VelocityFMG = 0.;
      auto errorRelL2PressureFMG = 0.;
      if ( !calculateDiscretizationError && discretizationErrorWasCalculated && currentLevel == maxLevel )
      {
         errorRelL2VelocityFMG = normL2Velocity( errDiscr, MVelocity, tmp, maxLevel, errorFlag ) / solutionNormVelocity;
         errorRelL2PressureFMG = normL2Scalar( errDiscr.p(), MPressure, tmp.p(), maxLevel, errorFlag ) / solutionNormVelocity;
      }

      auto errorL2GridIncrementVelocity = normL2Velocity( errGridIncrement, MVelocity, tmp, currentLevel, errorFlag );
      auto errorL2GridIncrementPressure = normL2Scalar( errGridIncrement.p(), MPressure, tmp.p(), currentLevel, errorFlag );

      writeDataRow( iteration,
                    "F",
                    currentLevel,
                    errorL2VelocityFMG,
                    errorRelL2VelocityFMG,
                    residualL2VelocityFMG,
                    errorL2PressureFMG,
                    errorRelL2PressureFMG,
                    residualL2PressureFMG,
                    errorL2GridIncrementVelocity,
                    errorL2GridIncrementPressure,
                    db );
      iteration++;
   };

   FullMultigridSolver< StokesOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, multigridSettings.fmgInnerIterations, postCycle );

   for ( uint_t prolongationSourceLevel = maxLevel; prolongationSourceLevel < errorLevel; prolongationSourceLevel++ )
   {
      prolongationOperator->prolongate( u, prolongationSourceLevel, All );
   }

   residual( u, f, A, tmp, maxLevel, errorFlag, r );
   error( u, exact, errorLevel, errorFlag, err );
   error( u, uSolution, maxLevel, errorFlag, errDiscr );

   auto errorL2Velocity    = normL2Velocity( err, MVelocity, tmp, errorLevel, errorFlag );
   auto residualL2Velocity = normL2Velocity( r, MVelocity, tmp, maxLevel, errorFlag );
   auto errorL2Pressure    = normL2Scalar( err.p(), MPressure, tmp.p(), errorLevel, errorFlag );
   auto residualL2Pressure = normL2Scalar( r.p(), MPressure, tmp.p(), maxLevel, errorFlag );

   auto errorAlgL2Velocity = 0.;
   auto errorAlgL2Pressure = 0.;
   if ( !calculateDiscretizationError && discretizationErrorWasCalculated )
   {
      errorAlgL2Velocity = normL2Velocity( errDiscr, MVelocity, tmp, maxLevel, errorFlag ) / solutionNormVelocity;
      errorAlgL2Pressure = normL2Scalar( errDiscr.p(), MPressure, tmp.p(), maxLevel, errorFlag ) / solutionNormPressure;
   }

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

   db.setConstantEntry( "error_calculation_level_increment", errorCalculationLevelIncrement );
   db.setConstantEntry( "solve_with_coarse_grid_solver_on_each_fmg_level", solveWithCoarseGridSolverOnEachFMGLevel );

   writeDataRow( iteration,
                 "I",
                 0,
                 errorL2Velocity,
                 errorAlgL2Velocity,
                 residualL2Velocity,
                 errorL2Pressure,
                 errorAlgL2Pressure,
                 residualL2Pressure,
                 0,
                 0,
                 db );
   iteration++;

   VTKOutput vtkOutput( "vtk", "TME", storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( err );
   vtkOutput.add( exact );

   if ( vtk )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   if ( multigridSettings.fmgInnerIterations > 0 )
   {
      fullMultigridSolver.solve( A, u, f, maxLevel );
   }
   else
   {
      multigridSolver->solve( A, u, f, maxLevel );

      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), maxLevel );
      }

      for ( uint_t prolongationSourceLevel = maxLevel; prolongationSourceLevel < errorLevel; prolongationSourceLevel++ )
      {
         prolongationOperator->prolongate( u, prolongationSourceLevel, All );
      }

      residual( u, f, A, tmp, maxLevel, errorFlag, r );
      error( u, exact, errorLevel, errorFlag, err );
      error( u, uSolution, maxLevel, errorFlag, errDiscr );

      errorL2Velocity    = normL2Velocity( err, MVelocity, tmp, errorLevel, errorFlag );
      residualL2Velocity = normL2Velocity( r, MVelocity, tmp, maxLevel, errorFlag );
      errorL2Pressure    = normL2Scalar( err.p(), MPressure, tmp.p(), errorLevel, errorFlag );
      residualL2Pressure = normL2Scalar( r.p(), MPressure, tmp.p(), maxLevel, errorFlag );

      errorAlgL2Velocity = 0.;
      errorAlgL2Pressure = 0.;
      if ( !calculateDiscretizationError && discretizationErrorWasCalculated )
      {
         errorAlgL2Velocity = normL2Velocity( errDiscr, MVelocity, tmp, maxLevel, errorFlag ) / solutionNormVelocity;
         errorAlgL2Pressure = normL2Scalar( errDiscr.p(), MPressure, tmp.p(), maxLevel, errorFlag ) / solutionNormPressure;
      }

      writeDataRow( iteration,
                    "V",
                    0,
                    errorL2Velocity,
                    errorAlgL2Velocity,
                    residualL2Velocity,
                    errorL2Pressure,
                    errorAlgL2Pressure,
                    residualL2Pressure,
                    0,
                    0,
                    db );
      iteration++;
   }

   for ( uint_t i = 1; i < multigridSettings.numCycles && ( multigridSettings.absoluteResidualTolerance < residualL2Velocity ||
                                                            multigridSettings.absoluteResidualTolerance < residualL2Pressure );
         i++ )
   {
      multigridSolver->solve( A, u, f, maxLevel );

      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), maxLevel );
      }

      for ( uint_t prolongationSourceLevel = maxLevel; prolongationSourceLevel < errorLevel; prolongationSourceLevel++ )
      {
         prolongationOperator->prolongate( u, prolongationSourceLevel, All );
      }

      residual( u, f, A, tmp, maxLevel, errorFlag, r );
      error( u, exact, errorLevel, errorFlag, err );
      error( u, uSolution, maxLevel, errorFlag, errDiscr );

      errorL2Velocity    = normL2Velocity( err, MVelocity, tmp, errorLevel, errorFlag );
      residualL2Velocity = normL2Velocity( r, MVelocity, tmp, maxLevel, errorFlag );
      errorL2Pressure    = normL2Scalar( err.p(), MPressure, tmp.p(), errorLevel, errorFlag );
      residualL2Pressure = normL2Scalar( r.p(), MPressure, tmp.p(), maxLevel, errorFlag );

      errorAlgL2Velocity = 0.;
      errorAlgL2Pressure = 0.;
      if ( !calculateDiscretizationError && discretizationErrorWasCalculated )
      {
         errorAlgL2Velocity = normL2Velocity( errDiscr, MVelocity, tmp, maxLevel, errorFlag ) / solutionNormVelocity;
         errorAlgL2Pressure = normL2Scalar( errDiscr.p(), MPressure, tmp.p(), maxLevel, errorFlag ) / solutionNormPressure;
      }

      writeDataRow( iteration,
                    "V",
                    0,
                    errorL2Velocity,
                    errorAlgL2Velocity,
                    residualL2Velocity,
                    errorL2Pressure,
                    errorAlgL2Pressure,
                    residualL2Pressure,
                    0,
                    0,
                    db );
      iteration++;
   }

   if ( vtk )
   {
      vtkOutput.write( maxLevel, 1 );
   }

   if ( calculateDiscretizationError )
   {
      for ( uint_t l = minLevel; l <= maxLevel; l++ )
      {
         uSolution.assign( { 1.0 }, { u }, maxLevel, All );
      }
   }
}

void solve( const std::shared_ptr< PrimitiveStorage >& storage,
            Discretization                             discretization,
            bool,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionU,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionV,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionW,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionP,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsU,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsV,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsW,
            uint_t                                                  minLevel,
            uint_t                                                  maxLevel,
            MultigridSettings                                       multigridSettings,
            SmootherSettings                                        smootherSettings,
            CoarseGridSettings                                      coarseGridSettings,
            MultigridSettings                                       multigridSettingsDiscrError,
            SmootherSettings                                        smootherSettingsDiscrError,
            CoarseGridSettings                                      coarseGridSettingsDiscrError,
            bool                                                    projectPressure,
            bool                                                    projectPressurefterRestriction,
            bool                                                    calculateDiscretizationError,
            uint_t                                                  normCalculationLevelIncrement,
            bool                                                    solveWithCoarseGridSolverOnEachFMGLevel,
            bool                                                    vtk,
            const std::string&                                      benchmarkName,
            bool                                                    verbose,
            std::string                                             dbFile,
            std::string                                             dbFileDiscr )
{
   if ( discretization == Discretization::P2_P1 )
   {
      P2P1TaylorHoodFunction< real_t > uSolution( "uSolution", storage, minLevel, maxLevel + normCalculationLevelIncrement );

      if ( calculateDiscretizationError )
      {
         solveImplementation< P2P1TaylorHoodFunction,
                              P2P1TaylorHoodStokesOperator,
                              P2ConstantMassOperator,
                              P1ConstantMassOperator,
                              P2P1StokesToP2P1StokesRestriction,
                              P2P1StokesToP2P1StokesProlongation,
                              P2P1StokesToP2P1StokesProlongation >( storage,
                                                                    uSolution,
                                                                    solutionU,
                                                                    solutionV,
                                                                    solutionW,
                                                                    solutionP,
                                                                    rhsU,
                                                                    rhsV,
                                                                    rhsW,
                                                                    minLevel,
                                                                    maxLevel,
                                                                    multigridSettingsDiscrError,
                                                                    smootherSettingsDiscrError,
                                                                    coarseGridSettingsDiscrError,
                                                                    projectPressure,
                                                                    projectPressurefterRestriction,
                                                                    true,
                                                                    false,
                                                                    normCalculationLevelIncrement,
                                                                    solveWithCoarseGridSolverOnEachFMGLevel,
                                                                    false,
                                                                    benchmarkName,
                                                                    verbose,
                                                                    dbFileDiscr );
      }

      solveImplementation< P2P1TaylorHoodFunction,
                           P2P1TaylorHoodStokesOperator,
                           P2ConstantMassOperator,
                           P1ConstantMassOperator,
                           P2P1StokesToP2P1StokesRestriction,
                           P2P1StokesToP2P1StokesProlongation,
                           P2P1StokesToP2P1StokesProlongation >( storage,
                                                                 uSolution,
                                                                 solutionU,
                                                                 solutionV,
                                                                 solutionW,
                                                                 solutionP,
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
                                                                 false,
                                                                 calculateDiscretizationError,
                                                                 normCalculationLevelIncrement,
                                                                 solveWithCoarseGridSolverOnEachFMGLevel,
                                                                 vtk,
                                                                 benchmarkName,
                                                                 verbose,
                                                                 dbFile );
   }
   else
   {
      P1StokesFunction< real_t > uSolution( "uSolution", storage, minLevel, maxLevel + normCalculationLevelIncrement );

      if ( calculateDiscretizationError )
      {
         solveImplementation< P1StokesFunction,
                              P1P1StokesOperator,
                              P1ConstantMassOperator,
                              P1ConstantMassOperator,
                              P1P1StokesToP1P1StokesRestriction,
                              P1P1StokesToP1P1StokesProlongation,
                              P1P1StokesToP1P1StokesProlongation >( storage,
                                                                    uSolution,
                                                                    solutionU,
                                                                    solutionV,
                                                                    solutionW,
                                                                    solutionP,
                                                                    rhsU,
                                                                    rhsV,
                                                                    rhsW,
                                                                    minLevel,
                                                                    maxLevel,
                                                                    multigridSettingsDiscrError,
                                                                    smootherSettingsDiscrError,
                                                                    coarseGridSettingsDiscrError,
                                                                    projectPressure,
                                                                    projectPressurefterRestriction,
                                                                    true,
                                                                    false,
                                                                    normCalculationLevelIncrement,
                                                                    solveWithCoarseGridSolverOnEachFMGLevel,
                                                                    false,
                                                                    benchmarkName,
                                                                    verbose,
                                                                    dbFileDiscr );
      }

      solveImplementation< P1StokesFunction,
                           P1P1StokesOperator,
                           P1ConstantMassOperator,
                           P1ConstantMassOperator,
                           P1P1StokesToP1P1StokesRestriction,
                           P1P1StokesToP1P1StokesProlongation,
                           P1P1StokesToP1P1StokesProlongation >( storage,
                                                                 uSolution,
                                                                 solutionU,
                                                                 solutionV,
                                                                 solutionW,
                                                                 solutionP,
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
                                                                 false,
                                                                 calculateDiscretizationError,
                                                                 normCalculationLevelIncrement,
                                                                 solveWithCoarseGridSolverOnEachFMGLevel,
                                                                 vtk,
                                                                 benchmarkName,
                                                                 verbose,
                                                                 dbFile );
   }
}

} // namespace tme_benchmarks
} // namespace hyteg
