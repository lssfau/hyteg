/*
 * Copyright (c) 2024 Michael Zikeli.
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
#pragma once

#ifndef DEACTIVATE_FLOATING_CHECKS
#include <cfenv>
#endif // n-def DEACTIVATE_FLOATING_CHECKS
#include <functional>
#include <iostream>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/Timer.h"
#include "core/timing/TimingTree.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/LaTeX/Table.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "Setups.h"
#include "SolverRegister.h"
#include "Tree.h"

namespace hyteg::mixedPrecisionMT {

/// @brief Update the residual res = Ax - b
template < class Operator_t, class Function_t >
void updateResidual( Function_t& res, const Operator_t& A, const Function_t& x, const Function_t& b, const uint_t level )
{
   A.apply( x, res, level, DoFType::Inner, Replace );
   res.assign( { 1, -1 }, { res, b }, level, DoFType::Inner );
}

/// @brief Compute the point wise norm ||r|| = sqrt( <r,r> / #DOFs )
template < class Function_t >
[[nodiscard]] static typename FunctionTrait< Function_t >::ValueType
    computeVectorNorm( const Function_t& r, const uint_t level, const bool computeSquaredNorm = false )
{
   using ValueType    = typename FunctionTrait< Function_t >::ValueType;
   auto squaredVector = static_cast< ValueType >( r.dotGlobal( r, level, DoFType::Inner ) );
   if ( computeSquaredNorm )
   {
      return squaredVector;
   }
   else
   {
      return static_cast< ValueType >( std::sqrt( squaredVector ) );
   }
}

/// @brief Compute the point wise norm ||r|| = sqrt( <r,r> / #DOFs )
template < class Function_t >
[[nodiscard]] static typename FunctionTrait< Function_t >::ValueType computeNormalizedVectorNorm( const Function_t& r,
                                                                                                  const uint_t      level )
{
   using ValueType = typename FunctionTrait< Function_t >::ValueType;
   auto vectorNorm = computeVectorNorm( r, level, true );
   return static_cast< ValueType >(
       std::sqrt( vectorNorm / numeric_cast< ValueType >( r.getNumberOfGlobalInnerDoFs( level ) ) ) );
}

///// @brief This function is called, when the program terminates, either because it is done or the outer solve run into an exeption.
template < class TimingTree,
           class SetupType
#ifndef PERFORMANCE_RUN
           ,
           class TreeNode,
           class TreeNodeReal
#endif // n-def PERFORMANCE_RUN
           >
static inline void writeBackAllInformation( const TimingTree&  tt,
                                            const SetupType&   config,
                                            const std::string& fileName
#ifndef PERFORMANCE_RUN
                                            ,
                                            const TreeNode&     pointWiseResidualTree,
                                            const TreeNode&     relativeResidualTree,
                                            const TreeNodeReal& pointWiseErrorTree
#endif // n-def PERFORMANCE_RUN
)
{
   hyteg::writeTimingTreeJSON( *tt, ( std::stringstream( "" ) << config.timingPath << "/timing_" << fileName << ".json" ).str() );

#ifndef PERFORMANCE_RUN
   WALBERLA_ROOT_SECTION()
   {
      std::ofstream pwFile( "./output/point-wise-residuals_" + fileName + ".json", std::ios::trunc );
      pwFile << pointWiseResidualTree.toJSON();
      std::ofstream relFile( "./output/relative-residuals_" + fileName + ".json", std::ios::trunc );
      relFile << relativeResidualTree.toJSON();
      std::ofstream errFile( "./output/point-wise-errors_" + fileName + ".json", std::ios::trunc );
      errFile << pointWiseErrorTree.toJSON();
   }
#endif // n-def PERFORMANCE_RUN
}

template < class SetupType, uint_t Quad, SolverRegister SolverKind, class SmootherOperator_t, class ResidualOperator_t >
static inline void solveLSE( const SetupType& config )
{
   // +++ Initialize Output +++
#ifdef WRITE_BENCHMARK
   hyteg::latex::Table< 17 > benchmarkTable( { "Dim",
                                               "bar",
                                               "dot",
                                               "Level",
                                               "Err-Conv-Rate",
                                               "Max-Error",
                                               "Point-Wise-Error",
                                               "Pseudo-Relative-Error",
                                               "Weak-Conv-Rate",
                                               "Error",
                                               "Relative-Error",
                                               "Point-Wise-Residual",
                                               "Relative-Residual",
                                               "Res-Conv-Rate",
                                               "Runtime",
                                               "Solver-Iterations",
                                               "Termination-Criteria" } );
#endif

   // +++ Define Names +++
   using Function_bar_t  = typename ResidualOperator_t::srcType;
   using Function_real_t = typename Function_bar_t::template FunctionType< real_t >;
   using Value_bar_t     = typename Function_bar_t::valueType;
   using Value_dot_t     = typename SmootherOperator_t::srcType::valueType;

   SetupPrimitiveStorage setupStorage( MeshInfo::fromGmshFile( config.problemSetup.meshFile ),
                                       uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // +++ Setup Mesh +++
   const uint_t minLevel  = config.minLevel;
   const uint_t maxLevel  = config.maxLevel;
   const uint_t errLevels = config.errorCalcLevels;
   auto         storage   = std::make_shared< PrimitiveStorage >( setupStorage );

   // +++ Get timing tree ++
   auto tt = storage->getTimingTree();

   // +++ Initialize Operator | From here the different Operators are compared with each other +++
   WALBERLA_LOG_INFO_ON_ROOT(
       "\n Start Benchmark of "
       << "\n\tName of PDE              :" << config.problemSetup.descriptionPDE
       << "\n\tName of Solver           :" << ( SolverTrait< SolverKind, ResidualOperator_t, SetupType >::toString )
       << "\n\tName of SmootherOperator :" << typeid( SmootherOperator_t ).name() << "\n\tName of ResidualOperator :"
       << typeid( ResidualOperator_t ).name() << "\n\tName of SmootherType     :" << AccuracyTrade< Value_dot_t >::name()
       << "\n\tName of ResidualType     :" << AccuracyTrade< Value_bar_t >::name()
       << "\n\tMin Level                :" << config.minLevel << "\n\tMax Level                :" << config.maxLevel
       << "\n\tConfig                   :" << config.configDescriptor
       << ".\n==========================================================================================\n" );

   std::string nameOfRun =
       ( std::stringstream( "" ) << ( ( config.problemSetup.dim3 ) ? 3 : 2 ) << "D"
                                 << "_R-" << AccuracyTrade< Value_bar_t >::name() << "_S-" << AccuracyTrade< Value_dot_t >::name()
                                 << "_" << config.problemSetup.fileDescription << "_"
                                 << SolverTrait< SolverKind, ResidualOperator_t, SetupType >::toString )
           .str();

   const uint_t minAllocationLevel = std::min( minLevel, config.coarseLevelGMG );

   tt->start( nameOfRun + config.configDescriptor );

   tt->start( "Initial Memory Allocation" );
   // +++ Initialize Functions +++
   // Since this benchmark should work for all kinds of solvers, all must be allocated on all levels. (GMG uses x and b on all levels)
   Function_bar_t x( "x", storage, minAllocationLevel, maxLevel );
   Function_bar_t b( "b", storage, minAllocationLevel, maxLevel );
   Function_bar_t res( "res", storage, minAllocationLevel, maxLevel );
   // <Below> is necessary, since I compare this on different levels to get the conRate // For output reasons, I want this in the current level as well.
   Function_real_t xExact( "xExact", storage, minLevel, uint_c( maxLevel + errLevels ) );
   Function_real_t err( "err", storage, minLevel, uint_c( maxLevel + errLevels ) );
   tt->stop( "Initial Memory Allocation" );

#ifndef PERFORMANCE_RUN
   // +++ Create structure to store residuals in +++
   Tree< Value_bar_t > pointWiseResidualTree        = {};
   Tree< Value_bar_t > relativeResidualTree         = {};
   Tree< real_t >      pointWiseErrorTree           = {};
   auto                pointWiseResidualRootNodePtr = pointWiseResidualTree.insert(
       static_cast< Value_bar_t >( 0 ), nameOfRun + config.configDescriptor, "point-wise", pointWiseResidualTree.getRoot() );
   auto relativeResidualRootNodePtr = relativeResidualTree.insert(
       static_cast< Value_bar_t >( 0 ), nameOfRun + config.configDescriptor, "relative", relativeResidualTree.getRoot() );
   auto pointWiseErrorRootNodePtr = pointWiseErrorTree.insert(
       real_c( 0 ), nameOfRun + config.configDescriptor, "point-wise-err", pointWiseErrorTree.getRoot() );
#endif // n-def PERFORMANCE_RUN

   tt->start( "Initialize Problem" );
   // Note x is initialized within the loop
   for ( uint_t level = minLevel; level <= uint_c( maxLevel + errLevels ); ++level )
   {
      if ( level <= maxLevel ) // skip error calculation levels
      {
         // The inner DoFs of x are initialized within the max level loop
         Function_real_t bInitL2Space( "bInitL2Space", storage, level, level );
         Function_real_t xRealForCopy( "xRealForCopy", storage, level, level );
         xRealForCopy.interpolate( config.problemSetup.initExact, level, DoFType::DirichletBoundary );
         x.copyFrom( xRealForCopy, level );
         // Initialize the weak form of the RHS function and copy it to the one in the right precision.
         L2Space< Quad, Function_real_t, real_t > L2( storage, level );
         L2.dot( config.problemSetup.initRHS, bInitL2Space ); // compute rhs using quadrature
         b.copyFrom( bInitL2Space, level );
      }
      xExact.interpolate( config.problemSetup.initExact, level, DoFType::All ); // Done at all levels
   }
   tt->stop( "Initialize Problem" );

#ifdef WRITE_VTK
   // +++ create temporary functions to convert function to ValueTypes that can be printed by the vtk printer +++
   Function_real_t printable_x( "x", storage, minLevel, maxLevel );
   Function_real_t printable_xExact( "xExact", storage, minLevel, maxLevel );
   Function_real_t printable_b( "b", storage, minLevel, maxLevel );
   Function_real_t printable_res( "res", storage, minLevel, maxLevel );
#endif // WRITE_VTK

   real_t                  pointWiseErrorNorm = -std::numeric_limits< real_t >::max();
   real_t                  lastErrorNorm      = -std::numeric_limits< real_t >::max();
   real_t                  errorRate          = -std::numeric_limits< real_t >::max();
   [[maybe_unused]] real_t lastWeakNorm       = -std::numeric_limits< real_t >::max();
   real_t                  weakRate           = -std::numeric_limits< real_t >::max();
   [[maybe_unused]] uint_t divergenceCounter  = 0;
   for ( uint_t currentMaxLevel = minLevel; currentMaxLevel <= maxLevel; ++currentMaxLevel )
   {
      tt->start( "maxLevel_" + std::to_string( currentMaxLevel ) );

      tt->start( "Initialization Outer solver" );
      uint_t operatorMinLvl = minAllocationLevel;
      if constexpr ( SolverKind == SolverRegister::IR )
      {
         operatorMinLvl = currentMaxLevel;
      }
      ResidualOperator_t laplaceOperatorGenBar( storage, operatorMinLvl, currentMaxLevel );
      SmootherOperator_t laplaceOperatorGenDot( storage, minAllocationLevel, currentMaxLevel );

      // +++ Initialize Solver +++
      auto solver = SolverTrait< SolverKind, ResidualOperator_t, SetupType, SmootherOperator_t >::getInitializedSolver(
          laplaceOperatorGenBar, config, storage, minAllocationLevel, currentMaxLevel, laplaceOperatorGenDot );

      // +++ Initialize Functions +++
      x.interpolate( AccuracyTrade< Value_bar_t >::initZero,
                     currentMaxLevel,
                     DoFType::Inner ); // All algorithms only work on the finest level of x.

      updateResidual( res, laplaceOperatorGenBar, x, b, currentMaxLevel );
      auto initResidualVectorNorm = computeVectorNorm( res, currentMaxLevel );
      if ( std::numeric_limits< Value_bar_t >::min() > initResidualVectorNorm )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "\n\t\033[33m"
                                       << "On level " << currentMaxLevel
                                       << " the initial residual is already to small to continue!"
                                       << "\033[0m\n" );
         tt->stop( "Initialization Outer solver" );
         tt->stop( "maxLevel_" + std::to_string( currentMaxLevel ) );
         break;
      }

      // +++ Create holder variable for residual to calculate residual convergence +++
      Value_bar_t relativeResidualNorm   = computeVectorNorm( res, currentMaxLevel ) / initResidualVectorNorm;
      Value_bar_t pointWiseResidualNorm  = computeNormalizedVectorNorm( res, currentMaxLevel );
      auto        getCurrentResidualNorm = [&pointWiseResidualNorm, &relativeResidualNorm, norm_t = config.residualNormType]() {
         return chooseRightNorm( norm_t, relativeResidualNorm, pointWiseResidualNorm );
      };
      Value_bar_t lastResidualNorm         = getCurrentResidualNorm();
      auto        lastErrorNormForStopping = real_c( std::sqrt( xExact.dotGlobal( xExact, currentMaxLevel, DoFType::Inner ) /
                                                         real_c( xExact.getNumberOfGlobalInnerDoFs( currentMaxLevel ) ) ) );

#ifndef PERFORMANCE_RUN
      WALBERLA_ROOT_SECTION()
      {
         // +++ Info for Chebychev Coefficients +++
         std::ofstream ofCheby( "ChebyCoeffs.dat", std::ios::app );
         ofCheby << "|\t" << nameOfRun << config.configDescriptor << "\t" << currentMaxLevel << "\n";
      }
      WALBERLA_MPI_BARRIER();

      auto pointWiseCurrentLevelInitialResidualNodePtr =
          pointWiseResidualTree.insert( static_cast< Value_bar_t >( pointWiseResidualNorm ),
                                        "level_" + std::to_string( currentMaxLevel ),
                                        "initial-residual",
                                        pointWiseResidualRootNodePtr );
      auto relativeCurrentLevelInitialResidualNodePtr =
          relativeResidualTree.insert( static_cast< Value_bar_t >( relativeResidualNorm ),
                                       "level_" + std::to_string( currentMaxLevel ),
                                       "initial-residual",
                                       relativeResidualRootNodePtr );
      [[maybe_unused]] auto pointWiseIterationResidualNodePtr =
          pointWiseResidualTree.insert( static_cast< Value_bar_t >( pointWiseResidualNorm ),
                                        "iteration_" + std::to_string( 0 ),
                                        "residual-per-iteration",
                                        pointWiseCurrentLevelInitialResidualNodePtr );
      [[maybe_unused]] auto relativeIterationResidualNodePtr =
          relativeResidualTree.insert( static_cast< Value_bar_t >( relativeResidualNorm ),
                                       "iteration_" + std::to_string( 0 ),
                                       "residual-per-iteration",
                                       relativeCurrentLevelInitialResidualNodePtr );

      auto pointWiseCurrentLevelInitialErrorNodePtr = pointWiseErrorTree.insert( real_c( lastErrorNormForStopping ),
                                                                                 "level_" + std::to_string( currentMaxLevel ),
                                                                                 "initial-error",
                                                                                 pointWiseErrorRootNodePtr );
      [[maybe_unused]] auto pointWiseIterationErrorNodePtr =
          pointWiseErrorTree.insert( real_c( lastErrorNormForStopping ),
                                     "iteration_" + std::to_string( 0 ),
                                     "error-per-iteration",
                                     pointWiseCurrentLevelInitialErrorNodePtr );

      if constexpr ( SolverKind == SolverRegister::IR )
      // Set the residual node on which the correction norm should be appended
      {
         solver->setResidualNode( pointWiseIterationResidualNodePtr, &pointWiseResidualTree );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "\nFor " << currentMaxLevel << " being the finest level, "
                                          << "the initial residual is: " << getCurrentResidualNorm() << ",\t"
                                          << "the initial error is: " << lastErrorNormForStopping );
#endif // n-def PERFORMANCE_RUN

#ifdef WRITE_VTK
      VTKOutput vtkOutput( "./output", std::to_string( currentMaxLevel ) + "_" + nameOfRun + config.configDescriptor, storage );
      vtkOutput.add( printable_x );
      vtkOutput.add( printable_b );
      vtkOutput.add( printable_res );
      vtkOutput.add( printable_xExact );

      if ( config.printVTK )
      {
         // +++ convert functions to print in a value type, that can be printed by the vtk printer +++
         printable_x.copyFrom( x, currentMaxLevel );
         printable_b.copyFrom( b, currentMaxLevel );
         printable_res.copyFrom( res, currentMaxLevel );
         printable_xExact.copyFrom( xExact, currentMaxLevel );
         vtkOutput.write( currentMaxLevel, 0 );
      }
#endif // WRITE_VTK

      tt->stop( "Initialization Outer solver" );

      walberla::WcTimer timeToSolution;

      // +++ Solve PDE +++
      tt->start( "OuterLoop" );
      uint_t solverIterations    = 0;
      auto   InverseResidualRate = real_t( getCurrentResidualNorm() / lastResidualNorm ); // The initial value is 1
      auto   InverseErrorRate    = real_t( 1. );
      uint_t terminationCriteria = 0;
      do
      {
         // +++ Do one solver step +++
         tt->start( "OuterSolve" );
         LIKWID_MARKER_START( "OuterSolve" );
         solver->solve( laplaceOperatorGenBar, x, b, currentMaxLevel );
         LIKWID_MARKER_STOP( "OuterSolve" );
         tt->stop( "OuterSolve" );

         tt->start( "OuterResidual" );
         // +++ Calculate Residual Norm +++
         updateResidual( res, laplaceOperatorGenBar, x, b, currentMaxLevel );
         relativeResidualNorm  = computeVectorNorm( res, currentMaxLevel ) / initResidualVectorNorm;
         pointWiseResidualNorm = computeNormalizedVectorNorm( res, currentMaxLevel );

         // +++ Update Residual and Error Rate for current maxLevel +++
         InverseResidualRate = real_c( getCurrentResidualNorm() / lastResidualNorm );
         lastResidualNorm    = getCurrentResidualNorm();
         tt->stop( "OuterResidual" );

         tt->start( "OuterError" );
         Function_real_t xErrorForStopping( "xForErrorFS", storage, currentMaxLevel, currentMaxLevel );
         xErrorForStopping.copyFrom( x, currentMaxLevel );
         err.assign( { 1, -1 }, { xErrorForStopping, xExact }, currentMaxLevel, DoFType::Inner );
         pointWiseErrorNorm = computeNormalizedVectorNorm( err, currentMaxLevel );

         InverseErrorRate         = pointWiseErrorNorm / lastErrorNormForStopping;
         lastErrorNormForStopping = pointWiseErrorNorm;
         tt->stop( "OuterError" );

         solverIterations++;
#ifndef PERFORMANCE_RUN
         pointWiseIterationResidualNodePtr = pointWiseResidualTree.insert( static_cast< Value_bar_t >( pointWiseResidualNorm ),
                                                                           "iteration_" + std::to_string( solverIterations ),
                                                                           "residual-per-iteration",
                                                                           pointWiseCurrentLevelInitialResidualNodePtr );
         relativeIterationResidualNodePtr  = relativeResidualTree.insert( static_cast< Value_bar_t >( relativeResidualNorm ),
                                                                         "iteration_" + std::to_string( solverIterations ),
                                                                         "residual-per-iteration",
                                                                         relativeCurrentLevelInitialResidualNodePtr );
         pointWiseIterationErrorNodePtr    = pointWiseErrorTree.insert( real_c( pointWiseErrorNorm ),
                                                                     "iteration_" + std::to_string( solverIterations ),
                                                                     "error-per-iteration",
                                                                     pointWiseCurrentLevelInitialErrorNodePtr );

         WALBERLA_LOG_INFO_ON_ROOT( "IT " << solverIterations << ".\tInverse Residual Rate: " << InverseResidualRate
                                          << ",\tresidual : " << getCurrentResidualNorm() << ".\tInverse Error Rate: "
                                          << InverseErrorRate << ",\terror : " << pointWiseErrorNorm );
#endif // n-def PERFORMANCE_RUN

#ifdef WRITE_VTK
         if ( config.printVTK )
         {
            // +++ convert functions to print in a value type, that can be printed by the vtk printer +++
            printable_x.copyFrom( x, currentMaxLevel );
            printable_b.copyFrom( b, currentMaxLevel );
            printable_res.copyFrom( res, currentMaxLevel );
            printable_xExact.copyFrom( xExact, currentMaxLevel );
            vtkOutput.write( currentMaxLevel, solverIterations );
         }
#endif // WRITE_VTK
      } while ( config.needFurtherOuterSolveIterations( getCurrentResidualNorm(),
                                                        InverseResidualRate,
                                                        lastErrorNormForStopping,
                                                        InverseErrorRate,
                                                        solverIterations,
                                                        &terminationCriteria ) );

      timeToSolution.end();

      auto residualRate = 1. / ( abs( InverseResidualRate ) + std::numeric_limits< real_t >::epsilon() );
      errorRate         = 1. / ( abs( InverseErrorRate ) + std::numeric_limits< real_t >::epsilon() );
      tt->stop( "OuterLoop" );

#ifndef PERFORMANCE_RUN
      [[maybe_unused]] auto pointWiseCurrentLevelResidualNodePtr =
          pointWiseResidualTree.insert( static_cast< Value_bar_t >( pointWiseResidualNorm ),
                                        "level_" + std::to_string( currentMaxLevel ),
                                        "final-residual",
                                        pointWiseResidualRootNodePtr );
      [[maybe_unused]] auto relativeCurrentLevelResidualNodePtr =
          relativeResidualTree.insert( static_cast< Value_bar_t >( relativeResidualNorm ),
                                       "level_" + std::to_string( currentMaxLevel ),
                                       "final-residual",
                                       relativeResidualRootNodePtr );
      [[maybe_unused]] auto pointWiseCurrentLevelErrorNodePtr =
          pointWiseErrorTree.insert( real_c( lastErrorNormForStopping ),
                                     "level_" + std::to_string( currentMaxLevel ),
                                     "final-error",
                                     pointWiseErrorRootNodePtr );
#endif // n-def PERFORMANCE_RUN

      // +++ Calculate Error in Sobolev space +++
      real_t relativeErrorNorm = computeVectorNorm( err, currentMaxLevel );
      real_t relativeWeakNorm  = -std::numeric_limits< real_t >::max();
      real_t weakNorm          = -std::numeric_limits< real_t >::max();
#ifndef PERFORMANCE_RUN
      tt->start( "Compute error" );
      // For the Hockey-Stick Runs, the function is prolongated to reduce sideeffects effects.
#ifdef USE_WEAK_NORM
      {
         L2Space< Quad, Function_real_t, real_t > L2( storage, currentMaxLevel );
         Function_real_t                          xForError( "xForError", storage, currentMaxLevel, currentMaxLevel );
         xForError.copyFrom( x, currentMaxLevel );
         auto computeL2Error = [&]( const Point3D& localPoint ) {
            real_t ux = 0;
            WALBERLA_CHECK( xForError.evaluate( localPoint, currentMaxLevel, ux ) );
            real_t e = config.problemSetup.initExact( localPoint ) - ux;
            return e;
         };

         weakNorm          = L2.norm( computeL2Error );
         auto relationNorm = L2.norm( config.problemSetup.initExact );

         // Relative Weak Norm
         if ( abs( relationNorm ) < std::numeric_limits< real_t >::min() )
         {
            if ( abs( weakNorm ) < std::numeric_limits< real_t >::min() )
            {
               relativeWeakNorm = real_t( 0.0 );
            }
            else
            {
               relativeWeakNorm = std::numeric_limits< real_t >::infinity();
            }
         }
         else
         {
            relativeWeakNorm = real_c( weakNorm / relationNorm );
         }

         // Relative Vector Norm
         if ( abs( relationNorm ) < std::numeric_limits< real_t >::min() )
         {
            if ( abs( relativeErrorNorm ) < std::numeric_limits< real_t >::min() )
            {
               relativeErrorNorm = real_t( 0.0 );
            }
            else
            {
               relativeErrorNorm = std::numeric_limits< real_t >::infinity();
            }
         }
         else
         {
            relativeErrorNorm = real_c( weakNorm / relationNorm );
         }
      }

      errorRate     = lastErrorNorm / chooseRightNorm( config.errorNormType, relativeErrorNorm, pointWiseErrorNorm );
      weakRate      = lastWeakNorm / relativeWeakNorm;
      lastErrorNorm = chooseRightNorm( config.errorNormType, relativeErrorNorm, pointWiseErrorNorm );
      lastWeakNorm  = relativeWeakNorm;
#else  // n-def USE_WEAK_NORM
      errorRate     = lastErrorNorm / pointWiseErrorNorm;
      lastErrorNorm = pointWiseErrorNorm;
#endif // def USE_WEAK_NORM

      tt->stop( "Compute error" );
#else  // n-def PERFORMANCE_RUN
      errorRate     = lastErrorNorm / pointWiseErrorNorm;
      lastErrorNorm = pointWiseErrorNorm;
#endif // n-def PERFORMANCE_RUN
      real_t maxError = err.getMaxMagnitude( currentMaxLevel, DoFType::Inner );

      tt->stop( "maxLevel_" + std::to_string( currentMaxLevel ) );

      // +++ Write Output +++
#ifdef WRITE_BENCHMARK
      benchmarkTable.pushRow( ( config.problemSetup.dim3 ) ? 3 : 2,
                              AccuracyTrade< Value_bar_t >::name(),
                              AccuracyTrade< Value_dot_t >::name(),
                              currentMaxLevel,
                              errorRate,
                              maxError,
                              pointWiseErrorNorm,
                              relativeErrorNorm,
                              weakRate,
                              weakNorm,
                              relativeWeakNorm,
                              pointWiseResidualNorm,
                              relativeResidualNorm,
                              residualRate,
                              timeToSolution.last(),
                              solverIterations,
                              config.matchTermination( terminationCriteria ) );
      WALBERLA_ROOT_SECTION()
      {
         benchmarkTable.writeUpdate( config.benchPath, nameOfRun + config.configDescriptor );
      }
#endif

#ifndef PERFORMANCE_RUN
      WALBERLA_LOG_INFO_ON_ROOT(
          "\nResidualOperator                  :"
          << typeid( ResidualOperator_t ).name()
          << "\nValue_bar_t                       :" << AccuracyTrade< Value_bar_t >::name()
          << "\nSmootherOperator                  :" << typeid( SmootherOperator_t ).name()
          << "\nValue_dot_t                       :" << AccuracyTrade< Value_dot_t >::name()
          << "\nMax Level                         :" << currentMaxLevel << "\nMax                 error    is   :" << maxError
          << "\nPoint-Wise       L2 error    is   :" << pointWiseErrorNorm
          << "\nPseudo-Relative  L2 error    is   :" << relativeErrorNorm << "\nPoint-Wise       L2 error    rate :" << errorRate
          << "\nWeak             L2 error    is   :" << weakNorm << "\nRelative         L2 error    is   :" << relativeWeakNorm
          << "\n                 L2 error    rate :" << weakRate << "\nPoint-Wise       L2 residual is   :"
          << pointWiseResidualNorm << "\nRelative         L2 residual is   :" << relativeResidualNorm
          << "\n                 L2 residual rate :" << residualRate << "\nRuntime                           :"
          << timeToSolution.last() << "\n#Iterations are                   :" << solverIterations << "\n" );
#endif // n-def PERFORMANCE_RUN

      WALBERLA_MPI_BARRIER()

#ifndef PERFORMANCE_RUN
      // +++ If discretization error starts to increase (more than fluctuation), terminate run +++
      if ( errorRate < 3.5 )
      {
         divergenceCounter++;
      }
      else
      {
         divergenceCounter = 0;
      }
      // +++ If too many iterations were necessary to solve PDE on this level, stop refining +++
      if ( std::isnan( chooseRightNorm( config.errorNormType, relativeErrorNorm, pointWiseErrorNorm ) ) ||
           divergenceCounter > config.maxUnnecessaryRefinements )
      {
         break;
      }
#endif // n-def PERFORMANCE_RUN

   } // loop for currentMaxLevel

   tt->stop( nameOfRun + config.configDescriptor );

   writeBackAllInformation( tt,
                            config,
                            nameOfRun + config.configDescriptor
#ifndef PERFORMANCE_RUN
                            ,
                            pointWiseResidualTree,
                            relativeResidualTree,
                            pointWiseErrorTree
#endif // n-def PERFORMANCE_RUN
   );

   WALBERLA_MPI_BARRIER()

} // hyteg::mixedPrecisionMT::solveLSE()

// ==============================================================================================================================

template < class SetupType, uint_t Quad, SolverRegister SolverKind, class SmootherOperator_t, class... ResidualOperators_t >
inline void solveLSEWithOperatorFold( const SetupType& config )
{
#if defined( PERFORMANCE_RUN )
#ifdef WRITE_VTK
#warning WRITE_VTK is not compatible with PERFORMANCE_RUN! WRITE_VTK will be undefineddefined.
#undef WRITE_VTK
#endif // def WRITE_VTK
#ifndef DEACTIVATE_FLOATING_CHECKS
#warning For performance reasons, all arithmetic checks are disabled.
#define DEACTIVATE_FLOATING_CHECKS
#endif // n-def DEACTIVATE_FLOATING_CHECKS
#endif // n-def PERFORMANCE_RUN

#ifndef DEACTIVATE_FLOATING_CHECKS
   feenableexcept( FE_DIVBYZERO
#ifdef CHECK_FOR_UNDERFLOW
                   | FE_UNDERFLOW
#endif // def CHECK_FOR_UNDERFLOW
   );
#endif // n-def DEACTIVATE_FLOATING_CHECKS

//#pragma message ( "The defines are set as following:" )
#ifdef WRITE_BENCHMARK
   //   #pragma message( "'WRITE_BENCHMARK'            : ON" )
   WALBERLA_LOG_INFO_ON_ROOT( "'WRITE_BENCHMARK'            : ON" );
#else  // WRITE_BENCHMARK
   //   #pragma message( "'WRITE_BENCHMARK'            : OFF" )
   WALBERLA_LOG_INFO_ON_ROOT( "'WRITE_BENCHMARK'            : OFF" );
#endif // def WRITE_BENCHMARK
#ifdef WRITE_VTK
   //   #pragma message( "'WRITE_VTK'                  : ON" )
   WALBERLA_LOG_INFO_ON_ROOT( "'WRITE_VTK'                  : ON" );
#else  // n-def WRITE_VTK
   //   #pragma message( "'WRITE_VTK'                  : OFF" )
   WALBERLA_LOG_INFO_ON_ROOT( "'WRITE_VTK'                  : OFF" );
#endif // def WRITE_VTK
#ifdef PERFORMANCE_RUN
   //   #pragma message( "'PERFORMANCE_RUN'            : ON")
   WALBERLA_LOG_INFO_ON_ROOT( "'PERFORMANCE_RUN'            : ON" );
#else  // n-def PERFORMANCE_RUN
   //   #pragma message( "'PERFORMANCE_RUN'            : OFF")
   WALBERLA_LOG_INFO_ON_ROOT( "'PERFORMANCE_RUN'            : OFF" );
#endif // def PERFORMANCE_RUN
#ifdef USE_WEAK_NORM
   //   #pragma message( "'USE_WEAK_NORM'          : ON")
   WALBERLA_LOG_INFO_ON_ROOT( "'USE_WEAK_NORM'          : ON" );
#else  // n-def USE_WEAK_NORM
   //   #pragma message( "'USE_WEAK_NORM'          : OFF")
   WALBERLA_LOG_INFO_ON_ROOT( "'USE_WEAK_NORM'          : OFF" );
#endif // def USE_WEAK_NORM
#ifdef DEACTIVATE_FLOATING_CHECKS
   //   #pragma message( "'DEACTIVATE_FLOATING_CHECKS' : ON" )
   WALBERLA_LOG_INFO_ON_ROOT( "'DEACTIVATE_FLOATING_CHECKS' : ON" );
#else  // n-def DEACTIVATE_FLOATING_CHECKS
   //   #pragma message( "'DEACTIVATE_FLOATING_CHECKS' : OFF" )
   WALBERLA_LOG_INFO_ON_ROOT( "'DEACTIVATE_FLOATING_CHECKS' : OFF" );
#endif // def DEACTIVATE_FLOATING_CHECKS
#ifdef CHECK_FOR_UNDERFLOW
   //   #pragma message( "'CHECK_FOR_UNDERFLOW'        : ON" )
   WALBERLA_LOG_INFO_ON_ROOT( "'CHECK_FOR_UNDERFLOW'        : ON" );
#else  // n-def CHECK_FOR_UNDERFLOW
   //   #pragma message( "'CHECK_FOR_UNDERFLOW'        : OFF" )
   WALBERLA_LOG_INFO_ON_ROOT( "'CHECK_FOR_UNDERFLOW'        : OFF" );
#endif // def CHECK_FOR_UNDERFLOW

   LIKWID_MARKER_INIT;
   LIKWID_MARKER_THREADINIT;

   ( solveLSE< SetupType, Quad, SolverKind, SmootherOperator_t, ResidualOperators_t >( config ), ... );

   LIKWID_MARKER_CLOSE;
} // hyteg::mixedPrecisionMT::solveLSEWithOperatorFold()

} // namespace hyteg::mixedPrecisionMT