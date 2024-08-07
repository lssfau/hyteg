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

#include <functional>
#include <optional>

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "GridtransformationRegister.h"

// Solvers
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"

#include "IterativeRefinementSolver.hpp" // TODO make MR to add this to master
#include "globalHeader.h"

namespace hyteg::mixedPrecisionMT {

enum class SolverRegister
{
   ChebyshevSmoother,
   GMG,
   IR,
   CG,
   WeightedJacobiSmoother,
   SolverLoop
};
using Storage = std::shared_ptr< PrimitiveStorage >;

template < SolverRegister Solver, class ResidualOperator_t, class SetupType, class SmootherOperator_t = ResidualOperator_t >
struct SolverTrait;

template < class Solver_t, class Operator_t, class SetupType >
static inline std::shared_ptr< SolverLoop< Operator_t > > makeRepetitive( const std::shared_ptr< Solver_t >& solver,
                                                                          const SetupType&,
                                                                          Storage&     storage,
                                                                          const uint_t iterations,
                                                                          const real_t tolerance )
{
   using Function_t = Operator_t::srcType;
   auto stopIterationCallback =
       [&storage, tolerance]( const Operator_t& A, const Function_t& x, const Function_t& b, const uint_t level ) -> bool {
      Function_t res( "res", storage, level, level );
      A.apply( x, res, level, DoFType::Inner, Replace );
      res.assign( { 1, -1 }, { res, b }, level, DoFType::Inner );
      const real_t residual =
          std::sqrt( res.dotGlobal( res, level, DoFType::Inner ) /
                     walberla::numeric_cast< typename Function_t::valueType >( res.getNumberOfGlobalDoFs( level ) ) );

      return ( residual < tolerance );
   };
   return std::make_shared< SolverLoop< Operator_t > >( solver, iterations, stopIterationCallback );
}

// ==============================================================================================================================
// |                                           Definition of SolverTraits                                                     |
// ==============================================================================================================================

// ===== ChebyshevSmoother =====

template < class Operator_t, class SetupType >
struct SolverTrait< SolverRegister::ChebyshevSmoother, Operator_t, SetupType >
{
   typedef ChebyshevSmoother< Operator_t > Solver_t;
   static constexpr auto                   toString = "ChebyshevSmoother";
   using Function_t                                 = Operator_t::srcType;
   using ValueType                                  = Function_t::valueType;

   static std::shared_ptr< Solver_t > getInitializedSolver( Operator_t&      A,
                                                            const SetupType& config,
                                                            Storage&         storage,
                                                            uint_t           minLevel,
                                                            uint_t           maxLevel,
                                                            std::optional< std::reference_wrapper< Operator_t > > = std::nullopt )
   {
      const auto levelOfSpectralRadius =
          std::min( maxLevel, uint_t( 5 ) ); // At max -> level 5 or the minLevel if this is larger.
      Function_t initEV( "initEV", storage, levelOfSpectralRadius, levelOfSpectralRadius );
      Function_t tmp( "tmp", storage, levelOfSpectralRadius, levelOfSpectralRadius );

      initEV.interpolate( AccuracyTrade< ValueType >::initRandom, levelOfSpectralRadius, DoFType::All );
      A.computeInverseDiagonalOperatorValues();
      auto solver = std::make_shared< Solver_t >( storage, minLevel, maxLevel );

      //      {  // Check if Inverse Diagonal Values are zero or NaN
      //         auto inverseDiagonalValues = A.getInverseDiagonalValues();
      //         const real_t localNormSqr = inverseDiagonalValues->dotLocal( *inverseDiagonalValues, levelOfSpectralRadius, (Inner | NeumannBoundary) );
      //         if ( localNormSqr <= 0.0 || std::isnan(localNormSqr) )
      //         {
      //            WALBERLA_LOG_WARNING_ON_ROOT( "\n\t\033[33m" << "Underflow: InverseDiagonal of Dot Operator." << "\033[0m\n" );
      //            // TODO find a nice solution
      //         }
      //      }
      auto spectralRadius = chebyshev::estimateRadius( A, levelOfSpectralRadius, maxIterForInvPow_, storage, initEV, tmp );
      solver->setupCoefficients( config.polynomialOrderCheby, spectralRadius );

      return solver;
   }

 private:
   static constexpr auto maxIterForInvPow_ = uint_t( 20 );
};

// ===== CG =====

template < class Operator_t, class SetupType >
struct SolverTrait< SolverRegister::CG, Operator_t, SetupType >
{
   typedef CGSolver< Operator_t > Solver_t;
   static constexpr auto          toString = "CG";
   using FunctionType                      = typename Operator_t::srcType;
   using ValueType                         = typename FunctionTrait< FunctionType >::ValueType;

   static std::shared_ptr< Solver_t > getInitializedSolver( Operator_t&,
                                                            const SetupType& config,
                                                            Storage&         storage,
                                                            uint_t           minLevel,
                                                            uint_t           maxLevel,
                                                            std::optional< std::reference_wrapper< Operator_t > > = std::nullopt )
   {
      auto solver = std::make_shared< Solver_t >(
          storage,
          minLevel,
          maxLevel,
          config.cgIterations,
          //                                           numeric_cast<ValueType>( std::max( residualFactor_ * std::numeric_limits<ValueType>::min(),
          //                                                                                std::numeric_limits<walberla::float64>::epsilon() ) )
          numeric_cast< ValueType >( residualFactor_ * std::numeric_limits< ValueType >::min() ) );
#ifndef PERFORMANCE_RUN
      solver->setPrintInfo( true );
#endif // n-def PERFORMANCE_RUN
      return solver;
   }

 private:
   static constexpr auto residualFactor_ = real_t( 10.0 );
};

// ===== GMG =====

template < class Operator_t, class SetupType >
struct SolverTrait< SolverRegister::GMG, Operator_t, SetupType >
{
   typedef GeometricMultigridSolver< Operator_t > Solver_t;
   static constexpr auto                          toString = "GMG";
   using FunctionType                                      = typename Operator_t::srcType;
   using ValueType                                         = typename FunctionTrait< FunctionType >::ValueType;

   static std::shared_ptr< Solver_t > getInitializedSolver( Operator_t&      A,
                                                            const SetupType& config,
                                                            Storage&         storage,
                                                            uint_t           minLevel,
                                                            uint_t           maxLevel,
                                                            std::optional< std::reference_wrapper< Operator_t > > = std::nullopt )
   {
      // +++ Setup of Smoother +++
      auto smoother =
          SolverTrait< solver_, Operator_t, SetupType >::getInitializedSolver( A, config, storage, minLevel, maxLevel );

      const auto coarseLevel = std::min( config.coarseLevelGMG, maxLevel );

      std::shared_ptr< Solver< Operator_t > > coarseSolver;
      if constexpr ( useSmootherForCoarseGridAsWell_ )
      {
         auto CoarseSmoother =
             SolverTrait< solver_, Operator_t, SetupType >::getInitializedSolver( A, config, storage, coarseLevel, coarseLevel );
         using SmootherType = SolverTrait< solver_, Operator_t, SetupType >::Solver_t;
         coarseSolver       = makeRepetitive< SmootherType, Operator_t, SetupType >(
             smoother, config, storage, uint_c( 1000 ), config.residualThreshold );
      }
      else
      {
         coarseSolver = SolverTrait< coarseSolver_, Operator_t, SetupType >::getInitializedSolver(
             A, config, storage, coarseLevel, coarseLevel );
      }

      // +++ Setup of Grid Transformations +++
      auto prolongation =
          GridTransformationTrait< ProlongationRegister::P1Linear, RestrictionRegister::P1Linear, ValueType >::prolongation();
      auto restriction =
          GridTransformationTrait< ProlongationRegister::P1Linear, RestrictionRegister::P1Linear, ValueType >::restriction();

      // +++ Setup of actual Solver +++
      return std::make_shared< Solver_t >( storage,
                                           smoother,
                                           coarseSolver,
                                           restriction,
                                           prolongation,
                                           coarseLevel,
                                           maxLevel,
                                           config.preSmoothingStepsGMG,
                                           config.postSmoothingStepsGMG );
   }

 private:
   static constexpr auto useSmootherForCoarseGridAsWell_ = false;
   static constexpr auto solver_                         = SolverRegister::ChebyshevSmoother;
   static constexpr auto coarseSolver_                   = SolverRegister::CG;
};

// ===== IR =====

template < class ResidualOperator_t, class SetupType, class SmootherOperator_t >
struct SolverTrait< SolverRegister::IR, ResidualOperator_t, SetupType, SmootherOperator_t >
{
   typedef IterativeRefinementSolver< SmootherOperator_t, ResidualOperator_t > Solver_t;
   static constexpr auto                                                       toString = "IR";

   static std::shared_ptr< Solver_t > getInitializedSolver( ResidualOperator_t&,
                                                            const SetupType&    config,
                                                            Storage&            storage,
                                                            uint_t              minLevel,
                                                            uint_t              maxLevel,
                                                            SmootherOperator_t& dotA )
   {
      // +++ Setup of inner Solver +++
      auto innerSolver = SolverTrait< SetupType::IRSmoother, SmootherOperator_t, SetupType >::getInitializedSolver(
          dotA, config, storage, minLevel, maxLevel );

      auto solver =
          std::make_shared< Solver_t >( storage, dotA, innerSolver, minLevel, maxLevel, config.numerOfInnerSolveCyclesIR );

      return solver;
   }

 private:
};

// +++ Solvers that can't be used for elementwise operators +++

// ===== WeightedJacobiSmoother =====

template < class Operator_t, class SetupType >
struct SolverTrait< SolverRegister::WeightedJacobiSmoother, Operator_t, SetupType >
{
   typedef WeightedJacobiSmoother< Operator_t > Solver_t;
   static constexpr auto                        toString = "WeightedJacobiSmoother";

   static std::shared_ptr< Solver_t > getInitializedSolver( Operator_t& A,
                                                            const SetupType&,
                                                            Storage& storage,
                                                            uint_t   minLevel,
                                                            uint_t   maxLevel,
                                                            std::optional< std::reference_wrapper< Operator_t > > = std::nullopt )
   {
      A.computeInverseDiagonalOperatorValues();
      return std::make_shared< Solver_t >( storage, minLevel, maxLevel, jacobiRelaxationFactor );
   }

 private:
   static constexpr real_t jacobiRelaxationFactor =
       real_t( 2.0 / 3 ); // according to Wikipedia 2/3 is the usual choice for weighted Jacobi
};

} // namespace hyteg::mixedPrecisionMT
