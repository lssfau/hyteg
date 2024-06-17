/*
 * Copyright (c) 2024 Andreas Burkhart, Ponsuganth Ilangovan P.
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

#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/FullStokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

namespace hyteg {
namespace solvertemplates {

// clang-format off

enum class StokesGMGUzawaFSSolverParamKey
{
    NUM_POWER_ITERATIONS_SPECTRUM,
    NUM_COARSE_GRID_ITERATIONS,
    COARSE_GRID_TOLERANCE,
    UZAWA_OMEGA,
    MG_PRE_SMOOTH,
    MG_POST_SMOOTH,
    UZAWA_VELOCITY_ITER,
    SMOOTH_INCREMENT_COARSE_GRID
};

// clang-format on

// Things are still restricted to P2-P1 space
template < typename StokesOperatorType, typename ProjectionType >
inline std::shared_ptr< Solver< StokesOperatorType > >
    stokesGMGUzawaFSSolver( const std::shared_ptr< PrimitiveStorage >&         storage,
                            const uint_t&                                      minLevel,
                            const uint_t&                                      maxLevel,
                            const std::shared_ptr< StokesOperatorType >&       stokesOperatorFSSelf,
                            const std::shared_ptr< ProjectionType >&           projectionOperator,
                            BoundaryCondition                                  bcVelocity,
                            bool                                               verbose     = false,
                            std::map< StokesGMGUzawaFSSolverParamKey, real_t > extraParams = {} )
{
   std::map< StokesGMGUzawaFSSolverParamKey, real_t > defaultParams = {
       { StokesGMGUzawaFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM, 25.0 },
       { StokesGMGUzawaFSSolverParamKey::NUM_COARSE_GRID_ITERATIONS, 10 },
       { StokesGMGUzawaFSSolverParamKey::COARSE_GRID_TOLERANCE, 1e-6 },
       { StokesGMGUzawaFSSolverParamKey::UZAWA_OMEGA, 0.3 },
       { StokesGMGUzawaFSSolverParamKey::MG_PRE_SMOOTH, 3 },
       { StokesGMGUzawaFSSolverParamKey::MG_POST_SMOOTH, 3 },
       { StokesGMGUzawaFSSolverParamKey::UZAWA_VELOCITY_ITER, 1 },
       { StokesGMGUzawaFSSolverParamKey::SMOOTH_INCREMENT_COARSE_GRID, 2 } };

   for ( auto const& param : extraParams )
   {
      if ( defaultParams.find( param.first ) == defaultParams.end() )
      {
         WALBERLA_ABORT( "Passed parameter key not found" );
      }
      else
      {
         defaultParams[param.first] = real_c( param.second );
      }
   }

   uint_t numPowerIteration = uint_c( defaultParams[StokesGMGUzawaFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM] );

   uint_t coarseGridIter = uint_c( defaultParams[StokesGMGUzawaFSSolverParamKey::NUM_COARSE_GRID_ITERATIONS] );
   real_t coarseGridTol  = real_c( defaultParams[StokesGMGUzawaFSSolverParamKey::COARSE_GRID_TOLERANCE] );

   real_t uzawaSmootherOmega = real_c( defaultParams[StokesGMGUzawaFSSolverParamKey::UZAWA_OMEGA] );

   uint_t MGPreSmooth  = uint_c( defaultParams[StokesGMGUzawaFSSolverParamKey::MG_PRE_SMOOTH] );
   uint_t MGPostSmooth = uint_c( defaultParams[StokesGMGUzawaFSSolverParamKey::MG_POST_SMOOTH] );

   uint_t uzawaVelocityIter         = uint_c( defaultParams[StokesGMGUzawaFSSolverParamKey::UZAWA_VELOCITY_ITER] );
   uint_t smoothIncrementCoarseGrid = uint_c( defaultParams[StokesGMGUzawaFSSolverParamKey::SMOOTH_INCREMENT_COARSE_GRID] );

   uint_t powerIterations = numPowerIteration;

   auto smoother = std::make_shared< ChebyshevSmootherWithFreeSlipProjection< typename StokesOperatorType::VelocityOperator_T > >(
       storage, minLevel, maxLevel, projectionOperator );

   auto uSpec = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
       "uSpec_stokesGMGUzawaFSSolver_solvertemplate", storage, minLevel, maxLevel, bcVelocity );

   auto uTmpSpec = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
       "uTmpSpec_stokesGMGUzawaFSSolver_solvertemplate", storage, minLevel, maxLevel, bcVelocity );

   std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   WALBERLA_LOG_INFO_ON_ROOT( "Estimate spectral radius!" );
   // avoid that the startpoint of our poweriteration is in the kernel of the operator
   uSpec->uvw().interpolate( { randFuncA, randFuncA, randFuncA }, maxLevel, Inner );

   real_t spectralRadius = chebyshev::estimateRadius(
       stokesOperatorFSSelf->getA().viscousOperator, maxLevel, powerIterations, storage, uSpec->uvw(), uTmpSpec->uvw() );

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated spectral radius: " << spectralRadius );

   smoother->setupCoefficients( 3, spectralRadius );

   auto _coarseGridSolver = solvertemplates::stokesMinResSolver< StokesOperatorType >( storage,
                                                                                   minLevel,
                                                                                   coarseGridTol,
                                                                                   coarseGridIter,
                                                                                   true );

   auto uzawaVelocityPreconditioner =
       std::make_shared< FullStokesVelocityBlockBlockDiagonalPreconditioner< StokesOperatorType > >( storage, smoother );

   auto _UzawaSmoother =
       std::make_shared< UzawaSmootherWithFreeSlipProjection< StokesOperatorType > >( storage,
                                                                              uzawaVelocityPreconditioner,
                                                                              minLevel,
                                                                              maxLevel,
                                                                              uzawaSmootherOmega,
                                                                              Inner | NeumannBoundary | FreeslipBoundary,
                                                                              projectionOperator,
                                                                              uzawaVelocityIter );

   auto prolongationOperator =
       std::make_shared< P2P1StokesToP2P1StokesProlongationWithFreeSlipProjection >( uTmpSpec, projectionOperator );
   auto restrictionOperator = std::make_shared< P2P1StokesToP2P1StokesRestrictionWithFreeSlipProjection >( projectionOperator );

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperatorType > >( storage,
                                                                                              _UzawaSmoother,
                                                                                              _coarseGridSolver,
                                                                                              restrictionOperator,
                                                                                              prolongationOperator,
                                                                                              minLevel,
                                                                                              maxLevel,
                                                                                              MGPreSmooth,
                                                                                              MGPostSmooth,
                                                                                              smoothIncrementCoarseGrid,
                                                                                              CycleType::VCYCLE );

   return multigridSolver;
}

} // namespace solvertemplates
} // namespace hyteg