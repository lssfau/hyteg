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
    COARSE_GRID_PETSC,
    UZAWA_OMEGA,
    MG_PRE_SMOOTH,
    MG_POST_SMOOTH,
    UZAWA_VELOCITY_ITER,
    SMOOTH_INCREMENT_COARSE_GRID
};

// clang-format on

// Things are still restricted to P2-P1 space
template < typename StokesOperatorType, typename ProjectionType >
inline std::tuple< std::shared_ptr< Solver< StokesOperatorType > >,
                   std::shared_ptr< Solver< typename StokesOperatorType::ViscousOperatorFS_T > > >
    stokesGMGUzawaFSSolver( const std::shared_ptr< PrimitiveStorage >&                                 storage,
                            const uint_t&                                                              minLevel,
                            const uint_t&                                                              maxLevel,
                            const std::shared_ptr< StokesOperatorType >&                               stokesOperatorFSSelf,
                            const std::shared_ptr< ProjectionType >&                                   projectionOperator,
                            const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&                 tmp1,
                            const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&                 tmpProlongation,
                            bool                                                                       verbose,
                            std::map< StokesGMGUzawaFSSolverParamKey, std::variant< real_t, uint_t > > extraParams )
{
   std::shared_ptr< Solver< StokesOperatorType > >                            _coarseGridSolver;
   std::map< StokesGMGUzawaFSSolverParamKey, std::variant< real_t, uint_t > > defaultParams = {
       { StokesGMGUzawaFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM, uint_c( 25u ) },
       { StokesGMGUzawaFSSolverParamKey::NUM_COARSE_GRID_ITERATIONS, uint_c( 10u ) },
       { StokesGMGUzawaFSSolverParamKey::COARSE_GRID_TOLERANCE, real_c( 1e-6 ) },
       { StokesGMGUzawaFSSolverParamKey::COARSE_GRID_PETSC, uint_c( 0 ) },
       { StokesGMGUzawaFSSolverParamKey::UZAWA_OMEGA, real_c( 0.3 ) },
       { StokesGMGUzawaFSSolverParamKey::MG_PRE_SMOOTH, uint_c( 3u ) },
       { StokesGMGUzawaFSSolverParamKey::MG_POST_SMOOTH, uint_c( 3u ) },
       { StokesGMGUzawaFSSolverParamKey::UZAWA_VELOCITY_ITER, uint_c( 1u ) },
       { StokesGMGUzawaFSSolverParamKey::SMOOTH_INCREMENT_COARSE_GRID, uint_c( 2u ) } };

   for ( auto const& param : extraParams )
   {
      if ( defaultParams.find( param.first ) == defaultParams.end() )
      {
         WALBERLA_ABORT( "Passed parameter key not found" );
      }
      else
      {
         defaultParams[param.first] = param.second;
      }
   }

   uint_t numPowerIteration = std::get< uint_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM] );

   uint_t coarseGridIter        = std::get< uint_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::NUM_COARSE_GRID_ITERATIONS] );
   real_t coarseGridTol         = std::get< real_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::COARSE_GRID_TOLERANCE] );
   uint_t coarseGridSolverPETSc = std::get< uint_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::COARSE_GRID_PETSC] );

   real_t uzawaSmootherOmega = std::get< real_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::UZAWA_OMEGA] );

   uint_t MGPreSmooth  = std::get< uint_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::MG_PRE_SMOOTH] );
   uint_t MGPostSmooth = std::get< uint_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::MG_POST_SMOOTH] );

   uint_t uzawaVelocityIter = std::get< uint_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::UZAWA_VELOCITY_ITER] );
   uint_t smoothIncrementCoarseGrid =
       std::get< uint_t >( defaultParams[StokesGMGUzawaFSSolverParamKey::SMOOTH_INCREMENT_COARSE_GRID] );

   uint_t powerIterations = numPowerIteration;

   auto smoother =
       std::make_shared< ChebyshevSmoother< typename StokesOperatorType::VelocityOperator_T, P2ProjectNormalOperator > >(
           storage, minLevel, maxLevel, false, projectionOperator, FreeslipBoundary );

   std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   WALBERLA_LOG_INFO_ON_ROOT( "Estimate spectral radius!" );
   // avoid that the startpoint of our poweriteration is in the kernel of the operator
   tmp1->uvw().interpolate( { randFuncA, randFuncA, randFuncA }, maxLevel, Inner );

   real_t spectralRadius = chebyshev::estimateRadius(
       stokesOperatorFSSelf->getA(), maxLevel, powerIterations, storage, tmp1->uvw(), tmpProlongation->uvw() );

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated spectral radius: " << spectralRadius );

   smoother->setupCoefficients( 3, spectralRadius );

   if ( coarseGridSolverPETSc == 1 )
   {
#ifdef HYTEG_BUILD_WITH_PETSC
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: PETSc MinRes" );
      PETScManager petscManager;
      auto         coarseGridrelativeTolerance = coarseGridTol;
      auto         PETScCoarseGridSolver       = std::make_shared< PETScMinResSolver< StokesOperatorType > >(
          storage, minLevel, coarseGridIter, coarseGridrelativeTolerance, coarseGridTol );
      PETScCoarseGridSolver->reassembleMatrix( true );
      _coarseGridSolver = PETScCoarseGridSolver;
#else
      WALBERLA_LOG_INFO_ON_ROOT( "PETSc module not found or enabled. Switching to HyTeG MinRes." );
      // Fall back to HyTeG MinRes solver
#endif
   }
   if ( !coarseGridSolverPETSc || !_coarseGridSolver )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solver: HyTeG MinRes" );
      _coarseGridSolver =
          solvertemplates::stokesMinResSolver< StokesOperatorType >( storage, minLevel, coarseGridTol, coarseGridIter, verbose );
   }

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
       std::make_shared< P2P1StokesToP2P1StokesProlongationWithFreeSlipProjection >( tmpProlongation, projectionOperator );
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

   return { multigridSolver, smoother };
}

} // namespace solvertemplates
} // namespace hyteg
