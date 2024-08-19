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
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"

namespace hyteg {
namespace solvertemplates {

// clang-format off

enum class StokesGMGFSSolverParamKey
{
    NUM_POWER_ITERATIONS_SPECTRUM,
    FGMRES_UZAWA_PRECONDITIONED_OUTER_ITER,
    FGMRES_UZAWA_PRECONDITIONED_OUTER_TOLERANCE,   

        INEXACT_UZAWA_VELOCITY_ITER, 
        INEXACT_UZAWA_OMEGA,

            ABLOCK_CG_SOLVER_MG_PRECONDITIONED_ITER,
            ABLOCK_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE,
                ABLOCK_MG_PRESMOOTH,    
                ABLOCK_MG_POSTSMOOTH,
                    ABLOCK_COARSE_ITER,
                    ABLOCK_COARSE_TOLERANCE,
            
            SCHUR_CG_SOLVER_MG_PRECONDITIONED_ITER,
            SCHUR_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE,
                SCHUR_MG_PRESMOOTH,    
                SCHUR_MG_POSTSMOOTH,
                    SCHUR_COARSE_GRID_CG_ITER,
                    SCHUR_COARSE_GRID_CG_TOLERANCE
};

// clang-format on

// Things are still restricted to P2-P1 space
template < typename StokesOperatorType, typename ProjectionType >
inline std::tuple< std::shared_ptr< Solver< StokesOperatorType > >,
                   std::shared_ptr< Solver< typename StokesOperatorType::ViscousOperatorFS_T > > >
    stokesGMGFSSolver( const std::shared_ptr< PrimitiveStorage >&                            storage,
                       const uint_t&                                                         minLevel,
                       const uint_t&                                                         maxLevel,
                       const std::shared_ptr< StokesOperatorType >&                          stokesOperatorFSSelf,
                       const std::shared_ptr< ProjectionType >&                              projectionOperator,
                       const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&            temp1,
                       const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&            temp2,
                       const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&            temp3,
                       bool                                                                  estimateUzawaOmegaValue = false,
                       bool                                                                  verbose                 = false,
                       std::map< StokesGMGFSSolverParamKey, std::variant< real_t, uint_t > > extraParams             = {} )
{
   std::map< StokesGMGFSSolverParamKey, std::variant< real_t, uint_t > > defaultParams = {
       { StokesGMGFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM, 25u },
       { StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_ITER, 5u },
       { StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_TOLERANCE, 1e-6 },
       { StokesGMGFSSolverParamKey::INEXACT_UZAWA_VELOCITY_ITER, 1u },
       { StokesGMGFSSolverParamKey::INEXACT_UZAWA_OMEGA, 1.0 },
       { StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_ITER, 3u },
       { StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE, 1e-6 },
       { StokesGMGFSSolverParamKey::ABLOCK_MG_PRESMOOTH, 3u },
       { StokesGMGFSSolverParamKey::ABLOCK_MG_POSTSMOOTH, 3u },
       { StokesGMGFSSolverParamKey::ABLOCK_COARSE_ITER, 10u },
       { StokesGMGFSSolverParamKey::ABLOCK_COARSE_TOLERANCE, 1e-8 },
       { StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_ITER, 1u },
       { StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE, 1e-6 },
       { StokesGMGFSSolverParamKey::SCHUR_MG_PRESMOOTH, 3u },
       { StokesGMGFSSolverParamKey::SCHUR_MG_POSTSMOOTH, 3u },
       { StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_ITER, 1u },
       { StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_TOLERANCE, 1e-6 },
   };

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

   uint_t numPowerIteration = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM] );

   uint_t fGMRESOuterIter =
       std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_ITER] );
   real_t fGMRESTol = std::get< real_t >( defaultParams[StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_TOLERANCE] );

   uint_t uzawaVelocityIter  = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::INEXACT_UZAWA_VELOCITY_ITER] );
   real_t uzawaSmootherOmega = std::get< real_t >( defaultParams[StokesGMGFSSolverParamKey::INEXACT_UZAWA_OMEGA] );

   uint_t ABlockCGOuterIter =
       std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_ITER] );
   real_t ABlockCGOuterTol =
       std::get< real_t >( defaultParams[StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE] );

   uint_t ABlockCGPreSmooth  = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::ABLOCK_MG_PRESMOOTH] );
   uint_t ABlockCGPostSmooth = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::ABLOCK_MG_POSTSMOOTH] );

   uint_t ABlockCGCoarseIter = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::ABLOCK_COARSE_ITER] );
   real_t ABlockCGCoarseTol  = std::get< real_t >( defaultParams[StokesGMGFSSolverParamKey::ABLOCK_COARSE_TOLERANCE] );

   uint_t SchurCGOuterIter =
       std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_ITER] );
   real_t SchurCGOuterTol =
       std::get< real_t >( defaultParams[StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE] );

   uint_t SchurCGPreSmooth  = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::SCHUR_MG_PRESMOOTH] );
   uint_t SchurCGPostSmooth = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::SCHUR_MG_POSTSMOOTH] );

   uint_t SchurCGCoarseIter = std::get< uint_t >( defaultParams[StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_ITER] );
   real_t SchurCGCoarseTol  = std::get< real_t >( defaultParams[StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_TOLERANCE] );

   auto prolongationOperator =
       std::make_shared< P2P1StokesToP2P1StokesProlongationWithFreeSlipProjection >( temp1, projectionOperator );
   auto restrictionOperator = std::make_shared< P2P1StokesToP2P1StokesRestrictionWithFreeSlipProjection >( projectionOperator );

   // Multigridsolver for A
   typedef typename StokesOperatorType::ViscousOperatorFS_T SubstAType;
   typedef StokesOperatorType                               StokesOperatorFS;

   auto APrecOperator = stokesOperatorFSSelf->getA();

   APrecOperator.computeInverseDiagonalOperatorValues();

   auto ABlockProlongationOperator = std::make_shared< P2toP2QuadraticVectorProlongationWithFreeSlipProjection >(
       std::make_shared< P2VectorFunction< real_t > >( temp1->uvw() ), projectionOperator );
   auto ABlockRestrictionOperator =
       std::make_shared< P2toP2QuadraticVectorRestrictionWithFreeSlipProjection >( projectionOperator );

   // TODO: TO CHECK THIS!
   // auto Jacobi = std::make_shared< WeightedJacobiSmoother< SubstAType > >( storage, minLevel, maxLevel, 2.0 / 3.0 );

   // auto ABlockCoarseGridSolver = std::make_shared< PETScLUSolver< SubstAType > >( storage, minLevel );

   auto ABlockCoarseGridSolver = std::make_shared< hyteg::MinResSolver< SubstAType > >(
       storage, minLevel, maxLevel, ABlockCGCoarseIter, ABlockCGCoarseTol );

   // auto ABlockCoarseGridSolver =
   //     std::make_shared< hyteg::MinResSolver< SubstAType > >( storage, minLevel, maxLevel, 100, 1e-8 );
   ABlockCoarseGridSolver->setPrintInfo( verbose );
   // ABlockCoarseGridSolver->setSolverName("ABlockCoarseGridSolver");

   auto ABlockSmoother = std::make_shared< ChebyshevSmootherWithFreeSlipProjection< SubstAType > >(
       storage, minLevel, maxLevel, projectionOperator );

   std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   WALBERLA_LOG_INFO_ON_ROOT( "Estimate spectral radius!" );
   // avoid that the startpoint of our poweriteration is in the kernel of the operator
   temp2->uvw().interpolate( randFuncA, maxLevel, All );

   auto spectralRadiusA = chebyshev::estimateRadius(
       APrecOperator.viscousOperator, maxLevel, numPowerIteration, storage, temp2->uvw(), temp3->uvw() );
   temp2->uvw().interpolate( 0, maxLevel, All );
   temp3->interpolate( 0, maxLevel, All );

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated spectral radius: " << spectralRadiusA );

   ABlockSmoother->setupCoefficients( 3, spectralRadiusA );

   auto ABlockMultigridSolver = std::make_shared< GeometricMultigridSolver< SubstAType > >( storage,
                                                                                            ABlockSmoother,
                                                                                            ABlockCoarseGridSolver,
                                                                                            ABlockRestrictionOperator,
                                                                                            ABlockProlongationOperator,
                                                                                            minLevel,
                                                                                            maxLevel,
                                                                                            ABlockCGPreSmooth,
                                                                                            ABlockCGPostSmooth,
                                                                                            0,
                                                                                            CycleType::VCYCLE );

   auto ABlockSolver = std::make_shared< hyteg::CGSolver< SubstAType > >(
       storage, minLevel, maxLevel, ABlockCGOuterIter, ABlockCGOuterTol, ABlockMultigridSolver );
   ABlockSolver->setPrintInfo( verbose );
   // ABlockSolver->setSolverName("ABlockSolverOuter");

   // estimate sigma
   real_t estimatedSigma = uzawaSmootherOmega;
   if ( estimateUzawaOmegaValue )
   {
      WALBERLA_ABORT( "Uzawa estimation not available" );
      //   estimatedSigma = estimateUzawaSigma< StokesOperatorFS, SubstAType >( storage,
      //                                                                        minLevel,
      //                                                                        maxLevel,
      //                                                                        *stokesOperatorFSSelf,
      //                                                                        ABlockSolver,
      //                                                                        5,
      //                                                                        1,
      //                                                                        bcVelocity,
      //                                                                        bcVelocity,
      //                                                                        bcVelocity,
      //                                                                        42,
      //                                                                        projectionOperator );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated sigma: " << estimatedSigma );

   // Multigridsolver for Schur
   typedef typename StokesOperatorFS::SchurOperator_T SubstSType;

   auto SchurProlongationOperator = std::make_shared< P1toP1LinearProlongation< real_t > >();
   auto SchurRestrictionOperator  = std::make_shared< P1toP1LinearRestriction< real_t > >();

   auto SchurCoarseGridSolver =
       std::make_shared< hyteg::CGSolver< SubstSType > >( storage, minLevel, maxLevel, SchurCGCoarseIter, SchurCGCoarseTol );
   SchurCoarseGridSolver->setPrintInfo( verbose );
   // SchurCoarseGridSolver->setSolverName("SchurCoarseGridSolver");

   auto SchurSmoother = std::make_shared< ChebyshevSmoother< SubstSType > >( storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > randFuncS = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   // avoid that the startpoint of our poweriteration is in the kernel of the operator
   temp2->p().interpolate( randFuncS, maxLevel, All );
   auto spectralRadiusSchur = chebyshev::estimateRadius(
       stokesOperatorFSSelf->getSchur(), maxLevel, numPowerIteration, storage, temp2->p(), temp3->p() );
   temp2->p().interpolate( 0, maxLevel, All );
   temp3->interpolate( 0, maxLevel, All );

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated spectral radius Schur: " << spectralRadiusSchur );

   SchurSmoother->setupCoefficients( 2, spectralRadiusSchur );

   auto SchurMultigridSolver_ = std::make_shared< GeometricMultigridSolver< SubstSType > >( storage,
                                                                                            SchurSmoother,
                                                                                            SchurCoarseGridSolver,
                                                                                            SchurRestrictionOperator,
                                                                                            SchurProlongationOperator,
                                                                                            minLevel,
                                                                                            maxLevel,
                                                                                            SchurCGPreSmooth,
                                                                                            SchurCGPostSmooth,
                                                                                            0,
                                                                                            CycleType::VCYCLE );

   auto SchurSolver = std::make_shared< CGSolver< SubstSType > >(
       storage, minLevel, maxLevel, SchurCGOuterIter, SchurCGOuterTol, SchurMultigridSolver_ );
   SchurSolver->setPrintInfo( verbose );
   // SchurSolver->setSolverName("SchurSolver");

   real_t estimatedOmega = uzawaSmootherOmega;
   if ( estimateUzawaOmegaValue )
   {
      WALBERLA_ABORT( "Uzawa estimation not available" );
      //   ABlockCoarseGridSolver->setPrintInfo( verbose );
      //   estimatedOmega = estimateUzawaOmega< StokesOperatorFS, SubstAType, SubstSType >( storage,
      //                                                                                    minLevel,
      //                                                                                    maxLevel,
      //                                                                                    *stokesOperatorFSSelf,
      //                                                                                    stokesOperatorFSSelf->getSchur(),
      //                                                                                    ABlockCoarseGridSolver,
      //                                                                                    SchurSolver,
      //                                                                                    15,
      //                                                                                    1,
      //                                                                                    bcVelocity,
      //                                                                                    bcVelocity,
      //                                                                                    bcVelocity,
      //                                                                                    42,
      //                                                                                    projectionOperator );
      //   ABlockCoarseGridSolver->setPrintInfo( false );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated omega: " << estimatedOmega );

   WALBERLA_LOG_INFO_ON_ROOT( "Sigma: " << estimatedSigma << " ; "
                                        << "Omega: " << estimatedOmega );

   auto uzawaSmoother = std::make_shared< InexactUzawaPreconditioner< StokesOperatorFS, SubstAType, SubstSType > >(
       storage,
       minLevel,
       maxLevel,
       stokesOperatorFSSelf->getSchur(),
       ABlockSolver,
       SchurSolver,
       estimatedSigma,
       estimatedOmega,
       uzawaVelocityIter,
       projectionOperator );

   auto finalStokesSolver = std::make_shared< FGMRESSolver< StokesOperatorFS > >(
       storage, minLevel, maxLevel, fGMRESOuterIter, 50, fGMRESTol, fGMRESTol, 0, uzawaSmoother );
   finalStokesSolver->setPrintInfo( true );

   return { finalStokesSolver, ABlockSmoother };
}

} // namespace solvertemplates
} // namespace hyteg