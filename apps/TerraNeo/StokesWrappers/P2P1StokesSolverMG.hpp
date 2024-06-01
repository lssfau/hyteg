/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P, Andreas Burkhart.
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

#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"

#include "P2P1StokesOperatorProjection.hpp"

namespace hyteg {

inline std::shared_ptr< Solver< P2P1StokesFullIcosahedralShellMapOperatorFS > >
    stokesGMGFSSolver( const std::shared_ptr< PrimitiveStorage >&                            storage,
                       const uint_t&                                                         minLevel,
                       const uint_t&                                                         maxLevel,
                       const std::shared_ptr< P2P1StokesFullIcosahedralShellMapOperatorFS >& stokesOperatorFSSelf,
                       const std::shared_ptr< P2ProjectNormalOperator >&                     projectionOperator,
                       const uint_t&                                                         fGMRESOuterIter,
                       const real_t&                                                         uzawaSmootherOmega,
                       BoundaryCondition                                                     bcVelocity,
                       bool                                                                  estimateUzawaOmegaValue = false,
                       bool                                                                  verbose                 = false )
{
   auto uTmp = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
       "uTmp_stokesGMGFSSolver_solvertemplate", storage, minLevel, maxLevel, bcVelocity );

   auto tempFct = std::make_shared< P2VectorFunction< real_t > >(
       "tempFct_stokesGMGFSSolver_solvertemplate", storage, minLevel, maxLevel, bcVelocity );

   auto uSpec = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
       "uSpec_stokesGMGFSSolver_solvertemplate", storage, minLevel, maxLevel, bcVelocity );

   auto uTmpSpec = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
       "uTmpSpec_stokesGMGFSSolver_solvertemplate", storage, minLevel, maxLevel, bcVelocity );

   uint_t numPowerIteration = 25U;

   auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongationWithProjection >( uTmp, projectionOperator );
   auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestrictionWithProjection >( projectionOperator );

   // Multigridsolver for A
   typedef P2P1StokesFullIcosahedralShellMapOperatorFS::ViscousOperatorFS_T SubstAType;
   typedef P2P1StokesFullIcosahedralShellMapOperatorFS                      StokesOperatorFS;

   auto APrecOperator = stokesOperatorFSSelf->getA();

   APrecOperator.computeInverseDiagonalOperatorValues();

   auto ABlockProlongationOperator =
       std::make_shared< P2toP2QuadraticVectorProlongationWithProjection >( tempFct, projectionOperator );
   auto ABlockRestrictionOperator = std::make_shared< P2toP2QuadraticVectorRestrictionWithProjection >( projectionOperator );

   // TODO: TO CHECK THIS!
   // auto Jacobi = std::make_shared< WeightedJacobiSmoother< SubstAType > >( storage, minLevel, maxLevel, 2.0 / 3.0 );

   // auto ABlockCoarseGridSolver = std::make_shared< PETScLUSolver< SubstAType > >( storage, minLevel );

   auto ABlockCoarseGridSolver = std::make_shared< hyteg::CGSolver< SubstAType > >( storage, minLevel, maxLevel, 10, 1e-8 );

   // auto ABlockCoarseGridSolver =
   //     std::make_shared< hyteg::MinResSolver< SubstAType > >( storage, minLevel, maxLevel, 100, 1e-8 );
   ABlockCoarseGridSolver->setPrintInfo( verbose );
   // ABlockCoarseGridSolver->setSolverName("ABlockCoarseGridSolver");

   auto ABlockSmoother =
       std::make_shared< ChebyshevSmootherWithProjection< SubstAType > >( storage, minLevel, maxLevel, projectionOperator );

   std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   WALBERLA_LOG_INFO_ON_ROOT( "Estimate spectral radius!" );
   // avoid that the startpoint of our poweriteration is in the kernel of the operator
   uSpec->uvw().interpolate( randFuncA, maxLevel, All );

   auto spectralRadiusA = chebyshev::estimateRadius(
       APrecOperator.viscousOperator, maxLevel, numPowerIteration, storage, uSpec->uvw(), uTmpSpec->uvw() );
   uSpec->uvw().interpolate( 0, maxLevel, All );
   uTmpSpec->interpolate( 0, maxLevel, All );

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated spectral radius: " << spectralRadiusA );

   ABlockSmoother->setupCoefficients( 3, spectralRadiusA );

   auto ABlockMultigridSolver = std::make_shared< GeometricMultigridSolver< SubstAType > >( storage,
                                                                                            ABlockSmoother,
                                                                                            ABlockCoarseGridSolver,
                                                                                            ABlockRestrictionOperator,
                                                                                            ABlockProlongationOperator,
                                                                                            minLevel,
                                                                                            maxLevel,
                                                                                            3,
                                                                                            3,
                                                                                            0,
                                                                                            CycleType::VCYCLE );

   auto ABlockSolver =
       std::make_shared< hyteg::CGSolver< SubstAType > >( storage, minLevel, maxLevel, 3, 1e-6, ABlockMultigridSolver );
   ABlockSolver->setPrintInfo( verbose );
   // ABlockSolver->setSolverName("ABlockSolverOuter");

   // estimate sigma
   real_t estimatedSigma = 1.0;
   if ( estimateUzawaOmegaValue )
   {
      estimatedSigma =
          estimateUzawaSigma< StokesOperatorFS, typename StokesOperatorFS::ViscousOperatorFS_T >( storage,
                                                                                                  minLevel,
                                                                                                  maxLevel,
                                                                                                  *stokesOperatorFSSelf,
                                                                                                  ABlockSolver,
                                                                                                  5,
                                                                                                  1,
                                                                                                  bcVelocity,
                                                                                                  bcVelocity,
                                                                                                  bcVelocity,
                                                                                                  42,
                                                                                                  projectionOperator );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated sigma: " << estimatedSigma );

   // Multigridsolver for Schur
   typedef StokesOperatorFS::SchurOperator_T SubstSType;

   auto SchurProlongationOperator = std::make_shared< P1toP1LinearProlongation< real_t > >();
   auto SchurRestrictionOperator  = std::make_shared< P1toP1LinearRestriction< real_t > >();

   auto SchurCoarseGridSolver = std::make_shared< hyteg::CGSolver< SubstSType > >( storage, minLevel, maxLevel, 1 );
   SchurCoarseGridSolver->setPrintInfo( verbose );
   // SchurCoarseGridSolver->setSolverName("SchurCoarseGridSolver");

   auto SchurSmoother = std::make_shared< ChebyshevSmoother< SubstSType > >( storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > randFuncS = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   // avoid that the startpoint of our poweriteration is in the kernel of the operator
   uSpec->p().interpolate( randFuncS, maxLevel, All );
   auto spectralRadiusSchur = chebyshev::estimateRadius(
       stokesOperatorFSSelf->getSchur(), maxLevel, numPowerIteration, storage, uSpec->p(), uTmpSpec->p() );
   uSpec->p().interpolate( 0, maxLevel, All );
   uTmpSpec->interpolate( 0, maxLevel, All );

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated spectral radius Schur: " << spectralRadiusSchur );

   SchurSmoother->setupCoefficients( 2, spectralRadiusSchur );

   auto SchurMultigridSolver_ = std::make_shared< GeometricMultigridSolver< SubstSType > >( storage,
                                                                                            SchurSmoother,
                                                                                            SchurCoarseGridSolver,
                                                                                            SchurRestrictionOperator,
                                                                                            SchurProlongationOperator,
                                                                                            minLevel,
                                                                                            maxLevel,
                                                                                            3,
                                                                                            3,
                                                                                            0,
                                                                                            CycleType::VCYCLE );

   auto SchurSolver = std::make_shared< CGSolver< SubstSType > >( storage, minLevel, maxLevel, 10, 1e-8, SchurMultigridSolver_ );
   SchurSolver->setPrintInfo( verbose );
   // SchurSolver->setSolverName("SchurSolver");

   real_t estimatedOmega = 1.0;
   if ( estimateUzawaOmegaValue )
   {
      ABlockCoarseGridSolver->setPrintInfo( verbose );
      estimatedOmega =
          estimateUzawaOmega< StokesOperatorFS, StokesOperatorFS::ViscousOperatorFS_T, StokesOperatorFS::SchurOperator_T >(
              storage,
              minLevel,
              maxLevel,
              *stokesOperatorFSSelf,
              stokesOperatorFSSelf->getSchur(),
              ABlockCoarseGridSolver,
              SchurSolver,
              15,
              1,
              bcVelocity,
              bcVelocity,
              bcVelocity,
              42,
              projectionOperator );
      ABlockCoarseGridSolver->setPrintInfo( false );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Estimated omega: " << estimatedOmega );

   WALBERLA_LOG_INFO_ON_ROOT( "Sigma: " << estimatedSigma << " ; "
                                        << "Omega: " << estimatedOmega );

   auto uzawaSmoother = std::make_shared<
       InexactUzawaPreconditioner< StokesOperatorFS, StokesOperatorFS::ViscousOperatorFS_T, StokesOperatorFS::SchurOperator_T > >(
       storage,
       minLevel,
       maxLevel,
       stokesOperatorFSSelf->getSchur(),
       ABlockSolver,
       SchurSolver,
       estimatedSigma,
       estimatedOmega,
       1U,
       projectionOperator );

   auto finalStokesSolver = std::make_shared< FGMRESSolver< StokesOperatorFS > >(
       storage, minLevel, maxLevel, fGMRESOuterIter, 50, 1e-6, 1e-6, 0, uzawaSmoother );
   finalStokesSolver->setPrintInfo( true );

   return finalStokesSolver;
}
} // namespace hyteg