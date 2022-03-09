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

#pragma once

#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

namespace hyteg {

/// \brief Contains convenience functions that assemble solvers for common problems.
namespace solvertemplates {

/// \brief Returns a pressure preconditioned MINRES solver for the Stokes system.
///
/// The pressure is pre-multiplied with the inverse of the lumped mass matrix.
/// It is assumed that the pressure is discretized with P1 finite elements.
///
/// \tparam StokesOperatorType most types of Stokes operators should be possible to pass here
///                            since the MINRES solver does not require any knowledge about the
///                            structure of the A-block
/// \param storage the PrimitiveStorage that defines the domain
/// \param level the refinement level of the grid
/// \param absoluteTargetResidual absolute (as opposed to relative) residual as a stopping criterion for the iteration
/// \param maxIterations if not converged to the target residual, the iteration stops after this many iterations
///
template < typename StokesOperatorType >
std::shared_ptr< Solver< StokesOperatorType > > stokesMinResSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    const uint_t&                              level,
                                                                    const real_t& absoluteTargetResidual,
                                                                    const uint_t& maxIterations,
                                                                    bool printInfo = false )
{
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< StokesOperatorType, P1LumpedInvMassOperator > >(
           storage, level, level );
   auto pressurePreconditionedMinResSolver = std::make_shared< MinResSolver< StokesOperatorType > >(
       storage, level, level, maxIterations, absoluteTargetResidual, pressurePreconditioner );

   pressurePreconditionedMinResSolver->setPrintInfo( printInfo );

   return pressurePreconditionedMinResSolver;
}

/// \brief Returns a geometric multigrid solver for the constant-coefficient Stokes system.
///
/// This solver performs a v-cycle and employs as smoother the inexact Uzawa method.
/// The relaxation on the A-block is performed with a forward Gauss-Seidel.
/// As coarse grid solver a pressure-preconditioned MINRES is employed.
///
/// As for all solver templates, this should not be used if performance
/// is critical. This serves as an alternative to speed up the implementation,
/// not the solution of a Stokes system solver. For fine grained control
/// the implementation can serve as a starting point.
///
/// \param storage the PrimitiveStorage that defines the domain
/// \param minLevel the coarse grid solver level
/// \param maxLevel the finest level
/// \param preSmoothingSteps number of relaxation iterations performed with the inexact Uzawa iteration
///                          before coarsening
/// \param postSmoothingSteps number of relaxation iterations performed with the inexact Uzawa iteration
///                           before after refinement
/// \param uzawaSmootherOmega relaxation parameter for the pressure preconditioner in the inexact Uzawa smoother,
///                           can be estimated by solving a related eigenvalue problem,
///                           usually 0 < omega < 1
///
template < typename StokesOperatorType >
std::shared_ptr< Solver< StokesOperatorType > > stokesGMGUzawaSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                      const uint_t&                              minLevel,
                                                                      const uint_t&                              maxLevel,
                                                                      const uint_t& preSmoothingSteps,
                                                                      const uint_t& postSmoothingSteps,
                                                                      const real_t& uzawaSmootherOmega )
{
   WALBERLA_ABORT( "Uzawa GMG solver template is not implemented for this operator." );
}

template <>
std::shared_ptr< Solver< P1P1StokesOperator > >
    stokesGMGUzawaSolver< P1P1StokesOperator >( const std::shared_ptr< PrimitiveStorage >& storage,
                                              const uint_t&                              minLevel,
                                              const uint_t&                              maxLevel,
                                              const uint_t&                              preSmoothingSteps,
                                              const uint_t&                              postSmoothingSteps,
                                              const real_t&                              uzawaSmootherOmega )
{
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< P1P1StokesOperator, P1LumpedInvMassOperator > >(
           storage, minLevel, minLevel );
   auto pressurePreconditionedMinResSolver =
       std::make_shared< MinResSolver< P1P1StokesOperator > >( storage, minLevel, minLevel, 1000, 1e-12, pressurePreconditioner );

   auto stokesRestriction  = std::make_shared< P1P1StokesToP1P1StokesRestriction >();
   auto stokesProlongation = std::make_shared< P1P1StokesToP1P1StokesProlongation >();
   auto gaussSeidel        = std::make_shared< GaussSeidelSmoother< P1P1StokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< P1P1StokesOperator > >( storage, gaussSeidel );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< P1P1StokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaSmootherOmega );

   auto gmgSolver = std::make_shared< GeometricMultigridSolver< P1P1StokesOperator > >( storage,
                                                                                      uzawaSmoother,
                                                                                      pressurePreconditionedMinResSolver,
                                                                                      stokesRestriction,
                                                                                      stokesProlongation,
                                                                                      minLevel,
                                                                                      maxLevel,
                                                                                      preSmoothingSteps,
                                                                                      postSmoothingSteps,
                                                                                      2 );
   return gmgSolver;
}

template <>
std::shared_ptr< Solver< P2P1TaylorHoodStokesOperator > >
    stokesGMGUzawaSolver< P2P1TaylorHoodStokesOperator >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                          const uint_t&                              minLevel,
                                                          const uint_t&                              maxLevel,
                                                          const uint_t&                              preSmoothingSteps,
                                                          const uint_t&                              postSmoothingSteps,
                                                          const real_t&                              uzawaSmootherOmega )
{
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< P2P1TaylorHoodStokesOperator, P1LumpedInvMassOperator > >(
           storage, minLevel, minLevel );
   auto pressurePreconditionedMinResSolver = std::make_shared< MinResSolver< P2P1TaylorHoodStokesOperator > >(
       storage, minLevel, minLevel, 1000, 1e-12, pressurePreconditioner );

   auto stokesRestriction  = std::make_shared< P2P1StokesToP2P1StokesRestriction >();
   auto stokesProlongation = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto gaussSeidel        = std::make_shared< GaussSeidelSmoother< P2P1TaylorHoodStokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< P2P1TaylorHoodStokesOperator > >( storage, gaussSeidel );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< P2P1TaylorHoodStokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaSmootherOmega );

   auto gmgSolver =
       std::make_shared< GeometricMultigridSolver< P2P1TaylorHoodStokesOperator > >( storage,
                                                                                     uzawaSmoother,
                                                                                     pressurePreconditionedMinResSolver,
                                                                                     stokesRestriction,
                                                                                     stokesProlongation,
                                                                                     minLevel,
                                                                                     maxLevel,
                                                                                     preSmoothingSteps,
                                                                                     postSmoothingSteps,
                                                                                     2 );
   return gmgSolver;
}

} // namespace solvertemplates
} // namespace hyteg
