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

#include <cmath>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/Git.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1InjectionRestriction.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/SORSmoother.hpp"
#include "hyteg/solvers/SymmetricSORSmoother.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/controlflow/TimedSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

namespace hyteg {
namespace tme_benchmarks {

/// Calculates and returns
///
///     ||u||_L2 = sqrt( u^T M u )
///
template < template < typename > class FunctionType, typename MassOperator >
inline real_t normL2Scalar( const FunctionType< real_t >& u,
                            const MassOperator&           M,
                            const FunctionType< real_t >& tmp,
                            const uint_t&                 level,
                            const DoFType&                flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   auto dot = u.dotGlobal( tmp, level, flag );
   return std::sqrt( dot );
}

template < template < typename > class FunctionType, typename MassOperator >
inline real_t normL2ScalarSquare( const FunctionType< real_t >& u,
                                  const MassOperator&           M,
                                  const FunctionType< real_t >& tmp,
                                  const uint_t&                 level,
                                  const DoFType&                flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   auto dot = u.dotGlobal( tmp, level, flag );
   return dot;
}

template < template < typename > class StokesFunctionType, typename VelocityMassOperator >
inline real_t normL2Velocity( const StokesFunctionType< real_t >& u,
                              const VelocityMassOperator&         M,
                              const StokesFunctionType< real_t >& tmp,
                              uint_t                              level,
                              DoFType                             flag )
{
   auto normUSquare = normL2ScalarSquare( u.uvw()[0], M, tmp.uvw()[0], level, flag );
   auto normVSquare = normL2ScalarSquare( u.uvw()[1], M, tmp.uvw()[1], level, flag );
   auto normWSquare = normL2ScalarSquare( u.uvw()[2], M, tmp.uvw()[2], level, flag );

   return std::sqrt( normUSquare + normVSquare + normWSquare );
}

template < template < typename > class StokesFunctionType, typename StokesOperator >
inline void residual( const StokesFunctionType< real_t >& u,
                      const StokesFunctionType< real_t >& f,
                      const StokesOperator&               A,
                      const StokesFunctionType< real_t >& tmp,
                      uint_t                              level,
                      DoFType                             flag,
                      StokesFunctionType< real_t >&       r )
{
   A.apply( u, tmp, level, flag );
   r.assign( { 1.0, -1.0 }, { f, tmp }, level, flag );
}

template < template < typename > class StokesFunctionType >
inline void error( const StokesFunctionType< real_t >& u,
                   const StokesFunctionType< real_t >& exact,
                   uint_t                              level,
                   DoFType                             flag,
                   StokesFunctionType< real_t >&       error )
{
   error.assign( { 1.0, -1.0 }, { u, exact }, level, flag );
}

enum class Discretization
{
   P2_P1,
   P1_P1
};

struct MultigridSettings
{
   uint_t preSmooth;
   uint_t postSmooth;
   uint_t incSmooth;

   // 0   : no FMG
   // else: number of inner cycles
   uint_t fmgInnerIterations;

   uint_t numCycles;

   // stops cycling after both, pressure and velocity residual drop below this value
   real_t absoluteResidualTolerance;
};

struct SmootherSettings
{
   uint_t numGSVelocity;
   uint_t symmGSVelocity;

   bool   estimateOmega;
   uint_t omegaEstimationLevel;
   uint_t omegaEstimationIterations;
   real_t omega;
};

struct CoarseGridSettings
{
   // coarse grid solver type:
   // 0: MUMPS                          (PETSc)
   // 1: block preconditioned MINRES    (PETSc)
   uint_t solverType;

   real_t absoluteResidualTolerance;
   uint_t maxIterations;
};

void solve( const std::shared_ptr< PrimitiveStorage >&              storage,
            Discretization                                          discretization,
            bool                                                    hasAnalyticalSolution,
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
            std::string                                             dbFileDiscrError );

} // namespace tme_benchmarks
} // namespace hyteg
