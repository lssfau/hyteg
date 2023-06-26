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

#include <algorithm>
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
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

namespace hyteg {
namespace moc_benchmarks {

struct LoadBalancingOptions
{
   // 0: round robin
   // 1: greedy
   // 2: round robin volume
   uint_t type = 0;
};

/// Calculates and returns
///
///     ||u||_L2 = sqrt( u^T M u )
///
template < typename FunctionType, typename MassOperator >
inline real_t
    normL2( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, const uint_t& level, const DoFType& flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}

/// Calculates and returns difference in point-wise max peaks:
///
///     (maxU / maxUSolution) - 1
///
template < typename FunctionType >
inline real_t maxPeakDifference( const FunctionType& u, const FunctionType& uSolution, uint_t level, DoFType flag )
{
   const real_t maxTempApproximate = u.getMaxValue( level, flag );
   const real_t maxTempAnalytical  = uSolution.getMaxValue( level, flag );
   const real_t peakError          = maxTempApproximate / maxTempAnalytical - 1;
   return peakError;
}

/// Calculates and returns var(t) as an indicator for the amount of spurious oscillations:
///
///    var(t) := max(u_h) - min(u_h)
///
template < typename FunctionType >
inline real_t spuriousOscillations( const FunctionType& u, uint_t level, DoFType flag )
{
   const real_t maxTemp = u.getMaxValue( level, flag );
   const real_t minTemp = u.getMinValue( level, flag );
   return maxTemp - minTemp;
}

template < typename FunctionType, typename MassOperator >
inline real_t
    globalMass( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, const uint_t& level, const DoFType& flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   auto mass = tmp.sumGlobal( level );
   return mass;
}

class Solution
{
 public:
   Solution()
   : currentTime_( 0 )
   {}

   explicit Solution( real_t t0 )
   : currentTime_( t0 )
   {}

   /// Evaluates the solution at a specific point.
   virtual real_t operator()( const Point3D& ) const = 0;

   /// Increments the current time by dt.
   void incTime( real_t dt ) { currentTime_ += dt; };

   void setTime( real_t currentTime ) { currentTime_ = currentTime; }

   real_t currentTime() const { return currentTime_; }

 protected:
   real_t currentTime_;
};

class ZeroSolution : public Solution
{
 public:
   ZeroSolution()
   : Solution( 0 )
   {}
   real_t operator()( const Point3D& ) const override { return 0; }
};

template < typename FunctionType_T, typename LaplaceOperator_T, typename MassOperator_T, typename UnsteadyDiffusionOperator_T >
void solve( MeshInfo&               meshInfo,
            bool                    setBlendingMap,
            Solution&               solution,
            Solution&               velocityX,
            Solution&               velocityY,
            Solution&               velocityZ,
            real_t                  dt,
            real_t                  diffusivity,
            uint_t                  level,
            DiffusionTimeIntegrator diffusionTimeIntegrator,
            bool                    enableDiffusion,
            bool                    strangSplitting,
            bool                    resetParticles,
            uint_t                  resetParticlesInterval,
            bool                    adjustedAdvection,
            bool                    globalMaxLimiter,
            bool                    setParticlesOutsideDomainToZero,
            uint_t                  numTimeSteps,
            LoadBalancingOptions    lbOptions,
            bool                    vtk,
            bool                    vtkOutputVelocity,
            const std::string&      benchmarkName,
            uint_t                  printInterval,
            uint_t                  vtkInterval,
            bool                    verbose,
            std::string             dbFile );

} // namespace moc_benchmarks
} // namespace hyteg
