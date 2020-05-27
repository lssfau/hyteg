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

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/MMOCTransport.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

namespace hyteg {
namespace moc_benchmarks {

/// Calculates and returns
///
///     ||u||_L2 = sqrt( u^T M u )
///
template < typename MassOperator >
inline real_t normL2( const P2Function< real_t >& u,
               const P2Function< real_t >& tmp,
               const MassOperator&         M,
               const uint_t&               level,
               const DoFType&              flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}

/// Calculates and returns difference in point-wise max peaks:
///
///    | maxMagnitudePeakU / maxMagnitudePeakUSolution - 1 |
///
inline real_t maxPeakDifference( const P2Function< real_t >& u, const P2Function< real_t >& uSolution, uint_t level, DoFType flag )
{
   const real_t maxTempApproximate = u.getMaxMagnitude( level, flag );
   const real_t maxTempAnalytical  = uSolution.getMaxMagnitude( level, flag );
   const real_t peakError          = std::abs( maxTempApproximate / maxTempAnalytical - 1 );
   return peakError;
}

/// Calculates and returns var(t) as an indicator for the amount of spurious oscillations:
///
///    var(t) := max(u_h) - min(u_h)
///
inline real_t spuriousOscillations( const P2Function< real_t >& u, uint_t level, DoFType flag )
{
   const real_t maxTemp = u.getMaxValue( level, flag );
   const real_t minTemp = u.getMinValue( level, flag );
   return maxTemp - minTemp;
}

template < typename MassOperator >
inline real_t globalMass( const P2Function< real_t >& u,
                   const P2Function< real_t >& tmp,
                   const MassOperator&         M,
                   const uint_t&               level,
                   const DoFType&              flag )
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

void solve( const MeshInfo&         meshInfo,
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
            bool                    resetParticles,
            bool                    adjustedAdvection,
            uint_t                  numTimeSteps,
            bool                    vtk,
            bool                    vtkOutputVelocity,
            const std::string&      benchmarkName,
            uint_t                  printInterval,
            uint_t                  vtkInterval );


} // namespace moc_benchmarks
} // namespace hyteg