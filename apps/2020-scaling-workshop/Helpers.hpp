/*
 * Copyright (c) 2017-2020 Nils Kohl, Dominik Thoennes.
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
#include "core/OpenMP.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/Git.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/OpenMPManager.hpp"
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
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
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
namespace scaling_workshop {

inline std::string getDateTimeID()
{
   std::vector< char > cTimeString( 64 );
   WALBERLA_ROOT_SECTION()
   {
      std::time_t t;
      std::time( &t );
      std::strftime( cTimeString.data(), 64, "%F_%H-%M-%S", std::localtime( &t ) );
   }

   walberla::mpi::broadcastObject( cTimeString );

   std::string timeString( cTimeString.data() );
   return timeString;
}

template < template < typename > class FunctionType >
inline real_t pointwiseScaledL2Norm( const FunctionType< real_t >& u, uint_t level )
{
   auto numUnknowns = numberOfGlobalDoFs< typename FunctionType< real_t >::Tag >( *u.getStorage(), level );
   auto scalarProd  = u.dotGlobal( u, level, All );
   return std::sqrt( scalarProd / real_c( numUnknowns ) );
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
   r.assign( {1.0, -1.0}, {f, tmp}, level, flag );
}

template < template < typename > class StokesFunctionType, typename StokesOperator >
inline void residualNegativeRHS0( const StokesFunctionType< real_t >& u,
                                  const StokesOperator&               A,
                                  uint_t                              level,
                                  DoFType                             flag,
                                  StokesFunctionType< real_t >&       r )
{
   A.apply( u, r, level, flag );
}

template < template < typename > class StokesFunctionType >
inline void error( const StokesFunctionType< real_t >& u,
                   const StokesFunctionType< real_t >& exact,
                   uint_t                              level,
                   DoFType                             flag,
                   StokesFunctionType< real_t >&       e )
{
   e.assign( {1.0, -1.0}, {u, exact}, level, flag );
}

template < template < typename > class StokesFunctionType, typename StokesOperator >
inline void errorAndResidual( const StokesOperator&               A,
                              const StokesFunctionType< real_t >& u,
                              const StokesFunctionType< real_t >& f,
                              StokesFunctionType< real_t >& r,
                              StokesFunctionType< real_t >& tmp,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionU,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionV,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionW,
                              const std::function< real_t( const hyteg::Point3D& ) >& solutionP,
                              uint_t                                                  level,
                              DoFType                                                 flag,
                              bool                                                    RHSisZero,
                              real_t&                                                 residualL2Velocity,
                              real_t&                                                 residualL2Pressure,
                              real_t&                                                 errorL2Velocity,
                              real_t&                                                 errorL2Pressure )
{
   if ( RHSisZero )
   {
      tmp.interpolate( 0, level );
      residualNegativeRHS0( u, A, level, flag, tmp );
      residualL2Velocity = pointwiseScaledL2Norm( tmp.uvw(), level );
      residualL2Pressure = pointwiseScaledL2Norm( tmp.p(), level );

      auto errorU = [&]( const Point3D& p, const std::vector< real_t >& values )
      {
         auto sol = solutionU(p);
         return sol - values[0];
      };

      auto errorV = [&]( const Point3D& p, const std::vector< real_t >& values )
      {
        auto sol = solutionV(p);
        return sol - values[0];
      };

      auto errorW = [&]( const Point3D& p, const std::vector< real_t >& values )
      {
        auto sol = solutionW(p);
        return sol - values[0];
      };

      auto errorP = [&]( const Point3D& p, const std::vector< real_t >& values )
      {
        auto sol = solutionP(p);
        return sol - values[0];
      };

      tmp.uvw()[0].interpolate( errorU, {u.uvw()[0]}, level, flag );
      tmp.uvw()[1].interpolate( errorV, {u.uvw()[1]}, level, flag );
      tmp.uvw()[2].interpolate( errorW, {u.uvw()[2]}, level, flag );
      tmp.p().interpolate( errorP, {u.p()}, level, flag );

      errorL2Velocity = pointwiseScaledL2Norm( tmp.uvw(), level );
      errorL2Pressure = pointwiseScaledL2Norm( tmp.p(), level );
      tmp.interpolate( 0, level );
   }
   else
   {
      tmp.interpolate( 0, level );
      residual( u, f, A, tmp, level, flag, r );
      residualL2Velocity = pointwiseScaledL2Norm( r.uvw(), level );
      residualL2Pressure = pointwiseScaledL2Norm( r.p(), level );

      tmp.uvw()[0].interpolate( solutionU, level, All );
      tmp.uvw()[1].interpolate( solutionV, level, All );
      tmp.uvw()[2].interpolate( solutionW, level, All );
      tmp.p().interpolate( solutionP, level, All );
      error( u, tmp, level, flag, r );

      errorL2Velocity = pointwiseScaledL2Norm( r.uvw(), level );
      errorL2Pressure = pointwiseScaledL2Norm( r.p(), level );
      tmp.interpolate( 0, level );
   }


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
            const std::function< real_t( const hyteg::Point3D& ) >& solutionU,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionV,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionW,
            const std::function< real_t( const hyteg::Point3D& ) >& solutionP,
            const std::function< real_t( const hyteg::Point3D& ) >& initialU,
            const std::function< real_t( const hyteg::Point3D& ) >& initialV,
            const std::function< real_t( const hyteg::Point3D& ) >& initialW,
            const std::function< real_t( const hyteg::Point3D& ) >& initialP,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsU,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsV,
            const std::function< real_t( const hyteg::Point3D& ) >& rhsW,
            uint_t                                                  minLevel,
            uint_t                                                  maxLevel,
            MultigridSettings                                       multigridSettings,
            SmootherSettings                                        smootherSettings,
            CoarseGridSettings                                      coarseGridSettings,
            bool                                                    projectPressure,
            bool                                                    projectPressurefterRestriction,
            bool                                                    vtk,
            std::string                                             benchmarkName,
            FixedSizeSQLDB                                          dbFile,
            std::string                                             timingFile,
            bool                                                    RHSisZero );

void domainInfo( const std::shared_ptr< PrimitiveStorage >& storage,
                 Discretization                             discretization,
                 uint_t                                     minLevel,
                 uint_t                                     maxLevel,
                 bool                                       vtk,
                 std::string                                benchmarkName );

} // namespace scaling_workshop
} // namespace hyteg
