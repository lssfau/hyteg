/*
 * Copyright (c) 2025 Ponsuganth Ilangovan P
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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/solvertemplates/StokesFSGMGSolverTemplate.hpp"

#include "terraneo/helpers/TerraNeoParameters.hpp"
#include "terraneo/solvers/MCSolverBase.hpp"

namespace hyteg {
template < typename StokesOperator_T, typename ProjectionOperator_T >
class StokesMCFGMRESSolverWithProjection : public MCSolverBase< StokesOperator_T >
{
 public:
   // clang-format off
   enum class Parameters
   {
    POWER_ITERATIONS,
        
        FGMRES_OUTER_ITER,
        FGMRES_OUTER_TOL,
            ABLOCK_MG_PRESMOOTH,    
            ABLOCK_MG_POSTSMOOTH,
                ABLOCK_COARSE_ITER,
                ABLOCK_COARSE_TOL,
                ABLOCK_COARSE_GRID_PETSC,
            
            SCHUR_CG_SOLVER_ITER,
            SCHUR_CG_SOLVER_TOL,
   };
   // clang-format on

   StokesMCFGMRESSolverWithProjection( const std::shared_ptr< PrimitiveStorage >&                 storage,
                                       const uint_t                                               minLevel,
                                       const uint_t                                               maxLevel,
                                       const std::shared_ptr< StokesOperator_T >&                 stokesOperator,
                                       const std::shared_ptr< ProjectionOperator_T >&             projectionOperator,
                                       const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp1,
                                       const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp2,
                                       const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp3,
                                       terraneo::TerraNeoParameters&                              TN )
   : MCSolverBase< StokesOperator_T >( storage, minLevel, maxLevel )
   , TN_( TN )
   , stokesOperator_( stokesOperator )
   , projectionOperator_( projectionOperator )
   , temp1_( temp1 )
   , temp2_( temp2 )
   , temp3_( temp3 )
   {
      auto solverContainer = hyteg::solvertemplates::stokesGMGFSSolver(
          storage,
          minLevel,
          maxLevel,
          stokesOperator,
          projectionOperator,
          temp1,
          temp2,
          temp3,
          TN.solverParameters.estimateUzawaOmega,
          TN.simulationParameters.verbose,
          {
              { solvertemplates::StokesGMGFSSolverParamKey::NUM_POWER_ITERATIONS_SPECTRUM,
                TN.solverParameters.numPowerIterations },
              { solvertemplates::StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_ITER,
                TN.solverParameters.FGMRESOuterIterations },
              { solvertemplates::StokesGMGFSSolverParamKey::FGMRES_UZAWA_PRECONDITIONED_OUTER_TOLERANCE,
                TN.solverParameters.FGMRESTolerance },
              { solvertemplates::StokesGMGFSSolverParamKey::INEXACT_UZAWA_VELOCITY_ITER, TN.solverParameters.uzawaIterations },
              { solvertemplates::StokesGMGFSSolverParamKey::INEXACT_UZAWA_OMEGA, TN.solverParameters.uzawaOmega },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_ITER,
                TN.solverParameters.ABlockMGIterations },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE,
                TN.solverParameters.ABlockMGTolerance },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_MG_PRESMOOTH, TN.solverParameters.ABlockMGPreSmooth },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_MG_POSTSMOOTH, TN.solverParameters.ABlockMGPostSmooth },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_COARSE_ITER, TN.solverParameters.ABlockCoarseGridIterations },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_COARSE_TOLERANCE,
                TN.solverParameters.ABlockCoarseGridTolerance },
              { solvertemplates::StokesGMGFSSolverParamKey::ABLOCK_COARSE_GRID_PETSC, TN.solverParameters.solverPETSc },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_ITER,
                TN.solverParameters.SchurMGIterations },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_CG_SOLVER_MG_PRECONDITIONED_TOLERANCE,
                TN.solverParameters.SchurMGTolerance },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_MG_PRESMOOTH, TN.solverParameters.SchurMGPreSmooth },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_MG_POSTSMOOTH, TN.solverParameters.SchurMGPostSmooth },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_ITER,
                TN.solverParameters.SchurCoarseGridIterations },
              { solvertemplates::StokesGMGFSSolverParamKey::SCHUR_COARSE_GRID_CG_TOLERANCE,
                TN.solverParameters.SchurCoarseGridTolerance },
          } );

      this->solverPtr_ = std::get< 0 >( solverContainer );
   }

 private:
   terraneo::TerraNeoParameters& TN_;

   const std::shared_ptr< StokesOperator_T >&     stokesOperator_;
   const std::shared_ptr< ProjectionOperator_T >& projectionOperator_;

   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp1_;
   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp2_;
   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp3_;
};
} // namespace hyteg
