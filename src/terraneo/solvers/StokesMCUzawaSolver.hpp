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
#include "hyteg/solvers/solvertemplates/StokesFSGMGUzawaSolverTemplate.hpp"

#include "terraneo/helpers/TerraNeoParameters.hpp"
#include "terraneo/solvers/MCSolverBase.hpp"

namespace hyteg {
template < typename StokesOperator_T, typename ProjectionOperator_T >
class StokesMCUzawaSolver : public MCSolverBase< StokesOperator_T >
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

   StokesMCUzawaSolver( const std::shared_ptr< PrimitiveStorage >&                 storage,
                        const uint_t                                               minLevel,
                        const uint_t                                               maxLevel,
                        const std::shared_ptr< StokesOperator_T >&                 stokesOperator,
                        const std::shared_ptr< ProjectionOperator_T >&             projectionOperator,
                        const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp1,
                        const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp2,
                        terraneo::TerraNeoParameters&                              TN )
   : MCSolverBase< StokesOperator_T >( storage, minLevel, maxLevel )
   , TN_( TN )
   , stokesOperator_( stokesOperator )
   , projectionOperator_( projectionOperator )
   , temp1_( temp1 )
   , temp2_( temp2 )
   {
      auto stopIterationCallback = [this]( const StokesOperator_T&                 A,
                                           const P2P1TaylorHoodFunction< real_t >& u,
                                           const P2P1TaylorHoodFunction< real_t >& b,
                                           const uint_t                            level ) {
         projectionOperator_->project( b, TN_.domainParameters.maxLevel, FreeslipBoundary );
         A.apply( u, *temp1_, level, Inner | NeumannBoundary | FreeslipBoundary );
         temp1_->assign( { real_c( 1 ), real_c( -1 ) }, { *temp1_, b }, level, Inner | NeumannBoundary | FreeslipBoundary );

         real_t stokesResidual = std::sqrt( temp1_->dotGlobal( *temp1_, level, Inner | NeumannBoundary | FreeslipBoundary ) );

         if ( TN_.solverParameters.numVCycles == 0 )
         {
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
                "[Uzawa] iter %3d | residual: %10.5e | initial ", 0, TN_.solverParameters.vCycleResidualUPrev ) );
         }

         auto reductionRateU = stokesResidual / TN_.solverParameters.vCycleResidualUPrev;

         TN_.solverParameters.vCycleResidualUPrev = stokesResidual;

         TN_.solverParameters.numVCycles++;
         TN_.solverParameters.averageResidualReductionU += reductionRateU;

         if ( TN_.simulationParameters.verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "[Uzawa] iter %3d | residual: %10.5e | reduction: %10.5e ",
                                                         TN_.solverParameters.numVCycles,
                                                         stokesResidual,
                                                         reductionRateU ) );
         }

         if ( stokesResidual / TN_.solverParameters.initialResidualU < TN_.solverParameters.stokesRelativeResidualUTolerance )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[Uzawa] reached relative residual threshold" )
            return true;
         }

         if ( stokesResidual < TN_.solverParameters.stokesAbsoluteResidualUTolerance )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[Uzawa] reached absolute residual threshold" )
            return true;
         }
         return false;
      };

      auto solverContainer = hyteg::solvertemplates::stokesGMGUzawaFSSolver(
          storage, minLevel, maxLevel, stokesOperator, projectionOperator, temp1, temp2, true, {} );

      this->solverPtr_ = std::make_shared< SolverLoop< StokesOperator_T > >(
          std::get< 0 >( solverContainer ), TN_.solverParameters.stokesMaxNumIterations, stopIterationCallback );
   }

 private:
   terraneo::TerraNeoParameters& TN_;

   const std::shared_ptr< StokesOperator_T >&     stokesOperator_;
   const std::shared_ptr< ProjectionOperator_T >& projectionOperator_;

   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp1_;
   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp2_;
};
} // namespace hyteg
