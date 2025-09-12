/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"

#include "hyteg/operators/NoOperator.hpp"

namespace MantleConvection {

class OperatorUpdater
{
 public:
   virtual ~OperatorUpdater() = default;

   OperatorUpdater() {}

   virtual void updateAfterCheckPointLoad()
   {
      WALBERLA_ABORT( "updateAfterCheckPointLoad not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateBeforeSaddlePointSolve()
   {
      WALBERLA_ABORT( "updateBeforeSaddlePointSolve not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterSaddlePointSolve()
   {
      WALBERLA_ABORT( "updateAfterSaddlePointSolve not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterAdvectionDiffusionSolve()
   {
      WALBERLA_ABORT( "updateAfterAdvectionDiffusionSolve not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterViscosityRecalculation()
   {
      WALBERLA_ABORT( "updateAfterViscosityRecalculation not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterDensityRecalculation()
   {
      WALBERLA_ABORT( "updateAfterDensityRecalculation not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterViscosityExtrapolationRecalculation()
   {
      WALBERLA_ABORT( "updateAfterViscosityExtrapolationRecalculation not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterDensityExtrapolationRecalculation()
   {
      WALBERLA_ABORT( "updateAfterDensityExtrapolationRecalculation not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterUpExtrapolationRecalculation()
   {
      WALBERLA_ABORT( "updateAfterUpExtrapolationRecalculation not implemented for OperatorUpdater abstract base class!" );
   }

   virtual void updateAfterTemperatureExtrapolationRecalculation()
   {
      WALBERLA_ABORT(
          "updateAfterTemperatureExtrapolationRecalculation not implemented for OperatorUpdater abstract base class!" );
   }
};

class IdentityOperatorUpdater : public OperatorUpdater
{
 public:
   IdentityOperatorUpdater() {}

   void updateAfterSaddlePointSolve() override {}

   void updateAfterAdvectionDiffusionSolve() override {}

   void updateAfterViscosityRecalculation() override {}

   void updateAfterDensityRecalculation() override {}

   void updateAfterViscosityExtrapolationRecalculation() override {}

   void updateAfterDensityExtrapolationRecalculation() override {}

   void updateAfterUpExtrapolationRecalculation() override {}

   void updateAfterTemperatureExtrapolationRecalculation() override {}
};

template < class SaddlePointOperatorType_         = hyteg::NoOperator,
           class InvKMassOperatorType_            = hyteg::NoOperator,
           class ABlockChebyshevSmootherType_     = hyteg::NoOperator,
           class ABlockPETScCoarseGridSolverType_ = hyteg::NoOperator >
class MantleConvectionOperatorUpdater : public OperatorUpdater
{
 public:
   MantleConvectionOperatorUpdater( const std::shared_ptr< SaddlePointOperatorType_ >&         saddlePointOp,
                                    const std::shared_ptr< InvKMassOperatorType_ >&            invKMass,
                                    const std::shared_ptr< ABlockChebyshevSmootherType_ >&     ABlockChebyshevSmoother,
                                    const std::shared_ptr< ABlockPETScCoarseGridSolverType_ >& ABlockPETScCoarseGridSolver )
   : saddlePointOp_( saddlePointOp )
   , invKMass_( invKMass )
   , ABlockChebyshevSmoother_( ABlockChebyshevSmoother )
   , ABlockPETScCoarseGridSolver_( ABlockPETScCoarseGridSolver )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if constexpr ( !std::is_same< SaddlePointOperatorType_, hyteg::NoOperator >::value )
      {
         if ( saddlePointOp == nullptr )
         {
            WALBERLA_ABORT( "saddlePointOp set to nullptr but type is not NoOperator!" );
         }
      }

      if constexpr ( !std::is_same< InvKMassOperatorType_, hyteg::NoOperator >::value )
      {
         if ( invKMass == nullptr )
         {
            WALBERLA_ABORT( "invKMass set to nullptr but type is not NoOperator!" );
         }
      }

      if constexpr ( !std::is_same< ABlockChebyshevSmootherType_, hyteg::NoOperator >::value )
      {
         if ( ABlockChebyshevSmoother == nullptr )
         {
            WALBERLA_ABORT( "ABlockChebyshevSmoother set to nullptr but type is not NoOperator!" );
         }
      }

      if constexpr ( ( !std::is_same< ABlockPETScCoarseGridSolverType_, hyteg::NoOperator >::value ) &&
                     ( !std::is_same< SaddlePointOperatorType_, hyteg::NoOperator >::value ) )
      {
         if ( ABlockPETScCoarseGridSolver == nullptr )
         {
            WALBERLA_ABORT( "ABlockPETScCoarseGridSolver set to nullptr but type is not NoOperator!" );
         }
      }
   }

   void updateAfterSaddlePointSolve() override {}

   void updateAfterAdvectionDiffusionSolve() override {}

   void updateAfterViscosityRecalculation() override
   {
      if constexpr ( !std::is_same< SaddlePointOperatorType_, hyteg::NoOperator >::value )
      {
         saddlePointOp_->computeInverseDiagonalOperatorValues();
      }

      if constexpr ( !std::is_same< InvKMassOperatorType_, hyteg::NoOperator >::value )
      {
         invKMass_->computeInverseDiagonalOperatorValues();
      }

      if constexpr ( !std::is_same< ABlockChebyshevSmootherType_, hyteg::NoOperator >::value )
      {
         ABlockChebyshevSmoother_->regenerate();
      }

      if constexpr ( ( !std::is_same< ABlockPETScCoarseGridSolverType_, hyteg::NoOperator >::value ) &&
                     ( !std::is_same< SaddlePointOperatorType_, hyteg::NoOperator >::value ) )
      {
         ABlockPETScCoarseGridSolver_->reassembleMatrix( saddlePointOp_->getA() );
      }
   }

   void updateAfterDensityRecalculation() override
   {
      // TODO: update density derivative approximation if needed
   }

   void updateAfterViscosityExtrapolationRecalculation() override {}

   void updateAfterDensityExtrapolationRecalculation() override {}

   void updateAfterUpExtrapolationRecalculation() override {}

   void updateAfterTemperatureExtrapolationRecalculation() override {}

 private:
   std::shared_ptr< SaddlePointOperatorType_ >         saddlePointOp_;
   std::shared_ptr< InvKMassOperatorType_ >            invKMass_;
   std::shared_ptr< ABlockChebyshevSmootherType_ >     ABlockChebyshevSmoother_;
   std::shared_ptr< ABlockPETScCoarseGridSolverType_ > ABlockPETScCoarseGridSolver_;
};

} // namespace MantleConvection