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
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"

#include "../../OperatorTools/ApproximateSchurComplementOperator.hpp"
#include "../../OperatorTools/SchurOperator.hpp"
#include "../ABlock/ABlockSolver.hpp"
#include "SchurSolver.hpp"

namespace MantleConvection {

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
class SchurMultigridSolver : public SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >
{
 public:
   using SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >::prefix_;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef ApproximateSchurComplementOperator< OperatorType > ApproximateSchurOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   SchurMultigridSolver(
       walberla::Config::BlockHandle&                                                                               parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                                                            storage,
       uint_t                                                                                                       minLevel,
       uint_t                                                                                                       maxLevel,
       const std::shared_ptr< OperatorType >&                                                                       saddleOp,
       const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >&                                    ABlockSolver,
       hyteg::BoundaryCondition                                                                                     velocityBC,
       const std::shared_ptr< SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > > >& SchurSmoother,
       const std::shared_ptr< SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > > >&
                                                             SchurCoarseGridSolver,
       const std::shared_ptr< RestrictionOperatorType >&     SchurRestriction,
       const std::shared_ptr< ProlongationOperatorType >&    SchurProlongation,
       const std::shared_ptr< ApproximateSchurOperatorType > approximateSchurOperator = nullptr,
       bool                                                  lowMemoryMode            = false,
       std::string                                           prefix                   = "" )
   : SchurMultigridSolver( parameters,
                           storage,
                           minLevel,
                           maxLevel,
                           saddleOp,
                           ABlockSolver,
                           velocityBC,
                           velocityBC,
                           velocityBC,
                           SchurSmoother,
                           SchurCoarseGridSolver,
                           SchurRestriction,
                           SchurProlongation,
                           approximateSchurOperator,
                           lowMemoryMode,
                           prefix )
   {}

   SchurMultigridSolver(
       walberla::Config::BlockHandle&                                                                               parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                                                            storage,
       uint_t                                                                                                       minLevel,
       uint_t                                                                                                       maxLevel,
       const std::shared_ptr< OperatorType >&                                                                       saddleOp,
       const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >&                                    ABlockSolver,
       hyteg::BoundaryCondition                                                                                     velocityBCx,
       hyteg::BoundaryCondition                                                                                     velocityBCy,
       hyteg::BoundaryCondition                                                                                     velocityBCz,
       const std::shared_ptr< SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > > >& SchurSmoother,
       const std::shared_ptr< SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > > >&
                                                             SchurCoarseGridSolver,
       const std::shared_ptr< RestrictionOperatorType >&     SchurRestriction,
       const std::shared_ptr< ProlongationOperatorType >&    SchurProlongation,
       const std::shared_ptr< ApproximateSchurOperatorType > approximateSchurOperator = nullptr,
       bool                                                  lowMemoryMode            = false,
       std::string                                           prefix                   = "" )
   : SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >( prefix )
   , SchurSmoother_( SchurSmoother )
   , SchurCoarseGridSolver_( SchurCoarseGridSolver )
   , SchurRestriction_( SchurRestriction )
   , SchurProlongation_( SchurProlongation )
   , lowMemoryMode_( lowMemoryMode )
   {
      if ( approximateSchurOperator != nullptr )
      {
         approximateSchurOperator_ = approximateSchurOperator;
      }
      else
      {
         approximateSchurOperator_ = std::make_shared< ApproximateSchurOperatorType >(
             storage, minLevel, maxLevel, saddleOp, ABlockSolver, velocityBCx, velocityBCy, velocityBCz, lowMemoryMode_ );
      }

      preSmoothingSteps_  = parameters.getParameter< uint_t >( prefix + std::string( "SchurMultigridPreSmoothingSteps" ) );
      postSmoothingSteps_ = parameters.getParameter< uint_t >( prefix + std::string( "SchurMultigridPostSmoothingSteps" ) );
      smoothingStepsIncrement_ =
          parameters.getParameter< uint_t >( prefix + std::string( "SchurMultigridSmoothingStepsIncrement" ) );
      numberOfCycles_ = parameters.getParameter< uint_t >( prefix + std::string( "SchurMultigridNumberOfCycles" ) );

      std::string cycleTypeName = parameters.getParameter< std::string >( prefix + std::string( "SchurMultigridCycleType" ) );
      if ( cycleTypeName == "WCycle" )
      {
         multigridCycleType_ = hyteg::CycleType::WCYCLE;
         cycleTypeStr_       = "WCycle";
      }
      else
      {
         multigridCycleType_ = hyteg::CycleType::VCYCLE;
         cycleTypeStr_       = "VCycle";
      }

      substOperator_ = std::make_shared< SchurOperator< PressureFunctionType > >( storage, minLevel, maxLevel );

      substSmoother_ = std::make_shared<
          hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >(
          SchurSmoother_, substOperator_ );

      substCoarseGridSolver_ = std::make_shared<
          hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >(
          SchurCoarseGridSolver_, substOperator_ );

      SchurMultigridSolver_ =
          std::make_shared< hyteg::GeometricMultigridSolver< ApproximateSchurOperatorType > >( storage,
                                                                                               substSmoother_,
                                                                                               substCoarseGridSolver_,
                                                                                               SchurRestriction_,
                                                                                               SchurProlongation_,
                                                                                               minLevel,
                                                                                               maxLevel,
                                                                                               preSmoothingSteps_,
                                                                                               postSmoothingSteps_,
                                                                                               smoothingStepsIncrement_,
                                                                                               multigridCycleType_,
                                                                                               lowMemoryMode_ );
   }

   void solve( const SchurOperator< PressureFunctionType >& A,
               const PressureFunctionType&                  x,
               const PressureFunctionType&                  b,
               const walberla::uint_t                       level ) override
   {
      for ( uint_t i = 0; i < numberOfCycles_; i++ )
      {
         SchurMultigridSolver_->solve( *approximateSchurOperator_, x, b, level );
      }
   };

   std::shared_ptr< hyteg::GeometricMultigridSolver< ApproximateSchurOperatorType > > getSolver()
   {
      return SchurMultigridSolver_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "##################################################"                                                      << "\n";
      os << std::string( offset, ' ') << "############# Schur Multigrid Solver #############"                                                      << "\n";
      os << std::string( offset, ' ') << "##################################################"                                                      << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "preSmoothingSteps_: "         << preSmoothingSteps_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "postSmoothingSteps_: "        << postSmoothingSteps_         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "smoothingStepsIncrement_: "   << smoothingStepsIncrement_    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "multigridCycleType_: "        << cycleTypeStr_               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "numberOfCycles_: "            << numberOfCycles_             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";
      // clang-format on

      SchurSmoother_->print( os, offset + 3 );
      os << "\n";
      SchurCoarseGridSolver_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< OperatorType >                          saddleOp_;
   std::shared_ptr< ApproximateSchurOperatorType >          approximateSchurOperator_;
   std::shared_ptr< SchurOperator< PressureFunctionType > > substOperator_;
   std::shared_ptr< hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >
       substSmoother_;
   std::shared_ptr< hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >
       substCoarseGridSolver_;

   std::shared_ptr< SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > > > SchurSmoother_;
   std::shared_ptr< SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > > > SchurCoarseGridSolver_;
   std::shared_ptr< RestrictionOperatorType >                                                            SchurRestriction_;
   std::shared_ptr< ProlongationOperatorType >                                                           SchurProlongation_;
   std::shared_ptr< hyteg::GeometricMultigridSolver< ApproximateSchurOperatorType > >                    SchurMultigridSolver_;

   uint_t preSmoothingSteps_;
   uint_t postSmoothingSteps_;
   uint_t smoothingStepsIncrement_;
   uint_t numberOfCycles_;

   hyteg::CycleType multigridCycleType_;
   std::string      cycleTypeStr_;

   bool lowMemoryMode_;
};

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
inline std::ostream&
    operator<<( std::ostream&                                                                                  os,
                const SchurMultigridSolver< OperatorType, RestrictionOperatorType, ProlongationOperatorType >& sms )
{
   return sms.print( os );
}

} // namespace MantleConvection