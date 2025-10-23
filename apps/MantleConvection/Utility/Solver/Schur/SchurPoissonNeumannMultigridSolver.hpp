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
class SchurPoissonNeumannMultigridSolver : public SchurSolver< SchurOperator< typename OperatorType::srcType > >
{
 public:
   using SchurSolver< SchurOperator< typename OperatorType::srcType > >::prefix_;
   typedef typename OperatorType::srcType PressureFunctionType;

   SchurPoissonNeumannMultigridSolver(
       walberla::Config::BlockHandle&                                                 parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       uint_t                                                                         minLevel,
       uint_t                                                                         maxLevel,
       const std::shared_ptr< OperatorType >&                                         op,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& SchurSmoother,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& SchurCoarseGridSolver,
       const std::shared_ptr< RestrictionOperatorType >&                              SchurRestriction,
       const std::shared_ptr< ProlongationOperatorType >&                             SchurProlongation,
       bool                                                                           lowMemoryMode = false,
       std::string                                                                    prefix        = "" )
   : SchurSolver< SchurOperator< PressureFunctionType > >( prefix )
   , op_( op )
   , SchurSmoother_( SchurSmoother )
   , SchurCoarseGridSolver_( SchurCoarseGridSolver )
   , SchurRestriction_( SchurRestriction )
   , SchurProlongation_( SchurProlongation )
   , lowMemoryMode_( lowMemoryMode )
   {
      preSmoothingSteps_ =
          parameters.getParameter< uint_t >( prefix + std::string( "SchurPoissonNeumannMultigridPreSmoothingSteps" ) );
      postSmoothingSteps_ =
          parameters.getParameter< uint_t >( prefix + std::string( "SchurPoissonNeumannMultigridPostSmoothingSteps" ) );
      smoothingStepsIncrement_ =
          parameters.getParameter< uint_t >( prefix + std::string( "SchurPoissonNeumannMultigridSmoothingStepsIncrement" ) );
      numberOfCycles_ = parameters.getParameter< uint_t >( prefix + std::string( "SchurPoissonNeumannMultigridNumberOfCycles" ) );

      std::string cycleTypeName =
          parameters.getParameter< std::string >( prefix + std::string( "SchurPoissonNeumannMultigridCycleType" ) );
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

      substSmoother_ = std::make_shared< hyteg::SubstitutePreconditioner< OperatorType, SchurOperator< PressureFunctionType > > >(
          SchurSmoother_, substOperator_ );

      substCoarseGridSolver_ =
          std::make_shared< hyteg::SubstitutePreconditioner< OperatorType, SchurOperator< PressureFunctionType > > >(
              SchurCoarseGridSolver_, substOperator_ );

      SchurPoissonNeumannMultigridSolver_ =
          std::make_shared< hyteg::GeometricMultigridSolver< OperatorType > >( storage,
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
      WALBERLA_UNUSED( A );
      for ( uint_t i = 0; i < numberOfCycles_; i++ )
      {
         SchurPoissonNeumannMultigridSolver_->solve( *op_, x, b, level );
      }
   };

   std::shared_ptr< hyteg::GeometricMultigridSolver< OperatorType > > getSolver() { return SchurPoissonNeumannMultigridSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "############################################################"                                            << "\n";
      os << std::string( offset, ' ') << "########## Schur Poisson Neumann Multigrid Solver ##########"                                            << "\n";
      os << std::string( offset, ' ') << "############################################################"                                            << "\n";
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
   std::shared_ptr< OperatorType >                          op_;
   std::shared_ptr< SchurOperator< PressureFunctionType > > substOperator_;

   std::shared_ptr< hyteg::SubstitutePreconditioner< OperatorType, SchurOperator< PressureFunctionType > > > substSmoother_;
   std::shared_ptr< hyteg::SubstitutePreconditioner< OperatorType, SchurOperator< PressureFunctionType > > >
                                                                      substCoarseGridSolver_;
   std::shared_ptr< hyteg::GeometricMultigridSolver< OperatorType > > SchurPoissonNeumannMultigridSolver_;

   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurSmoother_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurCoarseGridSolver_;
   std::shared_ptr< RestrictionOperatorType >                              SchurRestriction_;
   std::shared_ptr< ProlongationOperatorType >                             SchurProlongation_;

   uint_t preSmoothingSteps_;
   uint_t postSmoothingSteps_;
   uint_t smoothingStepsIncrement_;
   uint_t numberOfCycles_;

   hyteg::CycleType multigridCycleType_;
   std::string      cycleTypeStr_;

   bool lowMemoryMode_;
};

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
inline std::ostream& operator<<(
    std::ostream&                                                                                                os,
    const SchurPoissonNeumannMultigridSolver< OperatorType, RestrictionOperatorType, ProlongationOperatorType >& spnmgs )
{
   return spnmgs.print( os );
}

} // namespace MantleConvection