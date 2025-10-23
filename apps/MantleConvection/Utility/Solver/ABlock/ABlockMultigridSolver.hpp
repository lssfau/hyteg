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

#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
class ABlockMultigridSolver : public ABlockSolver< OperatorType >
{
 public:
   using ABlockSolver< OperatorType >::prefix_;

   ABlockMultigridSolver( walberla::Config::BlockHandle&                                           parameters,
                          const std::shared_ptr< hyteg::PrimitiveStorage >&                        storage,
                          uint_t                                                                   minLevel,
                          uint_t                                                                   maxLevel,
                          const std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > >& ABlockSmoother,
                          const std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > >& ABlockCoarseGridSolver,
                          const std::shared_ptr< RestrictionOperatorType >&                        ABlockRestriction,
                          const std::shared_ptr< ProlongationOperatorType >&                       ABlockProlongation,
                          bool                                                                     lowMemoryMode = false,
                          std::string                                                              prefix        = "" )
   : ABlockSolver< OperatorType >( prefix )
   , ABlockSmoother_( ABlockSmoother )
   , ABlockCoarseGridSolver_( ABlockCoarseGridSolver )
   , ABlockRestriction_( ABlockRestriction )
   , ABlockProlongation_( ABlockProlongation )
   , lowMemoryMode_( lowMemoryMode )
   {
      preSmoothingSteps_  = parameters.getParameter< uint_t >( prefix + std::string( "ABlockMultigridPreSmoothingSteps" ) );
      postSmoothingSteps_ = parameters.getParameter< uint_t >( prefix + std::string( "ABlockMultigridPostSmoothingSteps" ) );
      smoothingStepsIncrement_ =
          parameters.getParameter< uint_t >( prefix + std::string( "ABlockMultigridSmoothingStepsIncrement" ) );
      numberOfCycles_ = parameters.getParameter< uint_t >( prefix + std::string( "ABlockMultigridNumberOfCycles" ) );

      std::string cycleTypeName = parameters.getParameter< std::string >( prefix + std::string( "ABlockMultigridCycleType" ) );
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

      ABlockMultigridSolver_ = std::make_shared< hyteg::GeometricMultigridSolver< OperatorType > >( storage,
                                                                                                    ABlockSmoother_,
                                                                                                    ABlockCoarseGridSolver_,
                                                                                                    ABlockRestriction_,
                                                                                                    ABlockProlongation_,
                                                                                                    minLevel,
                                                                                                    maxLevel,
                                                                                                    preSmoothingSteps_,
                                                                                                    postSmoothingSteps_,
                                                                                                    smoothingStepsIncrement_,
                                                                                                    multigridCycleType_,
                                                                                                    lowMemoryMode_ );
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      for ( uint_t i = 0; i < numberOfCycles_; i++ )
      {
         ABlockMultigridSolver_->solve( A, x, b, level );
      }
   };

   std::shared_ptr< hyteg::GeometricMultigridSolver< OperatorType > > getSolver() { return ABlockMultigridSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "############ ABlock Multigrid Solver ############"                                                       << "\n";
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "preSmoothingSteps_: "         << preSmoothingSteps_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "postSmoothingSteps_: "        << postSmoothingSteps_         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "smoothingStepsIncrement_: "   << smoothingStepsIncrement_    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "multigridCycleType_: "        << cycleTypeStr_               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "numberOfCycles_: "            << numberOfCycles_             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";
      // clang-format on

      ABlockSmoother_->print( os, offset + 3 );
      os << "\n";
      ABlockCoarseGridSolver_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > >  ABlockSmoother_;
   std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > >  ABlockCoarseGridSolver_;
   std::shared_ptr< RestrictionOperatorType >                         ABlockRestriction_;
   std::shared_ptr< ProlongationOperatorType >                        ABlockProlongation_;
   std::shared_ptr< hyteg::GeometricMultigridSolver< OperatorType > > ABlockMultigridSolver_;

   uint_t           preSmoothingSteps_;
   uint_t           postSmoothingSteps_;
   uint_t           smoothingStepsIncrement_;
   uint_t           numberOfCycles_;
   hyteg::CycleType multigridCycleType_;
   std::string      cycleTypeStr_;

   bool lowMemoryMode_;
};

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
inline std::ostream&
    operator<<( std::ostream&                                                                                   os,
                const ABlockMultigridSolver< OperatorType, RestrictionOperatorType, ProlongationOperatorType >& abmgs )
{
   return abmgs.print( os );
}

} // namespace MantleConvection