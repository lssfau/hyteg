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

#include "SaddlePointSolver.hpp"

namespace MantleConvection {

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
class SaddlePointMultigridSolver : public SaddlePointSolver< OperatorType >
{
 public:
   using SaddlePointSolver< OperatorType >::prefix_;

   SaddlePointMultigridSolver(
       walberla::Config::BlockHandle&                                                parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                             storage,
       uint_t                                                                        minLevel,
       uint_t                                                                        maxLevel,
       const std::shared_ptr< MantleConvection::SaddlePointSolver< OperatorType > >& SaddlePointSmoother,
       const std::shared_ptr< MantleConvection::SaddlePointSolver< OperatorType > >& SaddlePointCoarseGridSolver,
       const std::shared_ptr< RestrictionOperatorType >&                             SaddlePointRestriction,
       const std::shared_ptr< ProlongationOperatorType >&                            SaddlePointProlongation,
       bool                                                                          lowMemoryMode = false,
       std::string                                                                   prefix        = "" )
   : SaddlePointSolver< OperatorType >( prefix )
   , SaddlePointSmoother_( SaddlePointSmoother )
   , SaddlePointCoarseGridSolver_( SaddlePointCoarseGridSolver )
   , SaddlePointRestriction_( SaddlePointRestriction )
   , SaddlePointProlongation_( SaddlePointProlongation )
   , lowMemoryMode_( lowMemoryMode )
   {
      preSmoothingSteps_  = parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointMultigridPreSmoothingSteps" ) );
      postSmoothingSteps_ = parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointMultigridPostSmoothingSteps" ) );
      smoothingStepsIncrement_ =
          parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointMultigridSmoothingStepsIncrement" ) );
      numberOfCycles_ = parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointMultigridNumberOfCycles" ) );

      std::string cycleTypeName =
          parameters.getParameter< std::string >( prefix + std::string( "SaddlePointMultigridCycleType" ) );
      if ( cycleTypeName == "WCycle" )
      {
         multigridCycleType_ = CycleType::WCYCLE;
         cycleTypeStr_       = "WCycle";
      }
      else
      {
         multigridCycleType_ = CycleType::VCYCLE;
         cycleTypeStr_       = "VCycle";
      }

      SaddlePointMultigridSolver_ = std::make_shared< GeometricMultigridSolver< OperatorType > >( storage,
                                                                                                  SaddlePointSmoother_,
                                                                                                  SaddlePointCoarseGridSolver_,
                                                                                                  SaddlePointRestriction_,
                                                                                                  SaddlePointProlongation_,
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
         SaddlePointMultigridSolver_->solve( A, x, b, level );
      }
   };

   std::shared_ptr< hyteg::GeometricMultigridSolver< OperatorType > > getSolver() { return SaddlePointMultigridSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "######### Saddle Point Multigrid Solver #########"                                                       << "\n";
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

      SaddlePointSmoother_->print( os, offset + 3 );
      os << "\n";
      SaddlePointCoarseGridSolver_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< MantleConvection::SaddlePointSolver< OperatorType > > SaddlePointSmoother_;
   std::shared_ptr< MantleConvection::SaddlePointSolver< OperatorType > > SaddlePointCoarseGridSolver_;
   std::shared_ptr< RestrictionOperatorType >                             SaddlePointRestriction_;
   std::shared_ptr< ProlongationOperatorType >                            SaddlePointProlongation_;
   std::shared_ptr< hyteg::GeometricMultigridSolver< OperatorType > >     SaddlePointMultigridSolver_;

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
    operator<<( std::ostream&                                                                                        os,
                const SaddlePointMultigridSolver< OperatorType, RestrictionOperatorType, ProlongationOperatorType >& spmgs )
{
   return spmgs.print( os );
}

} // namespace MantleConvection