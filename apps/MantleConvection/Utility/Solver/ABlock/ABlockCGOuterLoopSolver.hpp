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
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "ABlockIdentityPreconditioner.hpp"
#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class ABlockCGOuterLoopSolver : public ABlockSolver< OperatorType >
{
 public:
   using ABlockSolver< OperatorType >::prefix_;

   ABlockCGOuterLoopSolver( walberla::Config::BlockHandle&                         parameters,
                            const std::shared_ptr< hyteg::PrimitiveStorage >&      storage,
                            uint_t                                                 minLevel,
                            uint_t                                                 maxLevel,
                            bool                                                   lowMemoryMode = false,
                            const std::shared_ptr< ABlockSolver< OperatorType > >& ABlockPreconditioner =
                                std::make_shared< ABlockIdentityPreconditioner< OperatorType > >(),
                            std::string prefix = "" )
   : ABlockSolver< OperatorType >( prefix )
   , ABlockPreconditioner_( ABlockPreconditioner )
   , lowMemoryMode_( lowMemoryMode )
   {
      maxIterations_     = parameters.getParameter< uint_t >( prefix + std::string( "ABlockCGOuterLoopMaxIterations" ) );
      relativeTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "ABlockCGOuterLoopRelativeTolerance" ) );
      absoluteTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "ABlockCGOuterLoopAbsoluteTolerance" ) );
      printInfo_         = parameters.getParameter< bool >( prefix + std::string( "ABlockCGOuterLoopPrintInfo" ) );

      ABlockCGOuterLoopSolver_ = std::make_shared< hyteg::CGSolver< OperatorType > >( storage,
                                                                                      minLevel,
                                                                                      maxLevel,
                                                                                      maxIterations_,
                                                                                      relativeTolerance_,
                                                                                      absoluteTolerance_,
                                                                                      ABlockPreconditioner_,
                                                                                      lowMemoryMode_ );
      ABlockCGOuterLoopSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         ABlockCGOuterLoopSolver_->setName( std::string( "CG ABlock Outer Loop prefix " ) + prefix_ );
      }
      else
      {
         ABlockCGOuterLoopSolver_->setName( "CG ABlock Outer Loop" );
      }
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      ABlockCGOuterLoopSolver_->solve( A, x, b, level );
   };

   std::shared_ptr< hyteg::CGSolver< OperatorType > > getSolver() { return ABlockCGOuterLoopSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "########## ABlock CG Outer Loop Solver ##########"                                                       << "\n";
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxIterations_: "             << maxIterations_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relativeTolerance_: "         << relativeTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "printInfo_: "                 << printInfo_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";
      // clang-format on

      ABlockPreconditioner_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::CGSolver< OperatorType > > ABlockCGOuterLoopSolver_;
   std::shared_ptr< ABlockSolver< OperatorType > >    ABlockPreconditioner_;
   uint_t                                             maxIterations_;
   real_t                                             relativeTolerance_;
   real_t                                             absoluteTolerance_;
   bool                                               printInfo_;

   bool lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const ABlockCGOuterLoopSolver< OperatorType >& abcgols )
{
   return abcgols.print( os );
}

} // namespace MantleConvection