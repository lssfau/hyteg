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
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "SaddlePointIdentityPreconditioner.hpp"
#include "SaddlePointSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class SaddlePointMinResCoarseGridSolver : public SaddlePointSolver< OperatorType >
{
 public:
   using SaddlePointSolver< OperatorType >::prefix_;

   SaddlePointMinResCoarseGridSolver( walberla::Config::BlockHandle&                              parameters,
                                      const std::shared_ptr< hyteg::PrimitiveStorage >&           storage,
                                      uint_t                                                      minLevel,
                                      uint_t                                                      maxLevel,
                                      bool                                                        lowMemoryMode = false,
                                      const std::shared_ptr< SaddlePointSolver< OperatorType > >& SaddlePointPreconditioner =
                                          std::make_shared< SaddlePointIdentityPreconditioner< OperatorType > >(),
                                      std::string prefix = "" )
   : SaddlePointSolver< OperatorType >( prefix )
   , SaddlePointPreconditioner_( SaddlePointPreconditioner )
   , lowMemoryMode_( lowMemoryMode )
   {
      maxIterations_ = parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointMinResCoarseGridMaxIterations" ) );
      relativeTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "SaddlePointMinResCoarseGridRelativeTolerance" ) );
      absoluteTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "SaddlePointMinResCoarseGridAbsoluteTolerance" ) );
      printInfo_ = parameters.getParameter< bool >( prefix + std::string( "SaddlePointMinResCoarseGridPrintInfo" ) );

      SaddlePointMinResCoarseGridSolver_ = std::make_shared< hyteg::MinResSolver< OperatorType > >( storage,
                                                                                                    minLevel,
                                                                                                    maxLevel,
                                                                                                    maxIterations_,
                                                                                                    relativeTolerance_,
                                                                                                    absoluteTolerance_,
                                                                                                    SaddlePointPreconditioner_,
                                                                                                    lowMemoryMode_ );
      SaddlePointMinResCoarseGridSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         SaddlePointMinResCoarseGridSolver_->setName( std::string( "MinRes Saddle Point Coarse Grid prefix " ) + prefix_ );
      }
      else
      {
         SaddlePointMinResCoarseGridSolver_->setName( "MinRes Saddle Point Coarse Grid" );
      }
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      SaddlePointMinResCoarseGridSolver_->solve( A, x, b, level );
   };

   std::shared_ptr< hyteg::MinResSolver< OperatorType > > getSolver() { return SaddlePointMinResCoarseGridSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "################################################"                                                        << "\n";
      os << std::string( offset, ' ') << "#### Saddle Point MinRes Coarse Grid Solver ####"                                                        << "\n";
      os << std::string( offset, ' ') << "################################################"                                                        << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxIterations_: "             << maxIterations_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relativeTolerance_: "         << relativeTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "printInfo_: "                 << printInfo_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";
      // clang-format on

      SaddlePointPreconditioner_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::MinResSolver< OperatorType > > SaddlePointMinResCoarseGridSolver_;
   std::shared_ptr< SaddlePointSolver< OperatorType > >   SaddlePointPreconditioner_;
   uint_t                                                 maxIterations_;
   real_t                                                 relativeTolerance_;
   real_t                                                 absoluteTolerance_;
   bool                                                   printInfo_;

   bool lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const SaddlePointMinResCoarseGridSolver< OperatorType >& spmrs )
{
   return spmrs.print( os );
}

} // namespace MantleConvection