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
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "SaddlePointIdentityPreconditioner.hpp"
#include "SaddlePointSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class SaddlePointFGMRESSolver : public SaddlePointSolver< OperatorType >
{
 public:
   using SaddlePointSolver< OperatorType >::prefix_;

   SaddlePointFGMRESSolver( walberla::Config::BlockHandle&                              parameters,
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
      maxIterations_     = parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointFGMRESMaxIterations" ) );
      relativeTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SaddlePointFGMRESRelativeTolerance" ) );
      absoluteTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SaddlePointFGMRESAbsoluteTolerance" ) );
      printInfo_         = parameters.getParameter< bool >( prefix + std::string( "SaddlePointFGMRESPrintInfo" ) );

      restartLength_    = parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointFGMRESRestartLength" ) );
      arnoldiTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SaddlePointFGMRESArnoldiTolerance" ) );
      doubleOrthoTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "SaddlePointFGMRESDoubleOrthoTolerance" ) );

      SaddlePointFGMRESSolver_ = std::make_shared< hyteg::FGMRESSolver< OperatorType > >( storage,
                                                                                          minLevel,
                                                                                          maxLevel,
                                                                                          maxIterations_,
                                                                                          relativeTolerance_,
                                                                                          absoluteTolerance_,
                                                                                          SaddlePointPreconditioner_,
                                                                                          restartLength_,
                                                                                          arnoldiTolerance_,
                                                                                          doubleOrthoTolerance_,
                                                                                          lowMemoryMode );

      SaddlePointFGMRESSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         SaddlePointFGMRESSolver_->setName( std::string( "FGMRES Saddle Point prefix " ) + prefix_ );
      }
      else
      {
         SaddlePointFGMRESSolver_->setName( "FGMRES Saddle Point" );
      }
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      SaddlePointFGMRESSolver_->solve( A, x, b, level );
   };

   std::shared_ptr< hyteg::FGMRESSolver< OperatorType > > getSolver() { return SaddlePointFGMRESSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "################################################"                                                        << "\n";
      os << std::string( offset, ' ') << "########## Saddle Point FGMRES Solver ##########"                                                        << "\n";
      os << std::string( offset, ' ') << "################################################"                                                        << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "maxIterations_: "             << maxIterations_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "relativeTolerance_: "         << relativeTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "printInfo_: "                 << printInfo_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "restartLength_: "             << restartLength_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "arnoldiTolerance_: "          << arnoldiTolerance_           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "doubleOrthoTolerance_: "      << doubleOrthoTolerance_       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 24 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";
      // clang-format on

      SaddlePointPreconditioner_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::FGMRESSolver< OperatorType > > SaddlePointFGMRESSolver_;
   std::shared_ptr< SaddlePointSolver< OperatorType > >   SaddlePointPreconditioner_;
   uint_t                                                 maxIterations_;
   real_t                                                 relativeTolerance_;
   real_t                                                 absoluteTolerance_;
   bool                                                   printInfo_;

   uint_t restartLength_;
   real_t arnoldiTolerance_;
   real_t doubleOrthoTolerance_;

   bool lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const SaddlePointFGMRESSolver< OperatorType >& spfgmrs )
{
   return spfgmrs.print( os );
}

} // namespace MantleConvection