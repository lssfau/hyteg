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

#include "AdvectionDiffusionIdentityPreconditioner.hpp"
#include "AdvectionDiffusionSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class AdvectionDiffusionFGMRESSolver : public AdvectionDiffusionSolver< OperatorType >
{
 public:
   using AdvectionDiffusionSolver< OperatorType >::prefix_;

   AdvectionDiffusionFGMRESSolver(
       walberla::Config::BlockHandle&                                     parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                  storage,
       uint_t                                                             minLevel,
       uint_t                                                             maxLevel,
       bool                                                               lowMemoryMode = false,
       const std::shared_ptr< AdvectionDiffusionSolver< OperatorType > >& AdvectionDiffusionPreconditioner =
           std::make_shared< AdvectionDiffusionIdentityPreconditioner< OperatorType > >(),
       std::string prefix = "" )
   : AdvectionDiffusionSolver< OperatorType >( prefix )
   , AdvectionDiffusionPreconditioner_( AdvectionDiffusionPreconditioner )
   , lowMemoryMode_( lowMemoryMode )
   {
      maxIterations_ = parameters.getParameter< uint_t >( prefix + std::string( "AdvectionDiffusionFGMRESMaxIterations" ) );
      relativeTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "AdvectionDiffusionFGMRESRelativeTolerance" ) );
      absoluteTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "AdvectionDiffusionFGMRESAbsoluteTolerance" ) );
      printInfo_ = parameters.getParameter< bool >( prefix + std::string( "AdvectionDiffusionFGMRESPrintInfo" ) );

      restartLength_    = parameters.getParameter< uint_t >( prefix + std::string( "AdvectionDiffusionFGMRESRestartLength" ) );
      arnoldiTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "AdvectionDiffusionFGMRESArnoldiTolerance" ) );
      doubleOrthoTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "AdvectionDiffusionFGMRESDoubleOrthoTolerance" ) );

      AdvectionDiffusionFGMRESSolver_ =
          std::make_shared< hyteg::FGMRESSolver< OperatorType > >( storage,
                                                                   minLevel,
                                                                   maxLevel,
                                                                   maxIterations_,
                                                                   relativeTolerance_,
                                                                   absoluteTolerance_,
                                                                   AdvectionDiffusionPreconditioner_,
                                                                   restartLength_,
                                                                   arnoldiTolerance_,
                                                                   doubleOrthoTolerance_,
                                                                   lowMemoryMode_ );
      AdvectionDiffusionFGMRESSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         AdvectionDiffusionFGMRESSolver_->setName( std::string( "FGMRES Advection Diffusion prefix " ) + prefix_ );
      }
      else
      {
         AdvectionDiffusionFGMRESSolver_->setName( "FGMRES Advection Diffusion" );
      }
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      AdvectionDiffusionFGMRESSolver_->solve( A, x, b, level );
   };

   std::shared_ptr< hyteg::FGMRESSolver< OperatorType > > getSolver() { return AdvectionDiffusionFGMRESSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "###############################################"                                                         << "\n";
      os << std::string( offset, ' ') << "###### Advection Diffusion FGMRES Solver ######"                                                         << "\n";
      os << std::string( offset, ' ') << "###############################################"                                                         << "\n";
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

      AdvectionDiffusionPreconditioner_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::FGMRESSolver< OperatorType > >      AdvectionDiffusionFGMRESSolver_;
   std::shared_ptr< AdvectionDiffusionSolver< OperatorType > > AdvectionDiffusionPreconditioner_;
   uint_t                                                      maxIterations_;
   real_t                                                      relativeTolerance_;
   real_t                                                      absoluteTolerance_;
   bool                                                        printInfo_;

   uint_t restartLength_;
   real_t arnoldiTolerance_;
   real_t doubleOrthoTolerance_;

   bool lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const AdvectionDiffusionFGMRESSolver< OperatorType >& adfgmrs )
{
   return adfgmrs.print( os );
}

} // namespace MantleConvection