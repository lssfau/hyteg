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
#include "hyteg/petsc/PETScKSPSolver.hpp"
#include "hyteg/petsc/PETScSolverOptions.hpp"

#include "SaddlePointSolver.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace MantleConvection {

template < class OperatorType >
class SaddlePointPETScCoarseGridSolver : public SaddlePointSolver< OperatorType >
{
 public:
   using SaddlePointSolver< OperatorType >::prefix_;

   SaddlePointPETScCoarseGridSolver( walberla::Config::BlockHandle&                    parameters,
                                     const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                     uint_t                                            level,
                                     PETScSolverOptions&                               solverOptions,
                                     std::string                                       prefix = "" )
   : SaddlePointSolver< OperatorType >( prefix )
   , level_( level )
   {
      maxIterations_ = parameters.getParameter< uint_t >( prefix + std::string( "SaddlePointPETScCoarseGridMaxIterations" ) );
      relativeTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "SaddlePointPETScCoarseGridRelativeTolerance" ) );
      absoluteTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "SaddlePointPETScCoarseGridAbsoluteTolerance" ) );
      printInfo_ = parameters.getParameter< bool >( prefix + std::string( "SaddlePointPETScCoarseGridPrintInfo" ) );
      alwaysReassembleMatrix_ =
          parameters.getParameter< bool >( prefix + std::string( "SaddlePointPETScCoarseGridAlwaysReassembleMatrix" ) );
      zeroInitialGuess_ = parameters.getParameter< bool >( prefix + std::string( "SaddlePointPETScCoarseGridZeroInitialGuess" ) );

      if ( printInfo_ )
      {
         solverOptions.addOption( "-ksp_monitor", "" );
      }

      SaddlePointPETScCoarseGridSolver_ =
          std::make_shared< hyteg::PETScKSPSolver< OperatorType > >( storage,
                                                                     level,
                                                                     (PetscInt) maxIterations_,
                                                                     relativeTolerance_,
                                                                     absoluteTolerance_,
                                                                     solverOptions,
                                                                     alwaysReassembleMatrix_,
                                                                     zeroInitialGuess_,
                                                                     prefix + std::string( "SaddlePointPETScCoarseGrid" ) );
   }

   void reassembleMatrix( const OperatorType& A ) { SaddlePointPETScCoarseGridSolver_->reassembleMatrix( A, level_ ); }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      SaddlePointPETScCoarseGridSolver_->solve( A, x, b, level );
   };

   std::shared_ptr< hyteg::PETScKSPSolver< OperatorType > > getSolver() { return SaddlePointPETScCoarseGridSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "################################################"                                                        << "\n";
      os << std::string( offset, ' ') << "#### Saddle Point PETSc  Coarse Grid Solver ####"                                                        << "\n";
      os << std::string( offset, ' ') << "################################################"                                                        << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 26 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 26 ) << std::left << "maxIterations_: "             << maxIterations_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 26 ) << std::left << "relativeTolerance_: "         << relativeTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 26 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 26 ) << std::left << "printInfo_: "                 << printInfo_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 26 ) << std::left << "alwaysReassembleMatrix_: "    << alwaysReassembleMatrix_     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 26 ) << std::left << "zeroInitialGuess_: "          << zeroInitialGuess_                  ;
      // clang-format on

      return os;
   }

 private:
   uint_t                                                   level_;
   std::shared_ptr< hyteg::PETScKSPSolver< OperatorType > > SaddlePointPETScCoarseGridSolver_;
   uint_t                                                   maxIterations_;
   real_t                                                   relativeTolerance_;
   real_t                                                   absoluteTolerance_;
   bool                                                     printInfo_;
   bool                                                     alwaysReassembleMatrix_;
   bool                                                     zeroInitialGuess_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const SaddlePointPETScCoarseGridSolver< OperatorType >& sppcgs )
{
   return sppcgs.print( os );
}

} // namespace MantleConvection

#endif