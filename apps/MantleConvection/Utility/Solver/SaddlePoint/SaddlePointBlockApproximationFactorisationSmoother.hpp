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
#include "hyteg/solvers/Solver.hpp"

#include "../../OperatorTools/BlockPreconditioners.hpp"
#include "SaddlePointSolver.hpp"

namespace MantleConvection {

template < class OperatorType,
           bool projectASolverRHS_     = true,
           bool projectSchurSolverRHS_ = true,
           bool postProjectPressure_   = true,
           bool postProjectVelocity_   = true >
class SaddlePointBlockApproximationFactorisationSmoother : public SaddlePointSolver< OperatorType >
{
 public:
   using SaddlePointSolver< OperatorType >::prefix_;
   typedef typename OperatorType::srcType                     SrcFunctionType;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   SaddlePointBlockApproximationFactorisationSmoother(
       walberla::Config::BlockHandle&                                                 parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       uint_t                                                                         minLevel,
       uint_t                                                                         maxLevel,
       const std::shared_ptr< OperatorType >&                                         SaddleOp,
       const std::shared_ptr< ABlockSolver< AOperatorType > >&                        ABlockSolver,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& SchurComplementSolver,
       bool                                                                           lowMemoryMode = false,
       hyteg::DoFType flag   = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
       std::string    prefix = "" )
   : SaddlePointSolver< OperatorType >( prefix )
   , ABlockSolver_( ABlockSolver )
   , SchurComplementSolver_( SchurComplementSolver )
   , lowMemoryMode_( lowMemoryMode )
   {
      VelocityIterations_ = parameters.getParameter< uint_t >(
          prefix + std::string( "SaddlePointBlockApproximationFactorisationSmootherVelocityIterations" ) );
      relaxParamA_ = parameters.getParameter< real_t >(
          prefix + std::string( "SaddlePointBlockApproximationFactorisationSmootherRelaxParamA" ) );
      relaxParamSchur_ = parameters.getParameter< real_t >(
          prefix + std::string( "SaddlePointBlockApproximationFactorisationSmootherRelaxParamSchur" ) );

      uzawaSmoother_ = std::make_shared< MantleConvection::BlockApproximateFactorisationPreconditioner< OperatorType,
                                                                                                        projectASolverRHS_,
                                                                                                        projectSchurSolverRHS_,
                                                                                                        postProjectPressure_,
                                                                                                        postProjectVelocity_ > >(
          storage,
          minLevel,
          maxLevel,
          SaddleOp,
          ABlockSolver_,
          SchurComplementSolver_,
          relaxParamA_,
          relaxParamSchur_,
          lowMemoryMode_,
          VelocityIterations_,
          flag );
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      uzawaSmoother_->solve( A, x, b, level );
   };

   std::shared_ptr< MantleConvection::BlockApproximateFactorisationPreconditioner< OperatorType,
                                                                                   projectASolverRHS_,
                                                                                   projectSchurSolverRHS_,
                                                                                   postProjectPressure_,
                                                                                   postProjectVelocity_ > >
       getSolver()
   {
      return uzawaSmoother_;
   }   

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#################################################################"                                       << "\n";
      os << std::string( offset, ' ') << "#### Saddle Point Block Approximation Factorisation Smoother ####"                                       << "\n";
      os << std::string( offset, ' ') << "#################################################################"                                       << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "relaxParamA_: "               << relaxParamA_                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "relaxParamSchur_: "           << relaxParamSchur_            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "VelocityIterations_: "        << VelocityIterations_         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";      
      // clang-format on

      ABlockSolver_->print( os, offset + 3 );
      os << "\n";
      SchurComplementSolver_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< ABlockSolver< AOperatorType > >                        ABlockSolver_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurComplementSolver_;
   std::shared_ptr< MantleConvection::BlockApproximateFactorisationPreconditioner< OperatorType,
                                                                                   projectASolverRHS_,
                                                                                   projectSchurSolverRHS_,
                                                                                   postProjectPressure_,
                                                                                   postProjectVelocity_ > >
       uzawaSmoother_;

   real_t relaxParamA_;
   real_t relaxParamSchur_;
   uint_t VelocityIterations_;

   bool lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream&                                                             os,
                                 const SaddlePointBlockApproximationFactorisationSmoother< OperatorType >& spbafp )
{
   return spbafp.print( os );
}

} // namespace MantleConvection