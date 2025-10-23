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
#include "hyteg/solvers/SubstitutePreconditioner.hpp"

#include "../../OperatorTools/ApproximateSchurComplementOperator.hpp"
#include "../../OperatorTools/SchurOperator.hpp"
#include "../ABlock/ABlockSolver.hpp"
#include "SchurIdentityPreconditioner.hpp"
#include "SchurSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class SchurCGCoarseGridSolver : public SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >
{
 public:
   using SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >::prefix_;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef ApproximateSchurComplementOperator< OperatorType > ApproximateSchurOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   SchurCGCoarseGridSolver( walberla::Config::BlockHandle&                                                 parameters,
                            const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
                            uint_t                                                                         minLevel,
                            uint_t                                                                         maxLevel,
                            const std::shared_ptr< OperatorType >&                                         saddleOp,
                            const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >&      ABlockSolver,
                            hyteg::BoundaryCondition                                                       velocityBC,
                            bool                                                                           lowMemoryMode = false,
                            const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& schurPreconditioner =
                                std::make_shared< SchurIdentityPreconditioner< PressureFunctionType > >(),
                            const std::shared_ptr< ApproximateSchurOperatorType > approximateSchurOperator = nullptr,
                            std::string                                           prefix                   = "" )
   : SchurCGCoarseGridSolver( parameters,
                              storage,
                              minLevel,
                              maxLevel,
                              saddleOp,
                              ABlockSolver,
                              velocityBC,
                              velocityBC,
                              velocityBC,
                              lowMemoryMode,
                              schurPreconditioner,
                              approximateSchurOperator,
                              prefix )
   {}

   SchurCGCoarseGridSolver( walberla::Config::BlockHandle&                                                 parameters,
                            const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
                            uint_t                                                                         minLevel,
                            uint_t                                                                         maxLevel,
                            const std::shared_ptr< OperatorType >&                                         saddleOp,
                            const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >&      ABlockSolver,
                            hyteg::BoundaryCondition                                                       velocityBCx,
                            hyteg::BoundaryCondition                                                       velocityBCy,
                            hyteg::BoundaryCondition                                                       velocityBCz,
                            bool                                                                           lowMemoryMode = false,
                            const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& schurPreconditioner =
                                std::make_shared< SchurIdentityPreconditioner< SchurOperator< PressureFunctionType > > >(),
                            const std::shared_ptr< ApproximateSchurOperatorType > approximateSchurOperator = nullptr,
                            std::string                                           prefix                   = "" )
   : SchurSolver< SchurOperator< PressureFunctionType > >( prefix )
   , saddleOp_( saddleOp )
   , schurPreconditioner_( schurPreconditioner )
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

      maxIterations_     = parameters.getParameter< uint_t >( prefix + std::string( "SchurCGCoarseGridMaxIterations" ) );
      relativeTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SchurCGCoarseGridRelativeTolerance" ) );
      absoluteTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SchurCGCoarseGridAbsoluteTolerance" ) );
      printInfo_         = parameters.getParameter< bool >( prefix + std::string( "SchurCGCoarseGridPrintInfo" ) );

      substOperator_ = std::make_shared< SchurOperator< PressureFunctionType > >( storage, minLevel, maxLevel );

      substSolver_ = std::make_shared<
          hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >(
          schurPreconditioner_, substOperator_ );

      SchurCGCoarseGridSolver_ = std::make_shared< hyteg::CGSolver< ApproximateSchurOperatorType > >(
          storage, minLevel, maxLevel, maxIterations_, relativeTolerance_, absoluteTolerance_, substSolver_, lowMemoryMode_ );
      SchurCGCoarseGridSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         SchurCGCoarseGridSolver_->setName( std::string( "CG Schur Coarse Grid prefix " ) + prefix_ );
      }
      else
      {
         SchurCGCoarseGridSolver_->setName( "CG Schur Coarse Grid" );
      }
   }

   void solve( const SchurOperator< PressureFunctionType >& A,
               const PressureFunctionType&                  x,
               const PressureFunctionType&                  b,
               const walberla::uint_t                       level ) override
   {
      SchurCGCoarseGridSolver_->solve( *approximateSchurOperator_, x, b, level );
   };

   std::shared_ptr< hyteg::CGSolver< ApproximateSchurOperatorType > > getSolver() { return SchurCGCoarseGridSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "###############################################"                                                       << "\n";
      os << std::string( offset, ' ') << "######### Schur CG Coarse Grid Solver #########"                                                       << "\n";
      os << std::string( offset, ' ') << "###############################################"                                                       << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "                    << prefix_                   << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxIterations_: "             << maxIterations_            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relativeTolerance_: "         << relativeTolerance_        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "printInfo_: "                 << printInfo_                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_            << "\n";
      // clang-format on

      schurPreconditioner_->print( os, offset + 3 );

      return os;
   }

   SchurOperator< PressureFunctionType >&                   getSchurOp() { return *substOperator_; };
   std::shared_ptr< SchurOperator< PressureFunctionType > > getSchurOpPtr() { return substOperator_; };
   ApproximateSchurOperatorType&                            getApproximateSchurOp() { return *approximateSchurOperator_; };
   std::shared_ptr< ApproximateSchurOperatorType >          getApproximateSchurOpPtr() { return approximateSchurOperator_; };

 private:
   std::shared_ptr< OperatorType >                                         saddleOp_;
   std::shared_ptr< ApproximateSchurOperatorType >                         approximateSchurOperator_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > schurPreconditioner_;
   std::shared_ptr< hyteg::CGSolver< ApproximateSchurOperatorType > >      SchurCGCoarseGridSolver_;
   std::shared_ptr< SchurOperator< PressureFunctionType > >                substOperator_;
   std::shared_ptr< hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >
       substSolver_;

   uint_t maxIterations_;
   real_t relativeTolerance_;
   real_t absoluteTolerance_;
   bool   printInfo_;
   bool   lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const SchurCGCoarseGridSolver< OperatorType >& scgcgs )
{
   return scgcgs.print( os );
}

} // namespace MantleConvection