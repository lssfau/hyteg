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
#include "hyteg/solvers/FunctionMultiplicationPreconditioner.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "../../OperatorTools/SchurOperator.hpp"
#include "SchurIdentityPreconditioner.hpp"
#include "SchurSolver.hpp"

namespace MantleConvection {

template < class MassOperatorType >
class SchurCGMassSolver : public SchurSolver< SchurOperator< typename MassOperatorType::dstType > >
{
 public:
   using SchurSolver< SchurOperator< typename MassOperatorType::dstType > >::prefix_;
   typedef typename MassOperatorType::dstType PressureFunctionType;

   SchurCGMassSolver( walberla::Config::BlockHandle&                                                 parameters,
                      const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
                      uint_t                                                                         minLevel,
                      uint_t                                                                         maxLevel,
                      const std::shared_ptr< MassOperatorType >&                                     MassOperator,
                      bool                                                                           lowMemoryMode = false,
                      const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& schurPreconditioner =
                          std::make_shared< SchurIdentityPreconditioner< PressureFunctionType > >(),
                      bool        computeInverse = true,
                      std::string prefix         = "" )
   : SchurSolver< SchurOperator< typename MassOperatorType::dstType > >( prefix )
   , schurPreconditioner_( schurPreconditioner )
   , MassOperator_( MassOperator )
   , lowMemoryMode_( lowMemoryMode )
   {
      if ( computeInverse )
      {
         MassOperator_->computeInverseDiagonalOperatorValues();
      }

      maxIterations_     = parameters.getParameter< uint_t >( prefix + std::string( "SchurCGMassMaxIterations" ) );
      relativeTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SchurCGMassRelativeTolerance" ) );
      absoluteTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SchurCGMassAbsoluteTolerance" ) );
      printInfo_         = parameters.getParameter< bool >( prefix + std::string( "SchurCGMassPrintInfo" ) );

      auto substOperator_ = std::make_shared< SchurOperator< PressureFunctionType > >( storage, minLevel, maxLevel );

      auto substSolver_ =
          std::make_shared< hyteg::SubstitutePreconditioner< MassOperatorType, SchurOperator< PressureFunctionType > > >(
              schurPreconditioner, substOperator_ );

      SchurCGMassSolver_ = std::make_shared< hyteg::CGSolver< MassOperatorType > >(
          storage, minLevel, maxLevel, maxIterations_, relativeTolerance_, absoluteTolerance_, substSolver_, lowMemoryMode_ );
      SchurCGMassSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         SchurCGMassSolver_->setName( std::string( "CG Schur CG Mass " ) + prefix_ );
      }
      else
      {
         SchurCGMassSolver_->setName( "CG Schur CG Mass" );
      }
   }

   void solve( const SchurOperator< typename MassOperatorType::dstType >& A,
               const typename MassOperatorType::dstType&                  x,
               const typename MassOperatorType::dstType&                  b,
               const walberla::uint_t                                     level ) override
   {
      WALBERLA_UNUSED( A );
      SchurCGMassSolver_->solve( *MassOperator_, x, b, level );
   };

   std::shared_ptr< hyteg::CGSolver< MassOperatorType > > getSolver() { return SchurCGMassSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "################################################"                                             << "\n";
      os << std::string( offset, ' ') << "############# Schur CG Mass Solver #############"                                             << "\n";
      os << std::string( offset, ' ') << "################################################"                                             << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "              << prefix_                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxIterations_: "       << maxIterations_         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relativeTolerance_: "   << relativeTolerance_     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "absoluteTolerance_: "   << absoluteTolerance_     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "printInfo_: "           << printInfo_             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "lowMemoryMode_: "       << lowMemoryMode_         << "\n";
      // clang-format on

      schurPreconditioner_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::CGSolver< MassOperatorType > >                  SchurCGMassSolver_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > schurPreconditioner_;
   std::shared_ptr< MassOperatorType >                                     MassOperator_;

   uint_t maxIterations_;
   real_t relativeTolerance_;
   real_t absoluteTolerance_;
   bool   printInfo_;

   bool lowMemoryMode_;
};

template < class MassOperatorType >
inline std::ostream& operator<<( std::ostream& os, const SchurCGMassSolver< MassOperatorType >& scgms )
{
   return scgms.print( os );
}

} // namespace MantleConvection