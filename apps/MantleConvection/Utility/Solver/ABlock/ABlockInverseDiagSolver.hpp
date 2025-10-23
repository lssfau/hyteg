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
#include "hyteg/solvers/FunctionMultiplicationPreconditioner.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class ABlockInverseDiagSolver : public ABlockSolver< OperatorType >
{
 public:
   using ABlockSolver< OperatorType >::prefix_;

   ABlockInverseDiagSolver( walberla::Config::BlockHandle&                    parameters,
                            const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                            uint_t                                            minLevel,
                            uint_t                                            maxLevel,
                            OperatorType&                                     ABlock,
                            bool                                              computeInverse = true,
                            std::string                                       prefix         = "" )
   : ABlockSolver< OperatorType >( prefix )
   {
      WALBERLA_UNUSED( parameters );
      if ( computeInverse )
      {
         ABlock.computeInverseDiagonalOperatorValues();
      }

      funcMult_ = std::make_shared< hyteg::FunctionMultiplicationPreconditioner< OperatorType > >(
          storage, minLevel, maxLevel, ABlock.getInverseDiagonalValues() );
   }

   template < class OtherOperatorType >
   void setInvDiag( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                    uint_t                                            minLevel,
                    uint_t                                            maxLevel,
                    OtherOperatorType&                                Op )
   {
      funcMult_ = std::make_shared< hyteg::FunctionMultiplicationPreconditioner< OperatorType > >(
          storage, minLevel, maxLevel, Op.getInverseDiagonalValues() );
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      funcMult_->solve( A, x, b, level );
   };

   std::shared_ptr< hyteg::FunctionMultiplicationPreconditioner< OperatorType > > getSolver() { return funcMult_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "###############################################"                                                        << "\n";
      os << std::string( offset, ' ') << "######### A Block Inverse Diag Solver #########"                                                        << "\n";
      os << std::string( offset, ' ') << "###############################################"                                                        << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "prefix_: "  << prefix_                                             ;
      // clang-format on

      return os;
   }

 private:
   std::shared_ptr< hyteg::FunctionMultiplicationPreconditioner< OperatorType > > funcMult_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const ABlockInverseDiagSolver< OperatorType >& abids )
{
   return abids.print( os );
}

} // namespace MantleConvection