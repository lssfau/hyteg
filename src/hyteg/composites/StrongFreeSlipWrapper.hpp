/*
 * Copyright (c) 2020 Andreas Wagner
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

#include <type_traits>

#include "hyteg/Operator.hpp"

namespace hyteg {

/// Concatenates a HyTeG operators with the strong free-slip boundary operator.
///
/// Usage Example:
/// auto laplace = std::make_shared< P1ConstantLaplaceOperator >( storage, level, level );
/// auto projection = std::make_shared< P1ProjectNormalOperator > ( storage, level, level, normal );
/// StrongFreeSlipWrapper< P1ConstantLaplaceOperator > wrapper( laplace, projection );
/// auto solver = std::make_shared< CGSolver< decltype(wrapper) > >( storage, minLevel, minLevel, max_coarse_iter, coarse_tolerance );
///
template < typename OpType, typename ProjOpType >
class StrongFreeSlipWrapper : public Operator< typename OpType::srcType, typename OpType::dstType >
{
 public:
   StrongFreeSlipWrapper( std::shared_ptr< OpType > op, std::shared_ptr< ProjOpType > projOp, DoFType projFlag )
   : Operator< typename OpType::srcType, typename OpType::dstType >( op->getStorage(), op->getMinLevel(), op->getMaxLevel() )
   , op_( op )
   , projOp_( projOp )
   , projFlag_( projFlag )
   , diagonalValues_( nullptr )
   , inverseDiagonalValues_( nullptr )
   {}

   void apply( const typename OpType::srcType& src,
               const typename OpType::dstType& dst,
               size_t                          level,
               DoFType                         flag,
               UpdateType                      updateType = Replace ) const
   {
      WALBERLA_CHECK( updateType == Replace, "Operator concatenation only supported for updateType Replace" );

      op_->apply( src, dst, level, flag );
      projOp_->apply( dst, level, projFlag_ );
   }

   std::shared_ptr< typename OpType::srcType > getDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          diagonalValues_,
          "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return diagonalValues_;
   };

   std::shared_ptr< typename OpType::srcType > getInverseDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          inverseDiagonalValues_,
          "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return inverseDiagonalValues_;
   };

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { WALBERLA_ABORT( "computeDiagonalOperatorValues is not implemented yet." ) }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues()
   {
      WALBERLA_ABORT( "computeInverseDiagonalOperatorValues is not implemented yet." )
   }

 private:
   std::shared_ptr< OpType > op_;
   std::shared_ptr< ProjOpType > projOp_;

   DoFType projFlag_;

   std::shared_ptr< typename OpType::dstType > diagonalValues_;
   std::shared_ptr< typename OpType::dstType > inverseDiagonalValues_;
};

} // namespace hyteg
