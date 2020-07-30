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

/// Concatenates two HyTeG operators to form a new operator.
///
/// Example usage:
/// auto laplace = std::make_shared< P1LaplaceOperatorType >( storage, minLevel, maxLevel );
/// ConcatenatedOperator < P1LaplaceOperatorType, P1LaplaceOperatorType > doubleLaplaceOperator(laplace, laplace);
/// auto solver = std::make_shared< CGSolver< decltype(doubleLaplaceOperator) > >( storage, minLevel, minLevel, max_coarse_iter,
///                                                                                coarse_tolerance );
///
/// \tparam OpType1     Type of the operator applied first.
/// \tparam OpType2     Type of the operator applied second.
template < typename OpType1, typename OpType2 >
class ConcatenatedOperator : public Operator< typename OpType1::srcType, typename OpType2::dstType >
{
 public:
   /// Creates a composite operator.
   /// \param op1  The operator which gets executed first.
   /// \param op2  The operator which gets executed second.
   ConcatenatedOperator( std::shared_ptr< OpType1 > op1, std::shared_ptr< OpType2 > op2 )
   : Operator< typename OpType1::srcType, typename OpType2::dstType >( op1->getStorage(), op1->getMinLevel(), op1->getMaxLevel() )
   , op1_( op1 )
   , op2_( op2 )
   , tmp_( "concatenated_operator_tmp", op1->getStorage(), op1->getMinLevel(), op1->getMaxLevel() )
   {
      checkStaticPreconditions();
   }

   /// Creates a composite operator.
   /// \param op1  The operator which gets executed first.
   /// \param op2  The operator which gets executed second.
   /// \param tmp  Auxiliary vector to execute the matrix vector multiplication.
   ConcatenatedOperator( std::shared_ptr< OpType1 > op1, std::shared_ptr< OpType2 > op2, typename OpType1::dstType& tmp )
   : Operator< typename OpType1::srcType, typename OpType2::dstType >( op1->getStorage(), op1->getMinLevel(), op1->getMaxLevel() )
   , op1_( op1 )
   , op2_( op2 )
   , tmp_( tmp )
   {
      checkStaticPreconditions();
   }

   void apply( const typename OpType1::srcType& src,
               const typename OpType2::dstType& dst,
               size_t                           level,
               DoFType                          flag,
               UpdateType                       updateType = Replace ) const
   {
      WALBERLA_CHECK( updateType == Replace, "Operator concatenation only supported for updateType Replace" );

      op1_->apply( src, tmp_, level, flag, updateType );
      op2_->apply( tmp_, dst, level, flag, updateType );
   }

   std::shared_ptr< typename OpType1::srcType > getDiagonalValues() const
   {
      checkStaticAllTypesEqualPrecondition();
      // TODO: Calculate this with the tmp_ function, op1_.getDiagonalValues() and op2_.getDiagonalValues().
      WALBERLA_ABORT( "ConcatenatedOperator::getDiagonalValues was not implemented yet" );
   };

   std::shared_ptr< typename OpType1::srcType > getInverseDiagonalValues() const
   {
      checkStaticAllTypesEqualPrecondition();
      // TODO: Calculate this with the tmp_ function, op1_.getDiagonalValues() and op2_.getDiagonalValues() and then inverting.
      WALBERLA_ABORT( "ConcatenatedOperator::getDiagonalValues was not implemented yet" );
   };

 private:
   void checkStaticPreconditions()
   {
      static_assert( std::is_same< typename OpType1::dstType, typename OpType2::srcType >::value,
                     "The function output type of the first operator must equal the output type of the second operator" );
   }

   void checkStaticAllTypesEqualPrecondition()
   {
      static_assert( std::is_same< typename OpType1::dstType, typename OpType1::srcType >::value,
                     "All function types must be equal" );
      static_assert( std::is_same< typename OpType2::dstType, typename OpType2::srcType >::value,
                     "All function types must be equal" );
   }

   std::shared_ptr< OpType1 > op1_;
   std::shared_ptr< OpType2 > op2_;
   typename OpType1::dstType  tmp_;
};

} // namespace hyteg
