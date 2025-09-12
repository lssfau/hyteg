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

#include <type_traits>

#include "core/DataTypes.h"

#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/types/types.hpp"

namespace MantleConvection {

// Additive combination of two operators
// Both should have the same source and destination function type
// Tolerates one or both operators to be excluded (use NoOperator in the template)
// Having a SrcFunctionType and DstFunctionType is a deliberate choice here, since the NoOperator class does not have these types.
// Please do not change this.
template < class OperatorType1, class OperatorType2, class SrcFunctionType, class DstFunctionType >
class AdditiveOperator : public hyteg::Operator< SrcFunctionType, DstFunctionType >
{
 public:
   AdditiveOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                     uint_t                                            minLevel,
                     uint_t                                            maxLevel,
                     std::shared_ptr< OperatorType1 >                  op1,
                     std::shared_ptr< OperatorType2 >                  op2 )
   : hyteg::Operator< SrcFunctionType, DstFunctionType >( storage, minLevel, maxLevel )
   , op1_( op1 )
   , op2_( op2 )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if constexpr ( !std::is_same< OperatorType1, hyteg::NoOperator >::value )
      {
         if ( op1 == nullptr )
         {
            WALBERLA_ABORT( "Operator 1 set to nullptr but type is not NoOperator!" );
         }
      }

      // op2 being nullptr and not of type NoOperator is allowed because it is needed for the Advection Diffusion SPD SubOperator

      checkStaticPreconditions();
   }

   void apply( const SrcFunctionType& src,
               const DstFunctionType& dst,
               uint_t                 level,
               hyteg::DoFType         flag,
               hyteg::UpdateType      updateType = hyteg::Replace ) const override
   {
      // First operator
      if constexpr ( !std::is_same< OperatorType1, hyteg::NoOperator >::value )
      {
         op1_->apply( src, dst, level, flag, updateType );
      }

      // Second operator
      if constexpr ( !std::is_same< OperatorType2, hyteg::NoOperator >::value )
      {
         if constexpr ( firstOperatorApplied )
         {
            if ( op2_ != nullptr )
            {
               op2_->apply( src, dst, level, flag, hyteg::UpdateType::Add );
            }
         }
         else
         {
            if ( op2_ != nullptr )
            {
               op2_->apply( src, dst, level, flag, updateType );
            }
         }
      }
   }

   // dst = alpha * Op1(src) + alpha * Op2(src) + gamma * dst
   void gemv( const real_t&          alpha,
              const SrcFunctionType& src,
              const real_t&          gamma,
              const DstFunctionType& dst,
              uint_t                 level,
              hyteg::DoFType         flag ) const override
   {
      // First operator
      if constexpr ( !std::is_same< OperatorType1, hyteg::NoOperator >::value )
      {
         hyteg::applyGEMV< OperatorType1 >( *op1_, alpha, src, gamma, dst, level, flag );
      }

      // Second operator
      if constexpr ( !std::is_same< OperatorType2, hyteg::NoOperator >::value )
      {
         if constexpr ( firstOperatorApplied )
         {
            if ( op2_ != nullptr )
            {
               hyteg::applyGEMV< OperatorType2 >( *op2_, alpha, src, real_c( 1 ), dst, level, flag );
            }
         }
         else
         {
            if ( op2_ != nullptr )
            {
               hyteg::applyGEMV< OperatorType2 >( *op2_, alpha, src, gamma, dst, level, flag );
            }
         }
      }
   }

   static constexpr bool firstOperatorApplied = ( !std::is_same< OperatorType1, hyteg::NoOperator >::value );

 private:
   void checkStaticPreconditions()
   {
      if constexpr ( !std::is_same< OperatorType1, hyteg::NoOperator >::value )
      {
         static_assert( std::is_same< typename OperatorType1::srcType, SrcFunctionType >::value,
                        "Operator 1 source function type must match the templated source function type." );
         static_assert( std::is_same< typename OperatorType1::dstType, DstFunctionType >::value,
                        "Operator 1 destination function type must match the templated destination function type." );
      }

      if constexpr ( !std::is_same< OperatorType2, hyteg::NoOperator >::value )
      {
         static_assert( std::is_same< typename OperatorType2::srcType, SrcFunctionType >::value,
                        "Operator 2 source function type must match the templated source function type." );
         static_assert( std::is_same< typename OperatorType2::dstType, DstFunctionType >::value,
                        "Operator 2 destination function type must match the templated destination function type." );
      }
   }

   std::shared_ptr< OperatorType1 > op1_;
   std::shared_ptr< OperatorType2 > op2_;
};

} // namespace MantleConvection
