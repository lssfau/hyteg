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

#include "hyteg/misc/SFINAE.hpp"
#include "hyteg/operators/GEMV.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

// Wrapper for a scaled operator using GEMV
// Having a SrcFunctionType and DstFunctionType is a deliberate choice here, since the NoOperator class does not have these types.
// Please do not change this.
template < class OperatorType, class SrcFunctionType, class DstFunctionType >
class ScaledOperator : public hyteg::Operator< SrcFunctionType, DstFunctionType >,
                       public hyteg::OperatorWithInverseDiagonal< SrcFunctionType >,
                       public hyteg::WeightedJacobiSmoothable< SrcFunctionType >
{
 public:
   ScaledOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                   uint_t                                            minLevel,
                   uint_t                                            maxLevel,
                   std::shared_ptr< OperatorType >                   op,
                   typename SrcFunctionType::valueType               scaling )
   : hyteg::Operator< SrcFunctionType, DstFunctionType >( storage, minLevel, maxLevel )
   , op_( op )
   , scaling_( scaling )
   {}

   void apply( const SrcFunctionType& src,
               const DstFunctionType& dst,
               uint_t                 level,
               hyteg::DoFType         flag,
               hyteg::UpdateType      updateType = hyteg::Replace ) const override
   {
      if constexpr ( !std::is_same< OperatorType, hyteg::NoOperator >::value )
      {
         hyteg::applyGEMV< OperatorType >(
             *op_, scaling_, src, ( updateType == hyteg::UpdateType::Replace ) ? real_c( 0 ) : real_c( 1 ), dst, level, flag );
      }
      else
      {
         WALBERLA_ABORT( "Cannot call apply on scaled NoOperator." );
      }
   }

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >&                     mat,
                  const typename SrcFunctionType::template FunctionType< hyteg::idx_t >& src,
                  const typename DstFunctionType::template FunctionType< hyteg::idx_t >& dst,
                  uint_t                                                                 level,
                  hyteg::DoFType                                                         flag ) const override
   {
      if constexpr ( !std::is_same< OperatorType, hyteg::NoOperator >::value )
      {
         hyteg::applyToMatrixScaled( *op_, scaling_, mat, src, dst, level, flag );
      }
      else
      {
         WALBERLA_ABORT( "Cannot call toMatrix on scaled NoOperator." );
      }
   }

   std::shared_ptr< SrcFunctionType > getInverseDiagonalValues() const override
   {
      if constexpr ( !std::is_same< OperatorType, hyteg::NoOperator >::value )
      {
         if constexpr ( hyteg::SFINAE::has_getInverseDiagonalValues< OperatorType >() )
         {
            return op_->getInverseDiagonalValues();
         }

         WALBERLA_ABORT(
             "Scaled operator does not support getInverseDiagonalValues. Maybe the operator does not have a matching source and destination type?" );
         return nullptr;
      }
      else
      {
         WALBERLA_ABORT( "Cannot call getInverseDiagonalValues on scaled NoOperator." );
         return nullptr;
      }
   }

   void computeInverseDiagonalOperatorValues() override
   {
      if constexpr ( !std::is_same< OperatorType, hyteg::NoOperator >::value )
      {
         if constexpr ( hyteg::SFINAE::has_computeInverseDiagonalOperatorValues< OperatorType >() )
         {
            hyteg::applyComputeInverseDiagonalOperatorValuesScaled( *op_, scaling_ );
         }
         else
         {
            WALBERLA_ABORT(
                "Scaled operator does not support computeInverseDiagonalOperatorValues. Maybe the operator does not have a matching source and destination type?" );
         }
      }
      else
      {
         WALBERLA_ABORT( "Cannot call computeInverseDiagonalOperatorValues on scaled NoOperator." );
      }
   }

   void smooth_jac( const SrcFunctionType& dst,
                    const SrcFunctionType& rhs,
                    const SrcFunctionType& src,
                    real_t                 omega,
                    size_t                 level,
                    DoFType                flag ) const override
   {
      if constexpr ( !std::is_same< OperatorType, hyteg::NoOperator >::value )
      {
         if constexpr ( hyteg::SFINAE::has_smooth_jac< OperatorType >() )
         {
            hyteg::applySmoothJacScaled( *op_, scaling_, dst, rhs, src, omega, level, flag );
         }
         else
         {
            WALBERLA_ABORT(
                "Scaled operator does not support smooth_jac. Maybe the operator does not have a matching source and destination type?" );
         }
      }
      else
      {
         WALBERLA_ABORT( "Cannot call smooth_jac on scaled NoOperator." );
      }
   }

   std::shared_ptr< OperatorType >     getOperatorPtr() { return op_; }
   typename SrcFunctionType::valueType getScaling() { return scaling_; }

   void setScaling( typename SrcFunctionType::valueType scaling ) { scaling_ = scaling; }

 private:
   std::shared_ptr< OperatorType >     op_;
   typename SrcFunctionType::valueType scaling_;
};

} // namespace hyteg
