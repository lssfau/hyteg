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

#include "hyteg/operators/Operator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

namespace MantleConvection {

template < class FineOperatorType, class CoarseOperatorType, class FunctionType = hyteg::P2VectorFunction< real_t > >
class CombinedABlockOperator : public hyteg::Operator< FunctionType, FunctionType >,
                               public hyteg::OperatorWithInverseDiagonal< FunctionType >,
                               public hyteg::WeightedJacobiSmoothable< FunctionType >
{
 public:
   CombinedABlockOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                           size_t                                            minLevel,
                           size_t                                            maxLevel,
                           const std::shared_ptr< FineOperatorType >&        FineOperator,
                           const std::shared_ptr< CoarseOperatorType >&      CoarseOperator )
   : hyteg::Operator< FunctionType, FunctionType >( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , FineOperator_( FineOperator )
   , CoarseOperator_( CoarseOperator )
   {}

   void apply( const FunctionType& src,
               const FunctionType& dst,
               uint_t              level,
               hyteg::DoFType      flag,
               hyteg::UpdateType   updateType = hyteg::Replace ) const override
   {
      if ( CoarseOperator_ != nullptr && level <= CoarseOperator_->getMaxLevel() )
      {
         CoarseOperator_->apply( src, dst, level, flag, updateType );
      }
      else
      {
         if ( FineOperator_ != nullptr )
         {
            FineOperator_->apply( src, dst, level, flag, updateType );
         }
      }
   };

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >&                        mat,
                  const typename FineOperatorType::srcType::template FunctionType< idx_t >& numeratorSrc,
                  const typename FineOperatorType::dstType::template FunctionType< idx_t >& numeratorDst,
                  uint_t                                                                    level,
                  hyteg::DoFType                                                            flag ) const override
   {
      if ( CoarseOperator_ != nullptr && level <= CoarseOperator_->getMaxLevel() )
      {
         CoarseOperator_->toMatrix( mat, numeratorSrc, numeratorDst, level, flag );
      }
      else
      {
         if ( FineOperator_ != nullptr )
         {
            FineOperator_->toMatrix( mat, numeratorSrc, numeratorDst, level, flag );
         }
      }
   }

   std::shared_ptr< FunctionType > getInverseDiagonalValues() const override final { return invDiag_; }

   void computeInverseDiagonalOperatorValues() override final
   {
      if ( invDiag_ == nullptr )
      {
         invDiag_ = std::make_shared< FunctionType >( "combined A diag", storage_, minLevel_, maxLevel_ );
      }
      if ( CoarseOperator_ != nullptr )
      {
         CoarseOperator_->computeInverseDiagonalOperatorValues();
      }
      if ( FineOperator_ != nullptr )
      {
         FineOperator_->computeInverseDiagonalOperatorValues();
      }

      if ( CoarseOperator_ == nullptr && FineOperator_ != nullptr )
      {
         for ( uint_t level = FineOperator_->getMinLevel(); level <= FineOperator_->getMaxLevel(); level++ )
         {
            invDiag_->assign( { real_c( 1 ) }, { *( FineOperator_->getInverseDiagonalValues() ) }, level, hyteg::All );
         }
      }
      else if ( CoarseOperator_ != nullptr && FineOperator_ == nullptr )
      {
         for ( uint_t level = CoarseOperator_->getMinLevel(); level <= CoarseOperator_->getMaxLevel(); level++ )
         {
            invDiag_->assign( { real_c( 1 ) }, { *( CoarseOperator_->getInverseDiagonalValues() ) }, level, hyteg::All );
         }
      }
      else if ( CoarseOperator_ != nullptr && FineOperator_ != nullptr )
      {
         for ( uint_t level = CoarseOperator_->getMinLevel(); level <= FineOperator_->getMaxLevel(); level++ )
         {
            if ( level <= CoarseOperator_->getMaxLevel() )
            {
               invDiag_->assign( { real_c( 1 ) }, { *( CoarseOperator_->getInverseDiagonalValues() ) }, level, hyteg::All );
            }
            else
            {
               invDiag_->assign( { real_c( 1 ) }, { *( FineOperator_->getInverseDiagonalValues() ) }, level, hyteg::All );
            }
         }
      }
   }

   void smooth_jac( const FunctionType& dst,
                    const FunctionType& rhs,
                    const FunctionType& src,
                    real_t              omega,
                    size_t              level,
                    hyteg::DoFType      flag ) const override
   {
      if ( CoarseOperator_ != nullptr && level <= CoarseOperator_->getMaxLevel() )
      {
         if constexpr ( hyteg::SFINAE::has_smooth_jac< CoarseOperatorType >() )
         {
            CoarseOperator_->smooth_jac( dst, rhs, src, omega, level, flag );
         }
         else
         {
            WALBERLA_ABORT( "Coarse operator does not support smooth_jac!" );
         }
      }
      else
      {
         if ( FineOperator_ != nullptr )
         {
            if constexpr ( hyteg::SFINAE::has_smooth_jac< FineOperatorType >() )
            {
               FineOperator_->smooth_jac( dst, rhs, src, omega, level, flag );
            }
            else
            {
               WALBERLA_ABORT( "Fine operator does not support smooth_jac!" );
            }
         }
      }
   }

   std::shared_ptr< FineOperatorType >   getFineOperator() { return FineOperator_; }
   std::shared_ptr< CoarseOperatorType > getCoarseOperator() { return CoarseOperator_; }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< FineOperatorType >   FineOperator_;
   std::shared_ptr< CoarseOperatorType > CoarseOperator_;

   std::shared_ptr< FunctionType > invDiag_;
};

} // namespace MantleConvection
