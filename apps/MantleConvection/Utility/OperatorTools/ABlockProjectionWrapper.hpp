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

#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/misc/SFINAE.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

namespace MantleConvection {

template < typename OperatorType,
           typename VelocityProjectionOperatorType_ = hyteg::NoOperator,
           bool preProjection                       = true,
           bool postProjection                      = true,
           bool allowPreProjectionToChangeSrc       = true,
           bool allowPostProjectionToChangeDst      = true,
           bool includeProjectionsToMatrix          = false >
class ABlockProjectionWrapper : public hyteg::Operator< typename OperatorType::srcType, typename OperatorType::dstType >,
                                public hyteg::OperatorWithInverseDiagonal< typename OperatorType::srcType >,
                                public hyteg::WeightedJacobiSmoothable< typename OperatorType::srcType >
{
 public:
   typedef ABlockProjectionWrapper< OperatorType,
                                    VelocityProjectionOperatorType_,
                                    preProjection,
                                    postProjection,
                                    allowPreProjectionToChangeSrc,
                                    allowPostProjectionToChangeDst,
                                    includeProjectionsToMatrix >
       AOperatorType;

   typedef VelocityProjectionOperatorType_ VelocityProjectionOperatorType;

   ABlockProjectionWrapper( std::shared_ptr< OperatorType >                    op,
                            std::shared_ptr< VelocityProjectionOperatorType_ > projection     = nullptr,
                            hyteg::DoFType                                     projectionFlag = hyteg::FreeslipBoundary,
                            bool                                               lowMemoryMode  = false )
   : hyteg::Operator< typename OperatorType::srcType, typename OperatorType::dstType >( op->getStorage(),
                                                                                        op->getMinLevel(),
                                                                                        op->getMaxLevel() )
   , op_( op )
   , projection_( projection )
   , projectionFlag_( projectionFlag )
   , lowMemoryMode_( lowMemoryMode )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if ( op == nullptr )
      {
         WALBERLA_ABORT( "Operator is nullptr!" );
      }

      if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
      {
         if ( projection == nullptr )
         {
            WALBERLA_ABORT( "Pre projection activated and not of type NoOperator but projection is nullptr!" );
         }
      }

      if constexpr ( ( postProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
      {
         if ( projection == nullptr )
         {
            WALBERLA_ABORT( "Post projection activated and not of type NoOperator but projection is nullptr!" );
         }
      }

      // ##################################
      // ###### Define tmp functions ######
      // ##################################

      if constexpr ( useTmpSrc )
      {
         if ( !lowMemoryMode_ )
         {
            tmpSrc_ = std::make_shared< typename OperatorType::srcType >(
                "ABlockProjectionWrapper tmpSrc", op->getStorage(), op->getMinLevel(), op->getMaxLevel() );
         }
      }
      if constexpr ( useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmpDst_ = std::make_shared< typename OperatorType::dstType >(
                "ABlockProjectionWrapper tmpDst", op->getStorage(), op->getMinLevel(), op->getMaxLevel() );
         }
      }
   }

   void apply( const typename OperatorType::srcType& src,
               const typename OperatorType::dstType& dst,
               uint_t                                level,
               hyteg::DoFType                        flag,
               hyteg::UpdateType                     updateType = hyteg::Replace ) const override
   {
      std::shared_ptr< typename OperatorType::srcType > tmpSrcApply;
      std::shared_ptr< typename OperatorType::dstType > tmpDstApply;

      // handle low memory mode
      if constexpr ( useTmpSrc )
      {
         if ( !lowMemoryMode_ )
         {
            tmpSrcApply = tmpSrc_;
         }
         else
         {
            tmpSrcApply = hyteg::getTemporaryFunction< typename OperatorType::srcType >(
                op_->getStorage(), op_->getMinLevel(), op_->getMaxLevel() );
         }
      }

      if constexpr ( useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmpDstApply = tmpDst_;
         }
         else
         {
            tmpDstApply = hyteg::getTemporaryFunction< typename OperatorType::dstType >(
                op_->getStorage(), op_->getMinLevel(), op_->getMaxLevel() );
         }
      }

      // prepare tmp usage
      if constexpr ( useTmpSrc )
      {
         tmpSrcApply->copyBoundaryConditionFromFunction( src );
         tmpSrcApply->assign( { real_c( 1 ) }, { src }, level, hyteg::All );
      }

      if constexpr ( useTmpDst )
      {
         tmpDstApply->copyBoundaryConditionFromFunction( dst );
      }

      // pre project
      if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
      {
         if constexpr ( allowPreProjectionToChangeSrc )
         {
            projection_->project( src, level, projectionFlag_ );
         }
         else
         {
            projection_->project( *tmpSrcApply, level, projectionFlag_ );
         }
      }

      // apply operator
      if constexpr ( !useTmpSrc && !useTmpDst )
      {
         op_->apply( src, dst, level, flag, updateType );
      }
      else if constexpr ( useTmpSrc && !useTmpDst )
      {
         op_->apply( *tmpSrcApply, dst, level, flag, updateType );
      }
      else if constexpr ( !useTmpSrc && useTmpDst )
      {
         op_->apply( src, *tmpDstApply, level, flag, hyteg::Replace );
      }
      else
      {
         op_->apply( *tmpSrcApply, *tmpDstApply, level, flag, hyteg::Replace );
      }

      // post project
      if constexpr ( ( postProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
      {
         if constexpr ( allowPostProjectionToChangeDst )
         {
            projection_->project( dst, level, projectionFlag_ );
         }
         else
         {
            projection_->project( *tmpDstApply, level, projectionFlag_ );
            if ( updateType == hyteg::Replace )
            {
               dst.assign( { real_c( 1 ) }, { *tmpDstApply }, level, flag );
            }
            else
            {
               dst.assign( { real_c( 1 ), real_c( 1 ) }, { *tmpDstApply, dst }, level, flag );
            }
         }
      }
   }

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >&                    mat,
                  const typename OperatorType::srcType::template FunctionType< idx_t >& numeratorSrc,
                  const typename OperatorType::dstType::template FunctionType< idx_t >& numeratorDst,
                  uint_t                                                                level,
                  hyteg::DoFType                                                        flag ) const override
   {
      if constexpr ( includeProjectionsToMatrix )
      {
         auto matProxyOp = mat->createCopy();
         op_->toMatrix( matProxyOp, numeratorSrc, numeratorDst, level, flag );

         std::shared_ptr< hyteg::SparseMatrixProxy > matProxyProjectionPre;
         if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
         {
            matProxyProjectionPre = mat->createEmptyCopy();
            projection_->toMatrix( matProxyProjectionPre, numeratorSrc, numeratorDst, level, projectionFlag_ );
         }

         std::shared_ptr< hyteg::SparseMatrixProxy > matProxyProjectionPost;
         if constexpr ( ( postProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
         {
            matProxyProjectionPost = mat->createEmptyCopy();
            projection_->toMatrix( matProxyProjectionPost, numeratorSrc, numeratorDst, level, projectionFlag_ );
         }

         std::vector< std::shared_ptr< hyteg::SparseMatrixProxy > > matrices;

         if constexpr ( ( postProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
         {
            matrices.push_back( matProxyProjectionPost );
         }

         matrices.push_back( matProxyOp );

         if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
         {
            matrices.push_back( matProxyProjectionPre );
         }

         mat->createFromMatrixProduct( matrices );
      }
      else
      {
         op_->toMatrix( mat, numeratorSrc, numeratorDst, level, flag );
      }
   }

   std::shared_ptr< typename OperatorType::srcType > getInverseDiagonalValues() const override
   {
      if constexpr ( hyteg::SFINAE::has_getInverseDiagonalValues< OperatorType >() )
      {
         // TODO: Involve projections into the diagonal computation.
         return op_->getInverseDiagonalValues();
      }
      else
      {
         WALBERLA_ABORT( "getInverseDiagonalValues not supported by operator!" );
         return nullptr;
      }
   }

   void computeInverseDiagonalOperatorValues() override
   {
      // TODO: Involve projections into the diagonal computation.
      if constexpr ( hyteg::SFINAE::has_computeInverseDiagonalOperatorValues< OperatorType >() )
      {
         op_->computeInverseDiagonalOperatorValues();
      }
      else
      {
         WALBERLA_ABORT( "computeInverseDiagonalOperatorValues not supported by operator!" );
      }
   }

   std::shared_ptr< typename OperatorType::srcType > getLumpedInverseDiagonalValues() const
   {
      if constexpr ( hyteg::SFINAE::has_getLumpedInverseDiagonalValues< OperatorType >() )
      {
         // TODO: Involve projections into the diagonal computation.
         return op_->getLumpedInverseDiagonalValues();
      }
      else
      {
         WALBERLA_ABORT( "getLumpedInverseDiagonalValues not supported by operator!" );
         return nullptr;
      }
   }

   void computeLumpedInverseDiagonalOperatorValues()
   {
      // TODO: Involve projections into the diagonal computation.
      if constexpr ( hyteg::SFINAE::has_computeLumpedInverseDiagonalOperatorValues< OperatorType >() )
      {
         op_->computeLumpedInverseDiagonalOperatorValues();
      }
      else
      {
         WALBERLA_ABORT( "computeLumpedInverseDiagonalOperatorValues not supported by operator!" );
      }
   }

   void smooth_jac( const typename OperatorType::srcType& dst,
                    const typename OperatorType::srcType& rhs,
                    const typename OperatorType::srcType& src,
                    real_t                                omega,
                    size_t                                level,
                    hyteg::DoFType                        flag ) const override
   {
      if constexpr ( hyteg::SFINAE::has_smooth_jac< OperatorType >() )
      {
         op_->smooth_jac( dst, rhs, src, omega, level, flag );
         // post project
         if constexpr ( ( postProjection ) && ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) )
         {
            projection_->project( dst, level, projectionFlag_ );
         }
      }
      else
      {
         WALBERLA_ABORT( "smooth_jac not supported by operator!" );
      }
   }

   std::shared_ptr< OperatorType > getOperatorPtr() { return op_; }

   const AOperatorType&                               getA() const { return *this; }
   std::shared_ptr< VelocityProjectionOperatorType_ > getProjPtr() const { return projection_; }
   hyteg::DoFType                                     getProjFlag() const { return projectionFlag_; }

   static constexpr bool useTmpSrc = ( !allowPreProjectionToChangeSrc ) && ( preProjection ) &&
                                     ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value );
   static constexpr bool useTmpDst = ( !allowPostProjectionToChangeDst ) && ( postProjection ) &&
                                     ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value );

 private:
   std::shared_ptr< OperatorType >                    op_;
   std::shared_ptr< VelocityProjectionOperatorType_ > projection_;

   hyteg::DoFType projectionFlag_;

   std::shared_ptr< typename OperatorType::srcType > tmpSrc_;
   std::shared_ptr< typename OperatorType::dstType > tmpDst_;

   bool lowMemoryMode_;
};

} // namespace MantleConvection
