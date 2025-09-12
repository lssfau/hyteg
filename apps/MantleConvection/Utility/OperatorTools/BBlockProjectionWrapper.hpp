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

#include "hyteg/functions/PressureMeanProjection.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

namespace MantleConvection {

template < typename OperatorType,
           typename VelocityProjectionOperatorType = hyteg::NoOperator,
           bool preProjection                      = true,
           bool postProjection                     = true,
           bool allowPreProjectionToChangeSrc      = true,
           bool allowPostProjectionToChangeDst     = true,
           bool includeProjectionsToMatrix         = false >
class BBlockProjectionWrapper : public hyteg::Operator< typename OperatorType::srcType, typename OperatorType::dstType >
{
 public:
   BBlockProjectionWrapper( std::shared_ptr< OperatorType >                   op,
                            std::shared_ptr< VelocityProjectionOperatorType > projection     = nullptr,
                            hyteg::DoFType                                    projectionFlag = hyteg::FreeslipBoundary,
                            bool                                              lowMemoryMode  = false )
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

      if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
      {
         if ( projection == nullptr )
         {
            WALBERLA_ABORT( "Pre projection activated and not of type NoOperator but projection is nullptr!" );
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
                "BBlockProjectionWrapper tmpSrc", op->getStorage(), op->getMinLevel(), op->getMaxLevel() );
         }
      }
      if constexpr ( useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmpDst_ = std::make_shared< typename OperatorType::dstType >(
                "BBlockProjectionWrapper tmpDst", op->getStorage(), op->getMinLevel(), op->getMaxLevel() );
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
      if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
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
      if constexpr ( postProjection )
      {
         if constexpr ( allowPostProjectionToChangeDst )
         {
            hyteg::projectPressureMean( dst, level );
         }
         else
         {
            hyteg::projectPressureMean( *tmpDstApply, level );
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
         if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
         {
            // B Operator goes from Velocity Vector Function to Pressure Function
            // Hence we need a matrix of correct size for a Velocity Vector Function to Velocity Vector Function operator
            const uint_t localRows  = numberOfLocalDoFs( numeratorSrc, level );
            const uint_t localCols  = numberOfLocalDoFs( numeratorSrc, level );
            const uint_t globalRows = numberOfGlobalDoFs( numeratorSrc, level, walberla::mpi::MPIManager::instance()->comm() );
            const uint_t globalCols = numberOfGlobalDoFs( numeratorSrc, level, walberla::mpi::MPIManager::instance()->comm() );

            matProxyProjectionPre =
                mat->createMatrix( localRows, localCols, globalRows, globalCols, walberla::mpi::MPIManager::instance()->comm() );
            projection_->toMatrix( matProxyProjectionPre, numeratorSrc, numeratorSrc, level, projectionFlag_ );
         }

         std::vector< std::shared_ptr< hyteg::SparseMatrixProxy > > matrices;

         // TODO: project mean matrix?

         matrices.push_back( matProxyOp );

         if constexpr ( ( preProjection ) && ( !std::is_same< VelocityProjectionOperatorType, hyteg::NoOperator >::value ) )
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

   std::shared_ptr< OperatorType > getOperatorPtr() { return op_; }

   static constexpr bool useTmpSrc = ( !allowPreProjectionToChangeSrc ) && ( preProjection ) &&
                                     ( !std::is_same< VelocityProjectionOperatorType, hyteg::NoOperator >::value );
   static constexpr bool useTmpDst = ( !allowPostProjectionToChangeDst ) && ( postProjection );

 private:
   std::shared_ptr< OperatorType >                   op_;
   std::shared_ptr< VelocityProjectionOperatorType > projection_;

   hyteg::DoFType projectionFlag_;

   std::shared_ptr< typename OperatorType::srcType > tmpSrc_;
   std::shared_ptr< typename OperatorType::dstType > tmpDst_;

   bool lowMemoryMode_;
};

} // namespace MantleConvection
