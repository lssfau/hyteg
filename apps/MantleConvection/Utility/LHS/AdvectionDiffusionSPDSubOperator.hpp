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

#include "AdvectionDiffusionOperator.hpp"

using hyteg::idx_t;
using walberla::real_t;
using walberla::uint_t;

namespace MantleConvection {

template < class OperatorType >
class AdvectionDiffusionSPDSubOperator
: public hyteg::Operator< typename OperatorType::TemperatureFunctionType, typename OperatorType::TemperatureFunctionType >
{
 public:
   typedef typename OperatorType::TemperatureFunctionType           TemperatureFunctionType;
   typedef typename OperatorType::AdvectionDiffusionIdxFunctionType AdvectionDiffusionIdxFunctionType;

   AdvectionDiffusionSPDSubOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                     uint_t                                            minLevel,
                                     uint_t                                            maxLevel,
                                     const std::shared_ptr< OperatorType >&            op = nullptr )
   : hyteg::Operator< typename OperatorType::TemperatureFunctionType, typename OperatorType::TemperatureFunctionType >( storage,
                                                                                                                        minLevel,
                                                                                                                        maxLevel )
   , op_( op )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      if constexpr ( !std::is_same< typename OperatorType::DiffusionOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         if ( op->getDiffusionPtr()->getOperatorPtr() == nullptr )
         {
            WALBERLA_ABORT( "Diffusion operator set to nullptr but type is not NoOperator!" );
         }
      }

      if constexpr ( ( !std::is_same< typename OperatorType::MassOperatorTypeInternal, hyteg::NoOperator >::value ) ||
                     ( !std::is_same< typename OperatorType::MassStabilisationOperatorTypeInternal, hyteg::NoOperator >::value ) )
      {
         if ( op->getTimeDiscretisationPtr() == nullptr )
         {
            WALBERLA_ABORT( "Time Discretisation set to nullptr but one of the mass types is not NoOperator!" );
         }
      }
   }

   void applyTimeIndependent( const TemperatureFunctionType& src,
                              const TemperatureFunctionType& dst,
                              const uint_t                   level,
                              const hyteg::DoFType           flag,
                              const hyteg::UpdateType        updateType = hyteg::UpdateType::Replace ) const
   {
      // ######################################
      // ################ Main ################
      // ######################################

      // Diffusion
      if constexpr ( !std::is_same< typename OperatorType::DiffusionOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         op_->getDiffusion().apply( src, dst, level, flag, updateType );
      }
   }

   void apply( const TemperatureFunctionType& src,
               const TemperatureFunctionType& dst,
               const uint_t                   level,
               const hyteg::DoFType           flag,
               const hyteg::UpdateType        updateType = hyteg::UpdateType::Replace ) const
   {
      applyTimeIndependent( src, dst, level, flag, updateType );

      // #######################################
      // ######### Time Discretisation #########
      // #######################################

      if constexpr ( ( !std::is_same< typename OperatorType::MassOperatorTypeInternal, hyteg::NoOperator >::value ) ||
                     ( !std::is_same< typename OperatorType::MassStabilisationOperatorTypeInternal, hyteg::NoOperator >::value ) )
      {
         op_->getTimeDiscretisation().applyLHS( op_->getFunctionHistory(), src, dst, op_->getAdditiveMassNoStab(), level, flag );
      }
   }

   void toMatrixTimeIndependent( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                                 const AdvectionDiffusionIdxFunctionType&           src,
                                 const AdvectionDiffusionIdxFunctionType&           dst,
                                 uint_t                                             level,
                                 hyteg::DoFType                                     flag ) const
   {
      // ######################################
      // ################ Main ################
      // ######################################

      // Diffusion
      if constexpr ( !std::is_same< typename OperatorType::DiffusionOperatorTypeInternal, hyteg::NoOperator >::value )
      {
         op_->getDiffusion().toMatrix( mat, src, dst, level, flag );
      }
   }

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                  const AdvectionDiffusionIdxFunctionType&           src,
                  const AdvectionDiffusionIdxFunctionType&           dst,
                  uint_t                                             level,
                  hyteg::DoFType                                     flag ) const
   {
      toMatrixTimeIndependent( mat, src, dst, level, flag );

      // #######################################
      // ######### Time Discretisation #########
      // #######################################

      if constexpr ( ( !std::is_same< typename OperatorType::MassOperatorTypeInternal, hyteg::NoOperator >::value ) ||
                     ( !std::is_same< typename OperatorType::MassStabilisationOperatorTypeInternal, hyteg::NoOperator >::value ) )
      {
         op_->getTimeDiscretisation().toMatrix(
             op_->getFunctionHistory(), op_->getAdditiveMassNoStab(), mat, src, dst, level, flag );
      }
   }

 private:
   std::shared_ptr< OperatorType > op_;
};

} // namespace MantleConvection
