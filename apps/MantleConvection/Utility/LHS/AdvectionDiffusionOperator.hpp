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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/functions/FunctionHistory.hpp"
#include "hyteg/functions/PressureMeanProjection.hpp"
#include "hyteg/numerictools/TimeDiscretisationScheme.hpp"
#include "hyteg/operators/GEMV.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/ScaledOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

#include "../OperatorTools/AdditiveOperator.hpp"

using hyteg::idx_t;
using walberla::real_t;
using walberla::uint_t;

namespace MantleConvection {

template < class MassOperatorType_                = hyteg::NoOperator,
           class AdvectionOperatorType_           = hyteg::NoOperator,
           class DiffusionOperatorType_           = hyteg::NoOperator,
           class DiffusionAdditionalOperatorType_ = hyteg::NoOperator,
           class AdiabaticHeatingOperatorType_    = hyteg::NoOperator,

           class MassStabilisationOperatorType_             = hyteg::NoOperator,
           class AdvectionStabilisationOperatorType_        = hyteg::NoOperator,
           class DiffusionStabilisationOperatorType_        = hyteg::NoOperator,
           class AdiabaticHeatingStabilisationOperatorType_ = hyteg::NoOperator,

           class TemperatureFunctionType_           = hyteg::P2Function< real_t >,
           class AdvectionDiffusionIdxFunctionType_ = hyteg::P2Function< idx_t >,
           class AdditionalDataType_                = real_t >
class AdvectionDiffusionOperator : public hyteg::Operator< TemperatureFunctionType_, TemperatureFunctionType_ >
{
 public:
   // clang-format off
   typedef AdditiveOperator< MassOperatorType_, MassStabilisationOperatorType_, TemperatureFunctionType_, TemperatureFunctionType_ > AdditiveMassType;

   typedef TemperatureFunctionType_           TemperatureFunctionType;
   typedef AdvectionDiffusionIdxFunctionType_ AdvectionDiffusionIdxFunctionType;
   typedef AdditionalDataType_                AdditionalDataType;

   typedef hyteg::ScaledOperator< AdvectionOperatorType_          , TemperatureFunctionType_, TemperatureFunctionType_ >  AdvectionOperatorType;
   typedef hyteg::ScaledOperator< DiffusionOperatorType_          , TemperatureFunctionType_, TemperatureFunctionType_ >  DiffusionOperatorType;
   typedef hyteg::ScaledOperator< DiffusionAdditionalOperatorType_, TemperatureFunctionType_, TemperatureFunctionType_ >  DiffusionAdditionalOperatorType;
   typedef hyteg::ScaledOperator< AdiabaticHeatingOperatorType_   , TemperatureFunctionType_, TemperatureFunctionType_ >  AdiabaticHeatingOperatorType;

   typedef hyteg::ScaledOperator< AdvectionStabilisationOperatorType_       , TemperatureFunctionType_, TemperatureFunctionType_ >  AdvectionStabilisationOperatorType;
   typedef hyteg::ScaledOperator< DiffusionStabilisationOperatorType_       , TemperatureFunctionType_, TemperatureFunctionType_ >  DiffusionStabilisationOperatorType;
   typedef hyteg::ScaledOperator< AdiabaticHeatingStabilisationOperatorType_, TemperatureFunctionType_, TemperatureFunctionType_ >  AdiabaticHeatingStabilisationOperatorType;

   typedef MassOperatorType_                             MassOperatorTypeInternal;
   typedef MassStabilisationOperatorType_                MassStabilisationOperatorTypeInternal;

   typedef AdvectionOperatorType_                        AdvectionOperatorTypeInternal;
   typedef DiffusionOperatorType_                        DiffusionOperatorTypeInternal;
   typedef DiffusionAdditionalOperatorType_              DiffusionAdditionalOperatorTypeInternal;
   typedef AdiabaticHeatingOperatorType_                 AdiabaticHeatingOperatorTypeInternal;
   typedef AdvectionStabilisationOperatorType_           AdvectionStabilisationOperatorTypeInternal;
   typedef DiffusionStabilisationOperatorType_           DiffusionStabilisationOperatorTypeInternal;
   typedef AdiabaticHeatingStabilisationOperatorType_    AdiabaticHeatingStabilisationOperatorTypeInternal;
   // clang-format on

   AdvectionDiffusionOperator(
       const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
       uint_t                                            minLevel,
       uint_t                                            maxLevel,

       const std::shared_ptr< hyteg::FunctionHistory< TemperatureFunctionType_, AdditionalDataType_ > >& functionHistory =
           nullptr,

       const std::shared_ptr< MassOperatorType_ >&                massOperator                = nullptr,
       const std::shared_ptr< AdvectionOperatorType_ >&           advectionOperator           = nullptr,
       const std::shared_ptr< DiffusionOperatorType_ >&           diffusionOperator           = nullptr,
       const std::shared_ptr< DiffusionAdditionalOperatorType_ >& diffusionAdditionalOperator = nullptr,
       const std::shared_ptr< AdiabaticHeatingOperatorType_ >&    adiabaticHeatingOperator    = nullptr,

       real_t advectionOperatorScaling           = real_c( 1.0 ),
       real_t diffusionOperatorScaling           = real_c( 1.0 ),
       real_t diffusionAdditionalOperatorScaling = real_c( 1.0 ),
       real_t adiabaticHeatingOperatorScaling    = real_c( 1.0 ),

       const std::shared_ptr< hyteg::TimeDiscretisationScheme< TemperatureFunctionType_,
                                                               AdditiveMassType,
                                                               AdditionalDataType_,
                                                               TemperatureFunctionType_ > >& timeDiscretisation = nullptr,

       const std::shared_ptr< MassStabilisationOperatorType_ >&             massStabOperator             = nullptr,
       const std::shared_ptr< AdvectionStabilisationOperatorType_ >&        advectionStabOperator        = nullptr,
       const std::shared_ptr< DiffusionStabilisationOperatorType_ >&        diffusionStabOperator        = nullptr,
       const std::shared_ptr< AdiabaticHeatingStabilisationOperatorType_ >& adiabaticHeatingStabOperator = nullptr,

       real_t advectionStabOperatorScaling        = real_c( 1.0 ),
       real_t diffusionStabOperatorScaling        = real_c( 1.0 ),
       real_t adiabaticHeatingStabOperatorScaling = real_c( 1.0 ) )
   : hyteg::Operator< TemperatureFunctionType_, TemperatureFunctionType_ >( storage, minLevel, maxLevel )
   , massOperator_( massOperator )
   , massStabOperator_( massStabOperator )
   , timeDiscretisation_( timeDiscretisation )
   , functionHistory_( functionHistory )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      // Advection
      if constexpr ( !std::is_same< AdvectionOperatorType_, hyteg::NoOperator >::value )
      {
         if ( advectionOperator == nullptr )
         {
            WALBERLA_ABORT( "Advection operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Diffusion
      if constexpr ( !std::is_same< DiffusionOperatorType_, hyteg::NoOperator >::value )
      {
         if ( diffusionOperator == nullptr )
         {
            WALBERLA_ABORT( "Diffusion operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Diffusion Additional
      if constexpr ( !std::is_same< DiffusionAdditionalOperatorType_, hyteg::NoOperator >::value )
      {
         if ( diffusionAdditionalOperator == nullptr )
         {
            WALBERLA_ABORT( "DiffusionAdditional operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Adiabatic Heating
      if constexpr ( !std::is_same< AdiabaticHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         if ( adiabaticHeatingOperator == nullptr )
         {
            WALBERLA_ABORT( "Adiabatic Heating operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Advection Stabilisation
      if constexpr ( !std::is_same< AdvectionStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         if ( advectionStabOperator == nullptr )
         {
            WALBERLA_ABORT( "Advection Stabilisation operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Diffusion Stabilisation
      if constexpr ( !std::is_same< DiffusionStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         if ( diffusionStabOperator == nullptr )
         {
            WALBERLA_ABORT( "Diffusion Stabilisation operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Adiabatic Heating Stabilisation
      if constexpr ( !std::is_same< AdiabaticHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         if ( adiabaticHeatingStabOperator == nullptr )
         {
            WALBERLA_ABORT( "Adiabatic Heating Stabilisation operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Time Discretisation
      if constexpr ( ( !std::is_same< MassOperatorType_, hyteg::NoOperator >::value ) ||
                     ( !std::is_same< MassStabilisationOperatorType_, hyteg::NoOperator >::value ) )
      {
         if ( timeDiscretisation == nullptr )
         {
            WALBERLA_ABORT( "Time Discretisation set to nullptr but one of the mass types is not NoOperator!" );
         }
      }

      // ##################################
      // #### Define wrapped operators ####
      // ##################################

      additiveMass_       = std::make_shared< AdditiveMassType >( storage, minLevel, maxLevel, massOperator, massStabOperator );
      additiveMassNoStab_ = std::make_shared< AdditiveMassType >( storage, minLevel, maxLevel, massOperator, nullptr );

      advectionOperator_ =
          std::make_shared< AdvectionOperatorType >( storage, minLevel, maxLevel, advectionOperator, advectionOperatorScaling );
      diffusionOperator_ =
          std::make_shared< DiffusionOperatorType >( storage, minLevel, maxLevel, diffusionOperator, diffusionOperatorScaling );
      diffusionAdditionalOperator_ = std::make_shared< DiffusionAdditionalOperatorType >(
          storage, minLevel, maxLevel, diffusionAdditionalOperator, diffusionAdditionalOperatorScaling );
      adiabaticHeatingOperator_ = std::make_shared< AdiabaticHeatingOperatorType >(
          storage, minLevel, maxLevel, adiabaticHeatingOperator, adiabaticHeatingOperatorScaling );

      advectionStabOperator_ = std::make_shared< AdvectionStabilisationOperatorType >(
          storage, minLevel, maxLevel, advectionStabOperator, advectionStabOperatorScaling );
      diffusionStabOperator_ = std::make_shared< DiffusionStabilisationOperatorType >(
          storage, minLevel, maxLevel, diffusionStabOperator, diffusionStabOperatorScaling );
      adiabaticHeatingStabOperator_ = std::make_shared< AdiabaticHeatingStabilisationOperatorType >(
          storage, minLevel, maxLevel, adiabaticHeatingStabOperator, adiabaticHeatingStabOperatorScaling );
   }

   void applyTimeIndependent( const TemperatureFunctionType_& src,
                              const TemperatureFunctionType_& dst,
                              const uint_t                    level,
                              const hyteg::DoFType            flag,
                              const hyteg::UpdateType         updateType = hyteg::UpdateType::Replace ) const
   {
      // #######################################
      // ###### Init First Operator Apply ######
      // #######################################

      bool firstOperatorApply = true;

      // ######################################
      // ################ Main ################
      // ######################################

      // Advection
      if constexpr ( !std::is_same< AdvectionOperatorType_, hyteg::NoOperator >::value )
      {
         advectionOperator_->apply( src, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Diffusion
      if constexpr ( !std::is_same< DiffusionOperatorType_, hyteg::NoOperator >::value )
      {
         diffusionOperator_->apply( src, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Diffusion Additional
      if constexpr ( !std::is_same< DiffusionAdditionalOperatorType_, hyteg::NoOperator >::value )
      {
         diffusionAdditionalOperator_->apply(
             src, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Adiabatic Heating
      if constexpr ( !std::is_same< AdiabaticHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         adiabaticHeatingOperator_->apply( src, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // #######################################
      // ############ Stabilisation ############
      // #######################################

      // Advection Stabilisation
      if constexpr ( !std::is_same< AdvectionStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         advectionStabOperator_->apply( src, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Diffusion Stabilisation
      if constexpr ( !std::is_same< DiffusionStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         diffusionStabOperator_->apply( src, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Adiabatic Heating Stabilisation
      if constexpr ( !std::is_same< AdiabaticHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         adiabaticHeatingStabOperator_->apply(
             src, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }
   }

   void apply( const TemperatureFunctionType_& src,
               const TemperatureFunctionType_& dst,
               const uint_t                    level,
               const hyteg::DoFType            flag,
               const hyteg::UpdateType         updateType = hyteg::UpdateType::Replace ) const
   {
      applyTimeIndependent( src, dst, level, flag, updateType );

      // #######################################
      // ######### Time Discretisation #########
      // #######################################

      if constexpr ( ( !std::is_same< MassOperatorType_, hyteg::NoOperator >::value ) ||
                     ( !std::is_same< MassStabilisationOperatorType_, hyteg::NoOperator >::value ) )
      {
         timeDiscretisation_->applyLHS( *functionHistory_, src, dst, *additiveMass_, level, flag );
      }
   }

   void toMatrixTimeIndependent( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                                 const AdvectionDiffusionIdxFunctionType_&          src,
                                 const AdvectionDiffusionIdxFunctionType_&          dst,
                                 uint_t                                             level,
                                 hyteg::DoFType                                     flag ) const
   {
      // ######################################
      // ################ Main ################
      // ######################################

      // Advection
      if constexpr ( !std::is_same< AdvectionOperatorType_, hyteg::NoOperator >::value )
      {
         advectionOperator_->toMatrix( mat, src, dst, level, flag );
      }

      // Diffusion
      if constexpr ( !std::is_same< DiffusionOperatorType_, hyteg::NoOperator >::value )
      {
         diffusionOperator_->toMatrix( mat, src, dst, level, flag );
      }

      // Diffusion Additional
      if constexpr ( !std::is_same< DiffusionAdditionalOperatorType_, hyteg::NoOperator >::value )
      {
         diffusionAdditionalOperator_->toMatrix( mat, src, dst, level, flag );
      }

      // Adiabatic Heating
      if constexpr ( !std::is_same< AdiabaticHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         adiabaticHeatingOperator_->toMatrix( mat, src, dst, level, flag );
      }

      // #######################################
      // ############ Stabilisation ############
      // #######################################

      // Advection Stabilisation
      if constexpr ( !std::is_same< AdvectionStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         advectionStabOperator_->toMatrix( mat, src, dst, level, flag );
      }

      // Diffusion Stabilisation
      if constexpr ( !std::is_same< DiffusionStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         diffusionStabOperator_->toMatrix( mat, src, dst, level, flag );
      }

      // Adiabatic Heating Stabilisation
      if constexpr ( !std::is_same< AdiabaticHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         adiabaticHeatingStabOperator_->toMatrix( mat, src, dst, level, flag );
      }
   }

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                  const AdvectionDiffusionIdxFunctionType_&          src,
                  const AdvectionDiffusionIdxFunctionType_&          dst,
                  uint_t                                             level,
                  hyteg::DoFType                                     flag ) const
   {
      toMatrixTimeIndependent( mat, src, dst, level, flag );

      // #######################################
      // ######### Time Discretisation #########
      // #######################################

      if constexpr ( ( !std::is_same< MassOperatorType_, hyteg::NoOperator >::value ) ||
                     ( !std::is_same< MassStabilisationOperatorType_, hyteg::NoOperator >::value ) )
      {
         timeDiscretisation_->toMatrix( *functionHistory_, *additiveMass_, mat, src, dst, level, flag );
      }
   }

   // direct get functions for the internal operators
   AdditiveMassType& getAdditiveMass() { return *additiveMass_; }
   AdditiveMassType& getAdditiveMassNoStab() { return *additiveMassNoStab_; }

   AdvectionOperatorType&           getAdvection() { return *advectionOperator_; }
   DiffusionOperatorType&           getDiffusion() { return *diffusionOperator_; }
   DiffusionAdditionalOperatorType& getDiffusionAdditional() { return *diffusionAdditionalOperator_; }
   AdiabaticHeatingOperatorType&    getAdiabaticHeating() { return *adiabaticHeatingOperator_; }

   AdvectionStabilisationOperatorType&        getAdvectionStabilisation() { return *advectionStabOperator_; }
   DiffusionStabilisationOperatorType&        getDiffusionStabilisation() { return *diffusionStabOperator_; }
   AdiabaticHeatingStabilisationOperatorType& getAdiabaticHeatingStabilisation() { return *adiabaticHeatingStabOperator_; }

   hyteg::TimeDiscretisationScheme< TemperatureFunctionType_, AdditiveMassType, AdditionalDataType_, TemperatureFunctionType_ >&
       getTimeDiscretisation()
   {
      return *timeDiscretisation_;
   }
   hyteg::FunctionHistory< TemperatureFunctionType_, AdditionalDataType_ >& getFunctionHistory() { return *functionHistory_; }

   // pointer get functions for the internal operators
   std::shared_ptr< AdditiveMassType > getAdditiveMassPtr() { return additiveMass_; }
   std::shared_ptr< AdditiveMassType > getAdditiveMassNoStabPtr() { return additiveMassNoStab_; }

   std::shared_ptr< AdvectionOperatorType >           getAdvectionPtr() { return advectionOperator_; }
   std::shared_ptr< DiffusionOperatorType >           getDiffusionPtr() { return diffusionOperator_; }
   std::shared_ptr< DiffusionAdditionalOperatorType > getDiffusionAdditionalPtr() { return diffusionAdditionalOperator_; }
   std::shared_ptr< AdiabaticHeatingOperatorType >    getAdiabaticHeatingPtr() { return adiabaticHeatingOperator_; }

   std::shared_ptr< AdvectionStabilisationOperatorType >        getAdvectionStabilisationPtr() { return advectionStabOperator_; }
   std::shared_ptr< DiffusionStabilisationOperatorType >        getDiffusionStabilisationPtr() { return diffusionStabOperator_; }
   std::shared_ptr< AdiabaticHeatingStabilisationOperatorType > getAdiabaticHeatingStabilisationPtr()
   {
      return adiabaticHeatingStabOperator_;
   }

   std::shared_ptr<
       hyteg::
           TimeDiscretisationScheme< TemperatureFunctionType_, AdditiveMassType, AdditionalDataType_, TemperatureFunctionType_ > >
       getTimeDiscretisationPtr()
   {
      return timeDiscretisation_;
   }
   std::shared_ptr< hyteg::FunctionHistory< TemperatureFunctionType_, AdditionalDataType_ > > getFunctionHistoryPtr()
   {
      return functionHistory_;
   }

 private:
   std::shared_ptr< MassOperatorType_ >              massOperator_;
   std::shared_ptr< MassStabilisationOperatorType_ > massStabOperator_;

   std::shared_ptr< AdvectionOperatorType >           advectionOperator_;
   std::shared_ptr< DiffusionOperatorType >           diffusionOperator_;
   std::shared_ptr< DiffusionAdditionalOperatorType > diffusionAdditionalOperator_;
   std::shared_ptr< AdiabaticHeatingOperatorType >    adiabaticHeatingOperator_;

   std::shared_ptr< AdvectionStabilisationOperatorType >        advectionStabOperator_;
   std::shared_ptr< DiffusionStabilisationOperatorType >        diffusionStabOperator_;
   std::shared_ptr< AdiabaticHeatingStabilisationOperatorType > adiabaticHeatingStabOperator_;

   std::shared_ptr< AdditiveMassType > additiveMass_;
   std::shared_ptr< AdditiveMassType > additiveMassNoStab_;

   std::shared_ptr<
       hyteg::
           TimeDiscretisationScheme< TemperatureFunctionType_, AdditiveMassType, AdditionalDataType_, TemperatureFunctionType_ > >
                                                                                              timeDiscretisation_;
   std::shared_ptr< hyteg::FunctionHistory< TemperatureFunctionType_, AdditionalDataType_ > > functionHistory_;
};

} // namespace MantleConvection
