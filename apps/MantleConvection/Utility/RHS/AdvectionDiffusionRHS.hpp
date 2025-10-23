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
#include "hyteg/memory/TempFunctionManager.hpp"
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

template < class MassOperatorType_             = hyteg::NoOperator,
           class InternalHeatingOperatorType_  = hyteg::NoOperator,
           class ShearHeatingOperatorType_     = hyteg::NoOperator,
           class AdiabaticHeatingOperatorType_ = hyteg::NoOperator,

           class MassStabilisationOperatorType_             = hyteg::NoOperator,
           class InternalHeatingStabilisationOperatorType_  = hyteg::NoOperator,
           class ShearHeatingStabilisationOperatorType_     = hyteg::NoOperator,
           class AdiabaticHeatingStabilisationOperatorType_ = hyteg::NoOperator,

           class TemperatureFunctionType_ = hyteg::P2Function< real_t >,
           class AdditionalDataType_      = real_t >
class AdvectionDiffusionOperatorRHS
{
 public:
   // clang-format off
   typedef AdditiveOperator< MassOperatorType_, MassStabilisationOperatorType_, TemperatureFunctionType_, TemperatureFunctionType_ > AdditiveMassType;

   typedef hyteg::ScaledOperator< InternalHeatingOperatorType_  , TemperatureFunctionType_ , TemperatureFunctionType_ >  InternalHeatingOperatorType;
   typedef hyteg::ScaledOperator< ShearHeatingOperatorType_     , TemperatureFunctionType_ , TemperatureFunctionType_ >  ShearHeatingOperatorType;
   typedef hyteg::ScaledOperator< AdiabaticHeatingOperatorType_ , TemperatureFunctionType_ , TemperatureFunctionType_ >  AdiabaticHeatingOperatorType;

   typedef hyteg::ScaledOperator< InternalHeatingStabilisationOperatorType_  , TemperatureFunctionType_ , TemperatureFunctionType_ >  InternalHeatingStabilisationOperatorType;
   typedef hyteg::ScaledOperator< ShearHeatingStabilisationOperatorType_     , TemperatureFunctionType_ , TemperatureFunctionType_ >  ShearHeatingStabilisationOperatorType;
   typedef hyteg::ScaledOperator< AdiabaticHeatingStabilisationOperatorType_ , TemperatureFunctionType_ , TemperatureFunctionType_ >  AdiabaticHeatingStabilisationOperatorType;

   typedef InternalHeatingOperatorType_               InternalHeatingOperatorTypeInternal;
   typedef ShearHeatingOperatorType_                  ShearHeatingOperatorTypeInternal;
   typedef AdiabaticHeatingOperatorType_              AdiabaticHeatingOperatorTypeInternal;
   typedef InternalHeatingStabilisationOperatorType_  InternalHeatingStabilisationOperatorTypeInternal;
   typedef ShearHeatingStabilisationOperatorType_     ShearHeatingStabilisationOperatorTypeInternal;
   typedef AdiabaticHeatingStabilisationOperatorType_ AdiabaticHeatingStabilisationOperatorTypeInternal;

   typedef TemperatureFunctionType_                   TemperatureFunctionType;
   typedef AdditionalDataType_                        AdditionalDataType;
   // clang-format on

   AdvectionDiffusionOperatorRHS(
       const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
       uint_t                                            minLevel,
       uint_t                                            maxLevel,

       const std::shared_ptr< hyteg::FunctionHistory< TemperatureFunctionType_, AdditionalDataType_ > >& functionHistory =
           nullptr,

       const std::shared_ptr< MassOperatorType_ >&             massOperator             = nullptr,
       const std::shared_ptr< InternalHeatingOperatorType_ >&  internalHeatingOperator  = nullptr,
       const std::shared_ptr< ShearHeatingOperatorType_ >&     shearHeatingOperator     = nullptr,
       const std::shared_ptr< AdiabaticHeatingOperatorType_ >& adiabaticHeatingOperator = nullptr,

       real_t internalHeatingOperatorScaling  = real_c( 1.0 ),
       real_t shearHeatingOperatorScaling     = real_c( 1.0 ),
       real_t adiabaticHeatingOperatorScaling = real_c( 1.0 ),

       const std::shared_ptr< hyteg::TimeDiscretisationScheme< TemperatureFunctionType_,
                                                               AdditiveMassType,
                                                               AdditionalDataType_,
                                                               TemperatureFunctionType_ > >& timeDiscretisation = nullptr,

       const std::shared_ptr< MassStabilisationOperatorType_ >&             massStabOperator             = nullptr,
       const std::shared_ptr< InternalHeatingStabilisationOperatorType_ >&  internalHeatingStabOperator  = nullptr,
       const std::shared_ptr< ShearHeatingStabilisationOperatorType_ >&     shearHeatingStabOperator     = nullptr,
       const std::shared_ptr< AdiabaticHeatingStabilisationOperatorType_ >& adiabaticHeatingStabOperator = nullptr,

       real_t internalHeatingStabOperatorScaling  = real_c( 1.0 ),
       real_t shearHeatingStabOperatorScaling     = real_c( 1.0 ),
       real_t adiabaticHeatingStabOperatorScaling = real_c( 1.0 ) )
   : storage_( storage )
   , massOperator_( massOperator )
   , massStabOperator_( massStabOperator )

   , additiveMass_( storage, minLevel, maxLevel, massOperator, massStabOperator )
   , timeDiscretisation_( timeDiscretisation )
   , functionHistory_( functionHistory )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      // Internal Heating
      if constexpr ( !std::is_same< InternalHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         if ( internalHeatingOperator == nullptr )
         {
            WALBERLA_ABORT( "Internal heating operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Shear Heating
      if constexpr ( !std::is_same< ShearHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         if ( shearHeatingOperator == nullptr )
         {
            WALBERLA_ABORT( "Shear heating operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Adiabatic Heating
      if constexpr ( !std::is_same< AdiabaticHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         if ( adiabaticHeatingOperator == nullptr )
         {
            WALBERLA_ABORT( "Adiabatic heating operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Internal Heating Stabilisation
      if constexpr ( !std::is_same< InternalHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         if ( internalHeatingStabOperator == nullptr )
         {
            WALBERLA_ABORT( "Internal heating stabilisation operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Shear Heating Stabilisation
      if constexpr ( !std::is_same< ShearHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         if ( shearHeatingStabOperator == nullptr )
         {
            WALBERLA_ABORT( "Shear heating stabilisation operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Adiabatic Heating Stabilisation
      if constexpr ( !std::is_same< AdiabaticHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         if ( adiabaticHeatingStabOperator == nullptr )
         {
            WALBERLA_ABORT( "Adiabatic heating stabilisation operator set to nullptr but type is not NoOperator!" );
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

      internalHeatingOperator_ = std::make_shared< InternalHeatingOperatorType >(
          storage, minLevel, maxLevel, internalHeatingOperator, internalHeatingOperatorScaling );
      shearHeatingOperator_ = std::make_shared< ShearHeatingOperatorType >(
          storage, minLevel, maxLevel, shearHeatingOperator, shearHeatingOperatorScaling );
      adiabaticHeatingOperator_ = std::make_shared< AdiabaticHeatingOperatorType >(
          storage, minLevel, maxLevel, adiabaticHeatingOperator, adiabaticHeatingOperatorScaling );

      internalHeatingStabOperator_ = std::make_shared< InternalHeatingStabilisationOperatorType >(
          storage, minLevel, maxLevel, internalHeatingStabOperator, internalHeatingStabOperatorScaling );
      shearHeatingStabOperator_ = std::make_shared< ShearHeatingStabilisationOperatorType >(
          storage, minLevel, maxLevel, shearHeatingStabOperator, shearHeatingStabOperatorScaling );
      adiabaticHeatingStabOperator_ = std::make_shared< AdiabaticHeatingStabilisationOperatorType >(
          storage, minLevel, maxLevel, adiabaticHeatingStabOperator, adiabaticHeatingStabOperatorScaling );
   }

   // temperature is usually an extrapolation of the full temperature ( without subtracting a profile )
   void applyTimeIndependent( const TemperatureFunctionType_& temperature,
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
      // ####### Linear Form Workaround #######
      // ######################################

      // TODO: Proper linear forms in the HOG
      std::shared_ptr< TemperatureFunctionType_ > O =
          hyteg::getTemporaryFunction< TemperatureFunctionType_ >( storage_, level, level );
      O->interpolate( real_c( 1 ), level, hyteg::All );

      // ######################################
      // ################ Main ################
      // ######################################

      // Internal Heating ( linear form )
      if constexpr ( !std::is_same< InternalHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         internalHeatingOperator_->apply( *O, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Shear Heating ( linear form, the velocity and viscosity go into the shear heating via an external FEM function )
      if constexpr ( !std::is_same< ShearHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         shearHeatingOperator_->apply( *O, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Adiabatic Heating ( temperature dependent, the velocity goes into the adiabatic heating via an external FEM function )
      if constexpr ( !std::is_same< AdiabaticHeatingOperatorType_, hyteg::NoOperator >::value )
      {
         adiabaticHeatingOperator_->apply(
             temperature, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // #######################################
      // ############ Stabilisation ############
      // #######################################

      // Internal Heating Stabilisation ( linear form, the velocity goes into the internal heating stabilisation via an external FEM function  )
      if constexpr ( !std::is_same< InternalHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         internalHeatingStabOperator_->apply(
             *O, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Shear Heating Stabilisation ( linear form, the velocity and viscosity go into the shear heating stabilisation via an external FEM function )
      if constexpr ( !std::is_same< ShearHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         shearHeatingStabOperator_->apply( *O, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }

      // Adiabatic Heating Stabilisation ( temperature dependent, the velocity goes into the adiabatic heating stabilisation via an external FEM function )
      if constexpr ( !std::is_same< AdiabaticHeatingStabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         adiabaticHeatingStabOperator_->apply(
             temperature, dst, level, flag, ( firstOperatorApply ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApply = false;
      }
   }

   // temperature is usually an extrapolation of the full temperature ( without subtracting a profile )
   void apply( const TemperatureFunctionType_& temperature,
               const TemperatureFunctionType_& dst,
               const uint_t                    level,
               const hyteg::DoFType            flag,
               const hyteg::UpdateType         updateType = hyteg::UpdateType::Replace ) const
   {
      applyTimeIndependent( temperature, dst, level, flag, updateType );

      // #######################################
      // ######### Time Discretisation #########
      // #######################################

      if constexpr ( ( !std::is_same< MassOperatorType_, hyteg::NoOperator >::value ) ||
                     ( !std::is_same< MassStabilisationOperatorType_, hyteg::NoOperator >::value ) )
      {
         timeDiscretisation_->applyRHS( *functionHistory_, dst, additiveMass_, level, flag );
      }
   }

   // direct get functions for the internal operators
   AdditiveMassType& getAdditiveMass() { return *additiveMass_; }

   InternalHeatingOperatorType&  getInternalHeating() { return *internalHeatingOperator_; }
   ShearHeatingOperatorType&     getShearHeating() { return *shearHeatingOperator_; }
   AdiabaticHeatingOperatorType& getAdiabaticHeating() { return *adiabaticHeatingOperator_; }

   InternalHeatingStabilisationOperatorType&  getInternalHeatingStabilisation() { return *internalHeatingStabOperator_; }
   ShearHeatingStabilisationOperatorType&     getShearHeatingStabilisation() { return *shearHeatingStabOperator_; }
   AdiabaticHeatingStabilisationOperatorType& getAdiabaticHeatingStabilisation() { return *adiabaticHeatingStabOperator_; }

   hyteg::TimeDiscretisationScheme< TemperatureFunctionType_, AdditiveMassType, AdditionalDataType_, TemperatureFunctionType_ >&
       getTimeDiscretisation()
   {
      return *timeDiscretisation_;
   }
   hyteg::FunctionHistory< TemperatureFunctionType_, AdditionalDataType_ >& getFunctionHistory() { return *functionHistory_; }

   // pointer get functions for the internal operators
   std::shared_ptr< AdditiveMassType > getAdditiveMassPtr() { return additiveMass_; }

   std::shared_ptr< InternalHeatingOperatorType >  getInternalHeatingPtr() { return internalHeatingOperator_; }
   std::shared_ptr< ShearHeatingOperatorType >     getShearHeatingPtr() { return shearHeatingOperator_; }
   std::shared_ptr< AdiabaticHeatingOperatorType > getAdiabaticHeatingPtr() { return adiabaticHeatingOperator_; }

   std::shared_ptr< InternalHeatingStabilisationOperatorType > getInternalHeatingStabilisationPtr()
   {
      return internalHeatingStabOperator_;
   }
   std::shared_ptr< ShearHeatingStabilisationOperatorType > getShearHeatingStabilisationPtr()
   {
      return shearHeatingStabOperator_;
   }
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
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;

   std::shared_ptr< MassOperatorType_ >              massOperator_;
   std::shared_ptr< MassStabilisationOperatorType_ > massStabOperator_;

   std::shared_ptr< InternalHeatingOperatorType >  internalHeatingOperator_;
   std::shared_ptr< ShearHeatingOperatorType >     shearHeatingOperator_;
   std::shared_ptr< AdiabaticHeatingOperatorType > adiabaticHeatingOperator_;

   std::shared_ptr< InternalHeatingStabilisationOperatorType >  internalHeatingStabOperator_;
   std::shared_ptr< ShearHeatingStabilisationOperatorType >     shearHeatingStabOperator_;
   std::shared_ptr< AdiabaticHeatingStabilisationOperatorType > adiabaticHeatingStabOperator_;

   AdditiveMassType additiveMass_;

   std::shared_ptr<
       hyteg::
           TimeDiscretisationScheme< TemperatureFunctionType_, AdditiveMassType, AdditionalDataType_, TemperatureFunctionType_ > >
                                                                                              timeDiscretisation_;
   std::shared_ptr< hyteg::FunctionHistory< TemperatureFunctionType_, AdditionalDataType_ > > functionHistory_;
};

} // namespace MantleConvection
