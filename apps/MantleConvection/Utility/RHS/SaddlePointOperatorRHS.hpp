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
#include "hyteg/functions/PressureMeanProjection.hpp"
#include "hyteg/operators/GEMV.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/ScaledOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

using hyteg::idx_t;
using walberla::real_t;
using walberla::uint_t;

namespace MantleConvection {

template < class VelocityToVelocityRHSOperatorType_    = hyteg::NoOperator,
           class PressureToVelocityRHSOperatorType_    = hyteg::NoOperator,
           class TemperatureToVelocityRHSOperatorType_ = hyteg::NoOperator,

           class VelocityToPressureRHSOperatorType_    = hyteg::NoOperator,
           class PressureToPressureRHSOperatorType_    = hyteg::NoOperator,
           class TemperatureToPressureRHSOperatorType_ = hyteg::NoOperator,

           class DensityDerivativeToPressureRHSOperatorType_ = hyteg::NoOperator,

           class VelocityProjectionOperatorType_ = hyteg::NoOperator,

           bool preProjectPressure_  = true,
           bool preProjectVelocity_  = true,
           bool postProjectPressure_ = true,
           bool postProjectVelocity_ = true,

           bool allowPreProjectionToChangePressureSrc_ = true,
           bool allowPreProjectionToChangeVelocitySrc_ = true,

           class VelocityFunctionType_          = hyteg::P2VectorFunction< real_t >,
           class PressureFunctionType_          = hyteg::P1Function< real_t >,
           class TemperatureFunctionType_       = hyteg::P2Function< real_t >,
           class DensityDerivativeFunctionType_ = hyteg::P1Function< real_t >,
           class DstFunctionType_               = hyteg::P2P1TaylorHoodFunction< real_t > >
class SaddlePointOperatorRHS
{
 public:
   // clang-format off
   typedef hyteg::ScaledOperator< VelocityToVelocityRHSOperatorType_    , VelocityFunctionType_   , VelocityFunctionType_ >  VelocityToVelocityRHSOperatorType;
   typedef hyteg::ScaledOperator< PressureToVelocityRHSOperatorType_    , PressureFunctionType_   , VelocityFunctionType_ >  PressureToVelocityRHSOperatorType;
   typedef hyteg::ScaledOperator< TemperatureToVelocityRHSOperatorType_ , TemperatureFunctionType_, VelocityFunctionType_ >  TemperatureToVelocityRHSOperatorType;

   typedef hyteg::ScaledOperator< VelocityToPressureRHSOperatorType_    , VelocityFunctionType_   , PressureFunctionType_ >  VelocityToPressureRHSOperatorType;
   typedef hyteg::ScaledOperator< PressureToPressureRHSOperatorType_    , PressureFunctionType_   , PressureFunctionType_ >  PressureToPressureRHSOperatorType;
   typedef hyteg::ScaledOperator< TemperatureToPressureRHSOperatorType_ , TemperatureFunctionType_, PressureFunctionType_ >  TemperatureToPressureRHSOperatorType;

   typedef hyteg::ScaledOperator< DensityDerivativeToPressureRHSOperatorType_ , DensityDerivativeFunctionType_, PressureFunctionType_ >  DensityDerivativeToPressureRHSOperatorType;

   typedef VelocityToVelocityRHSOperatorType_            VelocityToVelocityRHSOperatorTypeInternal;
   typedef PressureToVelocityRHSOperatorType_            PressureToVelocityRHSOperatorTypeInternal;
   typedef TemperatureToVelocityRHSOperatorType_         TemperatureToVelocityRHSOperatorTypeInternal;
   
   typedef VelocityToPressureRHSOperatorType_            VelocityToPressureRHSOperatorTypeInternal;
   typedef PressureToPressureRHSOperatorType_            PressureToPressureRHSOperatorTypeInternal;
   typedef TemperatureToPressureRHSOperatorType_         TemperatureToPressureRHSOperatorTypeInternal;

   typedef DensityDerivativeToPressureRHSOperatorType_   DensityDerivativeToPressureRHSOperatorTypeInternal;

   typedef VelocityFunctionType_                         VelocityFunctionType;
   typedef PressureFunctionType_                         PressureFunctionType;
   typedef TemperatureFunctionType_                      TemperatureFunctionType;
   typedef DensityDerivativeFunctionType_                DensityDerivativeFunctionType;
   typedef DstFunctionType_                              DstFunctionType;
   // clang-format on

   SaddlePointOperatorRHS(
       const std::shared_ptr< hyteg::PrimitiveStorage >&               storage,
       uint_t                                                          minLevel,
       uint_t                                                          maxLevel,
       const std::shared_ptr< VelocityToVelocityRHSOperatorType_ >&    velocityToVelocityRHSOperator    = nullptr,
       const std::shared_ptr< PressureToVelocityRHSOperatorType_ >&    pressureToVelocityRHSOperator    = nullptr,
       const std::shared_ptr< TemperatureToVelocityRHSOperatorType_ >& temperatureToVelocityRHSOperator = nullptr,

       const std::shared_ptr< VelocityToPressureRHSOperatorType_ >&    velocityToPressureRHSOperator    = nullptr,
       const std::shared_ptr< PressureToPressureRHSOperatorType_ >&    pressureToPressureRHSOperator    = nullptr,
       const std::shared_ptr< TemperatureToPressureRHSOperatorType_ >& temperatureToPressureRHSOperator = nullptr,

       const std::shared_ptr< DensityDerivativeToPressureRHSOperatorType_ >& densityDerivativeToPressureRHSOperator = nullptr,

       real_t velocityToVelocityRHSOperatorScaling    = real_c( 1.0 ),
       real_t pressureToVelocityRHSOperatorScaling    = real_c( 1.0 ),
       real_t temperatureToVelocityRHSOperatorScaling = real_c( 1.0 ),

       real_t velocityToPressureRHSOperatorScaling    = real_c( 1.0 ),
       real_t pressureToPressureRHSOperatorScaling    = real_c( 1.0 ),
       real_t temperatureToPressureRHSOperatorScaling = real_c( 1.0 ),

       real_t densityDerivativeToPressureRHSOperatorScaling = real_c( 1.0 ),

       const std::shared_ptr< VelocityProjectionOperatorType_ >& projection     = nullptr,
       hyteg::DoFType                                            projectionFlag = hyteg::FreeslipBoundary,
       bool                                                      lowMemoryMode  = false )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , projection_( projection )
   , projectionFlag_( projectionFlag )
   , lowMemoryMode_( lowMemoryMode )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      // Velocity To Velocity

      if constexpr ( !std::is_same< VelocityToVelocityRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if ( velocityToVelocityRHSOperator == nullptr )
         {
            WALBERLA_ABORT( "Velocity to velocity operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Pressure To Velocity
      if constexpr ( !std::is_same< PressureToVelocityRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if ( pressureToVelocityRHSOperator == nullptr )
         {
            WALBERLA_ABORT( "Pressure to velocity operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Temperature To Velocity
      if constexpr ( !std::is_same< TemperatureToVelocityRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if ( temperatureToVelocityRHSOperator == nullptr )
         {
            WALBERLA_ABORT( "Temperature to velocity operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Velocity To Pressure
      if constexpr ( !std::is_same< VelocityToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if ( velocityToPressureRHSOperator == nullptr )
         {
            WALBERLA_ABORT( "Velocity to pressure operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Pressure To Pressure
      if constexpr ( !std::is_same< PressureToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if ( pressureToPressureRHSOperator == nullptr )
         {
            WALBERLA_ABORT( "Pressure to pressure operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Temperature To Pressure
      if constexpr ( !std::is_same< TemperatureToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if ( temperatureToPressureRHSOperator == nullptr )
         {
            WALBERLA_ABORT( "Temperature to pressure operator set to nullptr but type is not NoOperator!" );
         }
      }

      // Density Derivative To Pressure
      if constexpr ( !std::is_same< TemperatureToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if ( densityDerivativeToPressureRHSOperator == nullptr )
         {
            WALBERLA_ABORT( "Density to pressure operator set to nullptr but type is not NoOperator!" );
         }
      }

      if constexpr ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value )
      {
         if constexpr ( preProjectVelocity_ || postProjectVelocity_ )
         {
            if ( projection == nullptr )
            {
               WALBERLA_ABORT( "Velocity projection set to nullptr but type is not NoOperator!" );
            }
         }
      }

      // ##################################
      // #### Define wrapped operators ####
      // ##################################

      velocityToVelocityRHSOperator_ = std::make_shared< VelocityToVelocityRHSOperatorType >(
          storage, minLevel, maxLevel, velocityToVelocityRHSOperator, velocityToVelocityRHSOperatorScaling );
      pressureToVelocityRHSOperator_ = std::make_shared< PressureToVelocityRHSOperatorType >(
          storage, minLevel, maxLevel, pressureToVelocityRHSOperator, pressureToVelocityRHSOperatorScaling );
      temperatureToVelocityRHSOperator_ = std::make_shared< TemperatureToVelocityRHSOperatorType >(
          storage, minLevel, maxLevel, temperatureToVelocityRHSOperator, temperatureToVelocityRHSOperatorScaling );

      velocityToPressureRHSOperator_ = std::make_shared< VelocityToPressureRHSOperatorType >(
          storage, minLevel, maxLevel, velocityToPressureRHSOperator, velocityToPressureRHSOperatorScaling );
      pressureToPressureRHSOperator_ = std::make_shared< PressureToPressureRHSOperatorType >(
          storage, minLevel, maxLevel, pressureToPressureRHSOperator, pressureToPressureRHSOperatorScaling );
      temperatureToPressureRHSOperator_ = std::make_shared< TemperatureToPressureRHSOperatorType >(
          storage, minLevel, maxLevel, temperatureToPressureRHSOperator, temperatureToPressureRHSOperatorScaling );

      densityDerivativeToPressureRHSOperator_ = std::make_shared< DensityDerivativeToPressureRHSOperatorType >(
          storage, minLevel, maxLevel, densityDerivativeToPressureRHSOperator, densityDerivativeToPressureRHSOperatorScaling );

      // ######################################
      // #### Init tmp functions if needed ####
      // ######################################

      if constexpr ( ( preProjectPressure_ ) && ( !allowPreProjectionToChangePressureSrc_ ) )
      {
         if ( !lowMemoryMode_ )
         {
            tmpPressure_ =
                std::make_shared< PressureFunctionType_ >( "SaddlePointRHS Tmp Pressure", storage, minLevel, maxLevel );
         }
      }

      if constexpr ( ( preProjectVelocity_ ) && ( !allowPreProjectionToChangeVelocitySrc_ ) )
      {
         if ( !lowMemoryMode_ )
         {
            tmpVelocity_ =
                std::make_shared< VelocityFunctionType_ >( "SaddlePointRHS Tmp Velocity", storage, minLevel, maxLevel );
         }
      }
   }

   // velocity is usually an extrapolation
   // pressure is usually the pressure deviation from a profile
   // temperature is usually the temperature deviation from a profile
   // densityDerivative usually is a approximation of the density time derivative
   void apply( const VelocityFunctionType_&          velocity,
               const PressureFunctionType_&          pressure,
               const TemperatureFunctionType_&       temperature,
               const DensityDerivativeFunctionType_& densityDerivative,
               const DstFunctionType_&               dst,
               const uint_t                          level,
               const hyteg::DoFType                  flag,
               const hyteg::UpdateType               updateType = hyteg::UpdateType::Replace ) const
   {
      // Handle low memory mode

      std::shared_ptr< PressureFunctionType_ > pressureApply;
      std::shared_ptr< VelocityFunctionType_ > velocityApply;

      if constexpr ( preProjectPressure_ )
      {
         if constexpr ( !allowPreProjectionToChangePressureSrc_ )
         {
            if ( !lowMemoryMode_ )
            {
               pressureApply = tmpPressure_;
            }
            else
            {
               pressureApply = hyteg::getTemporaryFunction< PressureFunctionType_ >( storage_, minLevel_, maxLevel_ );
            }
         }
      }

      if constexpr ( ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) && ( preProjectVelocity_ ) )
      {
         if constexpr ( !allowPreProjectionToChangeVelocitySrc_ )
         {
            if ( !lowMemoryMode_ )
            {
               velocityApply = tmpVelocity_;
            }
            else
            {
               velocityApply = hyteg::getTemporaryFunction< VelocityFunctionType_ >( storage_, minLevel_, maxLevel_ );
            }
         }
      }

      // #######################################
      // ###### Init First Operator Apply ######
      // #######################################

      bool firstOperatorApplyVelocity = true;
      bool firstOperatorApplyPressure = true;

      // ####################################
      // ########## Pre Projection ##########
      // ####################################

      if constexpr ( preProjectPressure_ )
      {
         if constexpr ( !allowPreProjectionToChangePressureSrc_ )
         {
            pressureApply->copyBoundaryConditionFromFunction( pressure );
            pressureApply->assign( { real_c( 1 ) }, { pressure }, level, hyteg::All );
            hyteg::projectPressureMean( *pressureApply, level );
         }
         else
         {
            hyteg::projectPressureMean( pressure, level );
         }
      }

      if constexpr ( ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) && ( preProjectVelocity_ ) )
      {
         if constexpr ( !allowPreProjectionToChangeVelocitySrc_ )
         {
            velocityApply->copyBoundaryConditionFromFunction( velocity );
            velocityApply->assign( { real_c( 1 ) }, { velocity }, level, hyteg::All );
            projection_->project( *velocityApply, level, projectionFlag_ );
         }
         else
         {
            projection_->project( velocity, level, projectionFlag_ );
         }
      }

      // ######################################
      // ################ Main ################
      // ######################################

      // -------------- Velocity --------------

      // Velocity To Velocity

      if constexpr ( !std::is_same< VelocityToVelocityRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if constexpr ( ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) &&
                        ( preProjectVelocity_ ) && ( !allowPreProjectionToChangeVelocitySrc_ ) )
         {
            velocityToVelocityRHSOperator_->apply(
                *velocityApply, dst.uvw(), level, flag, ( firstOperatorApplyVelocity ? updateType : hyteg::UpdateType::Add ) );
         }
         else
         {
            velocityToVelocityRHSOperator_->apply(
                velocity, dst.uvw(), level, flag, ( firstOperatorApplyVelocity ? updateType : hyteg::UpdateType::Add ) );
         }

         firstOperatorApplyVelocity = false;
      }

      // Pressure To Velocity
      if constexpr ( !std::is_same< PressureToVelocityRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if constexpr ( ( preProjectPressure_ ) && ( !allowPreProjectionToChangePressureSrc_ ) )
         {
            pressureToVelocityRHSOperator_->apply(
                *pressureApply, dst.uvw(), level, flag, ( firstOperatorApplyVelocity ? updateType : hyteg::UpdateType::Add ) );
         }
         else
         {
            pressureToVelocityRHSOperator_->apply(
                pressure, dst.uvw(), level, flag, ( firstOperatorApplyVelocity ? updateType : hyteg::UpdateType::Add ) );
         }

         firstOperatorApplyVelocity = false;
      }

      // Temperature To Velocity
      if constexpr ( !std::is_same< TemperatureToVelocityRHSOperatorType_, hyteg::NoOperator >::value )
      {
         temperatureToVelocityRHSOperator_->apply(
             temperature, dst.uvw(), level, flag, ( firstOperatorApplyVelocity ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApplyVelocity = false;
      }

      // ---------------- Pressure ----------------

      // Velocity To Pressure
      if constexpr ( !std::is_same< VelocityToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if constexpr ( ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) &&
                        ( preProjectVelocity_ ) && ( !allowPreProjectionToChangeVelocitySrc_ ) )
         {
            velocityToPressureRHSOperator_->apply(
                *velocityApply, dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         }
         else
         {
            velocityToPressureRHSOperator_->apply(
                velocity, dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         }

         firstOperatorApplyPressure = false;
      }

      // Pressure To Pressure
      if constexpr ( !std::is_same< PressureToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         if constexpr ( ( preProjectPressure_ ) && ( !allowPreProjectionToChangePressureSrc_ ) )
         {
            pressureToPressureRHSOperator_->apply(
                *pressureApply, dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         }
         else
         {
            pressureToPressureRHSOperator_->apply(
                pressure, dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         }

         firstOperatorApplyPressure = false;
      }

      // Temperature To Pressure
      if constexpr ( !std::is_same< TemperatureToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         temperatureToPressureRHSOperator_->apply(
             temperature, dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApplyPressure = false;
      }

      // Density Derivative To Pressure
      if constexpr ( !std::is_same< TemperatureToPressureRHSOperatorType_, hyteg::NoOperator >::value )
      {
         densityDerivativeToPressureRHSOperator_->apply(
             densityDerivative, dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApplyPressure = false;
      }

      // #####################################
      // ########## Post Projection ##########
      // #####################################

      if constexpr ( postProjectPressure_ )
      {
         hyteg::projectPressureMean( dst.p(), level );
      }

      if constexpr ( ( !std::is_same< VelocityProjectionOperatorType_, hyteg::NoOperator >::value ) && ( postProjectVelocity_ ) )
      {
         projection_->project( dst.uvw(), level, projectionFlag_ );
      }
   }

   // direct get functions for the internal operators
   VelocityToVelocityRHSOperatorType&    getVelocityToVelocity() const { return *velocityToVelocityRHSOperator_; }
   PressureToVelocityRHSOperatorType&    getPressureToVelocity() const { return *pressureToVelocityRHSOperator_; }
   TemperatureToVelocityRHSOperatorType& getTemperatureToVelocity() const { return *temperatureToVelocityRHSOperator_; }

   VelocityToPressureRHSOperatorType&    getVelocityToPressure() const { return *velocityToPressureRHSOperator_; }
   PressureToPressureRHSOperatorType&    getPressureToPressure() const { return *pressureToPressureRHSOperator_; }
   TemperatureToPressureRHSOperatorType& getTemperatureToPressure() const { return *temperatureToPressureRHSOperator_; }

   DensityDerivativeToPressureRHSOperatorType& getDensityDerivativeToPressure() const
   {
      return *densityDerivativeToPressureRHSOperator_;
   }

   // pointer get functions for the internal operators
   std::shared_ptr< VelocityToVelocityRHSOperatorType > getVelocityToVelocityPtr() const
   {
      return velocityToVelocityRHSOperator_;
   }
   std::shared_ptr< PressureToVelocityRHSOperatorType > getPressureToVelocityPtr() const
   {
      return pressureToVelocityRHSOperator_;
   }
   std::shared_ptr< TemperatureToVelocityRHSOperatorType > getTemperatureToVelocityPtr() const
   {
      return temperatureToVelocityRHSOperator_;
   }

   std::shared_ptr< VelocityToPressureRHSOperatorType > getVelocityToPressurePtr() const
   {
      return velocityToPressureRHSOperator_;
   }
   std::shared_ptr< PressureToPressureRHSOperatorType > getPressureToPressurePtr() const
   {
      return pressureToPressureRHSOperator_;
   }
   std::shared_ptr< TemperatureToPressureRHSOperatorType > getTemperatureToPressurePtr() const
   {
      return temperatureToPressureRHSOperator_;
   }
   std::shared_ptr< DensityDerivativeToPressureRHSOperatorType > getDensityDerivativeToPressurePtr() const
   {
      return densityDerivativeToPressureRHSOperator_;
   }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< VelocityToVelocityRHSOperatorType >    velocityToVelocityRHSOperator_;
   std::shared_ptr< PressureToVelocityRHSOperatorType >    pressureToVelocityRHSOperator_;
   std::shared_ptr< TemperatureToVelocityRHSOperatorType > temperatureToVelocityRHSOperator_;

   std::shared_ptr< VelocityToPressureRHSOperatorType >    velocityToPressureRHSOperator_;
   std::shared_ptr< PressureToPressureRHSOperatorType >    pressureToPressureRHSOperator_;
   std::shared_ptr< TemperatureToPressureRHSOperatorType > temperatureToPressureRHSOperator_;

   std::shared_ptr< DensityDerivativeToPressureRHSOperatorType > densityDerivativeToPressureRHSOperator_;

   std::shared_ptr< VelocityProjectionOperatorType_ > projection_;
   const hyteg::DoFType                               projectionFlag_;

   std::shared_ptr< PressureFunctionType_ > tmpPressure_;
   std::shared_ptr< VelocityFunctionType_ > tmpVelocity_;

   bool lowMemoryMode_;
};

} // namespace MantleConvection
