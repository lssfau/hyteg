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

#include "../OperatorTools/ABlockProjectionWrapper.hpp"
#include "../OperatorTools/BBlockProjectionWrapper.hpp"
#include "../OperatorTools/BTBlockProjectionWrapper.hpp"
#include "../OperatorTools/StabilisationProjectionWrapper.hpp"

using hyteg::idx_t;
using walberla::real_t;
using walberla::uint_t;

namespace MantleConvection {

template < class AOperatorType_                         = hyteg::NoOperator,
           class BTOperatorType_                        = hyteg::NoOperator,
           class BOperatorType_                         = hyteg::NoOperator,
           class VelocityProjectionOperatorType_        = hyteg::NoOperator,
           class StabilisationOperatorType_             = hyteg::NoOperator,
           bool preProjectPressure_                     = true,
           bool postProjectPressure_                    = true,
           bool preProjectVelocity_                     = true,
           bool postProjectVelocity_                    = true,
           bool allowPreProjectionToChangePressureSrc_  = true,
           bool allowPostProjectionToChangePressureDst_ = true,
           bool allowPreProjectionToChangeVelocitySrc_  = true,
           bool allowPostProjectionToChangeVelocityDst_ = true,
           bool includeABlockProjectionToMatrix         = true,
           bool includeBBlockProjectionToMatrix         = true,
           bool includeBTBlockProjectionToMatrix        = true,
           class StokesFunctionType_                    = hyteg::P2P1TaylorHoodFunction< real_t >,
           class StokesIdxFunctionType_                 = hyteg::P2P1TaylorHoodFunction< idx_t > >
class SaddlePointOperator : public hyteg::Operator< StokesFunctionType_, StokesFunctionType_ >
{
 public:
   // clang-format off
   typedef hyteg::ScaledOperator< AOperatorType_,             typename StokesFunctionType_::VelocityFunction_T, typename StokesFunctionType_::VelocityFunction_T > ScaledAOperatorType;
   typedef hyteg::ScaledOperator< BOperatorType_,             typename StokesFunctionType_::VelocityFunction_T, typename StokesFunctionType_::PressureFunction_T > ScaledBOperatorType;
   typedef hyteg::ScaledOperator< BTOperatorType_,            typename StokesFunctionType_::PressureFunction_T, typename StokesFunctionType_::VelocityFunction_T > ScaledBTOperatorType;
   typedef hyteg::ScaledOperator< StabilisationOperatorType_, typename StokesFunctionType_::PressureFunction_T, typename StokesFunctionType_::PressureFunction_T > ScaledStabilisationOperatorType;

   typedef ABlockProjectionWrapper< ScaledAOperatorType  , VelocityProjectionOperatorType_, preProjectVelocity_, postProjectVelocity_, allowPreProjectionToChangeVelocitySrc_, allowPostProjectionToChangeVelocityDst_, includeABlockProjectionToMatrix  > AOperatorType;
   typedef BBlockProjectionWrapper< ScaledBOperatorType  , VelocityProjectionOperatorType_, preProjectVelocity_, postProjectPressure_, allowPreProjectionToChangeVelocitySrc_, allowPostProjectionToChangePressureDst_, includeBBlockProjectionToMatrix  > BOperatorType;
   typedef BTBlockProjectionWrapper< ScaledBTOperatorType, VelocityProjectionOperatorType_, preProjectPressure_, postProjectVelocity_, allowPreProjectionToChangePressureSrc_, allowPostProjectionToChangeVelocityDst_, includeBTBlockProjectionToMatrix > BTOperatorType;
   typedef StabilisationProjectionWrapper< ScaledStabilisationOperatorType                , preProjectPressure_, postProjectPressure_, allowPreProjectionToChangePressureSrc_, allowPostProjectionToChangePressureDst_ > StabilisationOperatorType;

   typedef AOperatorType_             AOperatorTypeInternal;
   typedef BOperatorType_             BOperatorTypeInternal;
   typedef BTOperatorType_            BTOperatorTypeInternal;
   typedef StabilisationOperatorType_ StabilisationOperatorTypeInternal;

   typedef StokesFunctionType_        StokesFunctionType;
   typedef StokesIdxFunctionType_     StokesIdxFunctionType;
   // clang-format on

   typedef VelocityProjectionOperatorType_ VelocityProjectionOperatorType;

   SaddlePointOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                        uint_t                                            minLevel,
                        uint_t                                            maxLevel,

                        const std::shared_ptr< AOperatorType_ >&             ABlock        = nullptr,
                        const std::shared_ptr< BOperatorType_ >&             BBlock        = nullptr,
                        const std::shared_ptr< BTOperatorType_ >&            BTBlock       = nullptr,
                        const std::shared_ptr< StabilisationOperatorType_ >& stabilisation = nullptr,

                        real_t ABlockScaling        = real_c( 1.0 ),
                        real_t BBlockScaling        = real_c( 1.0 ),
                        real_t BTBlockScaling       = real_c( 1.0 ),
                        real_t stabilisationScaling = real_c( 1.0 ),

                        bool lowMemoryMode = false,

                        const std::shared_ptr< VelocityProjectionOperatorType_ >& projection     = nullptr,
                        hyteg::DoFType                                            projectionFlag = hyteg::FreeslipBoundary )
   : hyteg::Operator< StokesFunctionType_, StokesFunctionType_ >( storage, minLevel, maxLevel )
   , projection_( projection )
   , projectionFlag_( projectionFlag )
   , lowMemoryMode_( lowMemoryMode )
   {
      // #############################
      // #### Check preconditions ####
      // #############################

      // A Block
      if constexpr ( !std::is_same< AOperatorType_, hyteg::NoOperator >::value )
      {
         if ( ABlock == nullptr )
         {
            WALBERLA_ABORT( "A Block set to nullptr but type is not NoOperator!" );
         }
      }

      // BT Block
      if constexpr ( !std::is_same< BTOperatorType_, hyteg::NoOperator >::value )
      {
         if ( BTBlock == nullptr )
         {
            WALBERLA_ABORT( "BT Block set to nullptr but type is not NoOperator!" );
         }
      }

      // B Block
      if constexpr ( !std::is_same< BOperatorType_, hyteg::NoOperator >::value )
      {
         if ( BBlock == nullptr )
         {
            WALBERLA_ABORT( "B Block set to nullptr but type is not NoOperator!" );
         }
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< StabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         if ( stabilisation == nullptr )
         {
            WALBERLA_ABORT( "Stabilisation Block set to nullptr but type is not NoOperator!" );
         }
      }

      // ##################################
      // #### Define wrapped operators ####
      // ##################################

      auto ScaledABlock_  = std::make_shared< ScaledAOperatorType >( storage, minLevel, maxLevel, ABlock, ABlockScaling );
      auto ScaledBBlock_  = std::make_shared< ScaledBOperatorType >( storage, minLevel, maxLevel, BBlock, BBlockScaling );
      auto ScaledBTBlock_ = std::make_shared< ScaledBTOperatorType >( storage, minLevel, maxLevel, BTBlock, BTBlockScaling );
      auto ScaledStabilisation_ =
          std::make_shared< ScaledStabilisationOperatorType >( storage, minLevel, maxLevel, stabilisation, stabilisationScaling );

      ABlock_        = std::make_shared< AOperatorType >( ScaledABlock_, projection_, projectionFlag_, lowMemoryMode_ );
      BBlock_        = std::make_shared< BOperatorType >( ScaledBBlock_, projection_, projectionFlag_, lowMemoryMode_ );
      BTBlock_       = std::make_shared< BTOperatorType >( ScaledBTBlock_, projection_, projectionFlag_, lowMemoryMode_ );
      stabilisation_ = std::make_shared< StabilisationOperatorType >( ScaledStabilisation_, lowMemoryMode_ );
   }

   void apply( const StokesFunctionType_& src,
               const StokesFunctionType_& dst,
               const uint_t               level,
               const hyteg::DoFType       flag,
               const hyteg::UpdateType    updateType = hyteg::UpdateType::Replace ) const
   {
      // #######################################
      // ###### Init First Operator Apply ######
      // #######################################

      bool firstOperatorApplyVelocity = true;
      bool firstOperatorApplyPressure = true;

      // ######################################
      // ################ Main ################
      // ######################################

      // -------------- Momentum --------------

      // A Block
      if constexpr ( !std::is_same< AOperatorType_, hyteg::NoOperator >::value )
      {
         ABlock_->apply(
             src.uvw(), dst.uvw(), level, flag, ( firstOperatorApplyVelocity ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApplyVelocity = false;
      }

      // BT Block
      if constexpr ( !std::is_same< BTOperatorType_, hyteg::NoOperator >::value )
      {
         BTBlock_->apply( src.p(), dst.uvw(), level, flag, ( firstOperatorApplyVelocity ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApplyVelocity = false;
      }

      // ---------------- Mass ----------------

      // B Block
      if constexpr ( !std::is_same< BOperatorType_, hyteg::NoOperator >::value )
      {
         BBlock_->apply( src.uvw(), dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApplyPressure = false;
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< StabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         stabilisation_->apply(
             src.p(), dst.p(), level, flag, ( firstOperatorApplyPressure ? updateType : hyteg::UpdateType::Add ) );
         firstOperatorApplyPressure = false;
      }
   }

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                  const StokesIdxFunctionType_&                      src,
                  const StokesIdxFunctionType_&                      dst,
                  uint_t                                             level,
                  hyteg::DoFType                                     flag ) const
   {
      // ######################################
      // ################ Main ################
      // ######################################

      // -------------- Momentum --------------

      // A Block
      if constexpr ( !std::is_same< AOperatorType_, hyteg::NoOperator >::value )
      {
         ABlock_->toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      }

      // BT Block
      if constexpr ( !std::is_same< BTOperatorType_, hyteg::NoOperator >::value )
      {
         BTBlock_->toMatrix( mat, src.p(), dst.uvw(), level, flag );
      }

      // ---------------- Mass ----------------

      // B Block
      if constexpr ( !std::is_same< BOperatorType_, hyteg::NoOperator >::value )
      {
         BBlock_->toMatrix( mat, src.uvw(), dst.p(), level, flag );
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< StabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         stabilisation_->toMatrix( mat, src.p(), dst.p(), level, flag );
      }
   }

   void computeInverseDiagonalOperatorValues()
   {
      // A Block
      if constexpr ( !std::is_same< AOperatorType_, hyteg::NoOperator >::value )
      {
         ABlock_->computeInverseDiagonalOperatorValues();
      }

      // Stabilisation Block
      if constexpr ( !std::is_same< StabilisationOperatorType_, hyteg::NoOperator >::value )
      {
         stabilisation_->computeInverseDiagonalOperatorValues();
      }
   }

   // direct get functions for the internal operators
   AOperatorType&                   getA() const { return *ABlock_; }
   BTOperatorType&                  getBT() const { return *BTBlock_; }
   BOperatorType&                   getB() const { return *BBlock_; }
   StabilisationOperatorType&       getStab() const { return *stabilisation_; }
   VelocityProjectionOperatorType_& getProj() const { return *projection_; }
   hyteg::DoFType                   getProjFlag() const { return projectionFlag_; }

   // pointer get functions for the internal operators
   std::shared_ptr< AOperatorType >                   getAPtr() const { return ABlock_; }
   std::shared_ptr< BTOperatorType >                  getBTPtr() const { return BTBlock_; }
   std::shared_ptr< BOperatorType >                   getBPtr() const { return BBlock_; }
   std::shared_ptr< StabilisationOperatorType >       getStabPtr() const { return stabilisation_; }
   std::shared_ptr< VelocityProjectionOperatorType_ > getProjPtr() const { return projection_; }

   // static bool values
   static constexpr bool preProjectPressure  = preProjectPressure_;
   static constexpr bool postProjectPressure = postProjectPressure_;
   static constexpr bool preProjectVelocity  = preProjectVelocity_;
   static constexpr bool postProjectVelocity = postProjectVelocity_;

   static constexpr bool allowPreProjectionToChangePressureSrc  = allowPreProjectionToChangePressureSrc_;
   static constexpr bool allowPostProjectionToChangePressureDst = allowPostProjectionToChangePressureDst_;
   static constexpr bool allowPreProjectionToChangeVelocitySrc  = allowPreProjectionToChangeVelocitySrc_;
   static constexpr bool allowPostProjectionToChangeVelocityDst = allowPostProjectionToChangeVelocityDst_;

 private:
   std::shared_ptr< AOperatorType >             ABlock_;
   std::shared_ptr< BOperatorType >             BBlock_;
   std::shared_ptr< BTOperatorType >            BTBlock_;
   std::shared_ptr< StabilisationOperatorType > stabilisation_;

   std::shared_ptr< VelocityProjectionOperatorType_ > projection_;
   const hyteg::DoFType                               projectionFlag_;

   bool lowMemoryMode_;
};

} // namespace MantleConvection
