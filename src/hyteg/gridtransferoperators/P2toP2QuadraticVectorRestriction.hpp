/*
 * Copyright (c) 2023-2025 Andreas Burkhart.
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

class P2toP2QuadraticVectorRestriction : public RestrictionOperator< P2VectorFunction< real_t > >
{
 public:
   P2toP2QuadraticVectorRestriction() {}

   void restrict( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         quadraticRestrictionOperator_.restrict( function[k], sourceLevel, flag );
      }
   }

 private:
   P2toP2QuadraticRestriction quadraticRestrictionOperator_;
};

template < typename ProjectionOperatorType    = hyteg::NoOperator,
           bool preProject                    = false,
           bool postProject                   = true,
           bool allowPreProjectionToChangeSrc = true >
class P2toP2QuadraticVectorRestrictionWithProjection : public P2toP2QuadraticVectorRestriction
{
 public:
   P2toP2QuadraticVectorRestrictionWithProjection( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                                   uint_t                                            minLevel,
                                                   uint_t                                            maxLevel,
                                                   std::shared_ptr< ProjectionOperatorType >         projection = nullptr,
                                                   DoFType projectionFlag = FreeslipBoundary,
                                                   bool    lowMemoryMode  = false )
   : P2toP2QuadraticVectorRestriction()
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , projection_( projection )
   , projectionFlag_( projectionFlag )
   , lowMemoryMode_( lowMemoryMode )
   {
      if constexpr ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value )
      {
         if ( projection == nullptr )
         {
            WALBERLA_ABORT( "projection set to nullptr but type is not NoOperator!" );
         }
      }

      if constexpr ( useTmp )
      {
         if ( !lowMemoryMode_ )
         {
            tmp_ = std::make_shared< P2VectorFunction< real_t > >(
                "P2toP2QuadraticVectorRestrictionWithProjection tmp", storage, minLevel, maxLevel );
         }
      }
   }

   void restrict( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      std::shared_ptr< P2VectorFunction< real_t > > tmpRestrict;

      if constexpr ( useTmp )
      {
         if ( !lowMemoryMode_ )
         {
            tmpRestrict = tmp_;
         }
         else
         {
            tmpRestrict = getTemporaryFunction< P2VectorFunction< real_t > >( storage_, minLevel_, maxLevel_ );
         }
      }

      // init tmp
      if constexpr ( useTmp )
      {
         tmpRestrict->copyBoundaryConditionFromFunction( function );

         tmpRestrict->assign( { real_c( 1 ) }, { function }, sourceLevel, All );
         // project and restrict are precommunicating -> no communication necessary here
      }

      // preproject
      if constexpr ( ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( preProject ) )
      {
         if constexpr ( !allowPreProjectionToChangeSrc )
         {
            projection_->project( *tmpRestrict, sourceLevel, projectionFlag_ );
         }
         else
         {
            projection_->project( function, sourceLevel, projectionFlag_ );
         }
      }

      // restrict
      if constexpr ( useTmp )
      {
         P2toP2QuadraticVectorRestriction::restrict( *tmpRestrict, sourceLevel, flag );
         function.assign( { real_c( 1 ) }, { *tmpRestrict }, sourceLevel - 1, flag );
      }
      else
      {
         P2toP2QuadraticVectorRestriction::restrict( function, sourceLevel, flag );
      }

      // postproject
      if constexpr ( ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( postProject ) )
      {
         projection_->project( function, sourceLevel - 1, projectionFlag_ );
      }
   }

   static constexpr bool useTmp = ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( preProject ) &&
                                  ( !allowPreProjectionToChangeSrc );

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< ProjectionOperatorType >     projection_;
   DoFType                                       projectionFlag_;
   std::shared_ptr< P2VectorFunction< real_t > > tmp_;

   bool lowMemoryMode_;
};

/***************************************************************************
NOTE: This restricts the FE function and calls the project function on it 
      so that the normal components are set to zero on the FreeslipBoundary
***************************************************************************/
class P2toP2QuadraticVectorRestrictionWithFreeSlipProjection : public P2toP2QuadraticVectorRestriction
{
 public:
   P2toP2QuadraticVectorRestrictionWithFreeSlipProjection( std::shared_ptr< P2ProjectNormalOperator > projection )
   : P2toP2QuadraticVectorRestriction()
   , projection_( projection )
   {}

   void restrict( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      P2toP2QuadraticVectorRestriction::restrict( function, sourceLevel, flag );
      projection_->project( function, sourceLevel - 1, FreeslipBoundary );
   }

 private:
   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

} // namespace hyteg
