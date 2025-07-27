/*
 * Copyright (c) 2020-2025 Daniel Drzisga, Andreas Burkhart.
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

#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

class P1VectorFunctionLinearProlongation : public ProlongationOperator< P1VectorFunction< real_t > >
{
 public:
   void prolongate( const P1VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         prolongationOperator_.prolongate( function[k], sourceLevel, flag );
      }
   }

   void prolongateAndAdd( const P1VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         prolongationOperator_.prolongateAndAdd( function[k], sourceLevel, flag );
      }
   }

 private:
   P1toP1LinearProlongation< real_t > prolongationOperator_;
};

template < typename ProjectionOperatorType     = hyteg::NoOperator,
           bool preProject                     = false,
           bool postProject                    = true,
           bool allowPreProjectionToChangeSrc  = true,
           bool allowPostProjectionToChangeDst = true >
class P1toP1LinearVectorProlongationWithProjection : public P1VectorFunctionLinearProlongation
{
 public:
   P1toP1LinearVectorProlongationWithProjection( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                                 uint_t                                            minLevel,
                                                 uint_t                                            maxLevel,
                                                 std::shared_ptr< ProjectionOperatorType >         projection = nullptr,
                                                 DoFType projectionFlag                                       = FreeslipBoundary,
                                                 bool    lowMemoryMode                                        = false )
   : P1VectorFunctionLinearProlongation()
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

      if constexpr ( useTmpSrc || useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmp_ = std::make_shared< P1VectorFunction< real_t > >(
                "P1toP1LinearVectorProlongationWithProjection tmp", storage, minLevel, maxLevel );
         }
      }
   }

   void prolongate( const P1VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      std::shared_ptr< P1VectorFunction< real_t > > tmpProlongate;

      if constexpr ( useTmpSrc )
      {
         if ( !lowMemoryMode_ )
         {
            tmpProlongate = tmp_;
         }
         else
         {
            tmpProlongate = getTemporaryFunction< P1VectorFunction< real_t > >( storage_, minLevel_, maxLevel_ );
         }
      }

      // init tmp
      if constexpr ( useTmpSrc )
      {
         tmpProlongate->copyBoundaryConditionFromFunction( function );
         tmpProlongate->assign( { real_c( 1 ) }, { function }, sourceLevel, All );
         // project and restrict are precommunicating -> no communication necessary here
      }

      // preproject
      if constexpr ( ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( preProject ) )
      {
         if constexpr ( !allowPreProjectionToChangeSrc )
         {
            projection_->project( *tmpProlongate, sourceLevel, projectionFlag_ );
         }
         else
         {
            projection_->project( function, sourceLevel, projectionFlag_ );
         }
      }

      // prolongate
      if constexpr ( useTmpSrc )
      {
         P1VectorFunctionLinearProlongation::prolongate( *tmpProlongate, sourceLevel, flag );
         function.assign( { real_c( 1 ) }, { *tmpProlongate }, sourceLevel + 1, flag );
      }
      else
      {
         P1VectorFunctionLinearProlongation::prolongate( function, sourceLevel, flag );
      }

      // postproject
      if constexpr ( ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( postProject ) )
      {
         projection_->project( function, sourceLevel + 1, projectionFlag_ );
      }
   }

   void prolongateAndAdd( const P1VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      std::shared_ptr< P1VectorFunction< real_t > > tmpProlongate;

      if constexpr ( useTmpSrc || useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmpProlongate = tmp_;
         }
         else
         {
            tmpProlongate = getTemporaryFunction< P1VectorFunction< real_t > >( storage_, minLevel_, maxLevel_ );
         }
      }

      // init tmp
      if constexpr ( useTmpSrc || useTmpDst )
      {
         tmpProlongate->copyBoundaryConditionFromFunction( function );
         tmpProlongate->assign( { real_c( 1 ) }, { function }, sourceLevel, All );
         tmpProlongate->assign( { real_c( 1 ) }, { function }, sourceLevel + 1, All );
         // project and restrict are precommunicating -> no communication necessary here
      }

      // preproject
      if constexpr ( preProject )
      {
         if constexpr ( useTmpSrc || useTmpDst )
         {
            projection_->project( *tmpProlongate, sourceLevel, projectionFlag_ );
         }
         else
         {
            projection_->project( function, sourceLevel, projectionFlag_ );
         }
      }

      // prolongateAndAdd
      if constexpr ( useTmpSrc || useTmpDst )
      {
         P1VectorFunctionLinearProlongation::prolongateAndAdd( *tmpProlongate, sourceLevel, flag );

         // postproject
         if constexpr ( ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( postProject ) )
         {
            if ( !allowPostProjectionToChangeDst )
            {
               projection_->project( *tmpProlongate, sourceLevel + 1, projectionFlag_ );
               function.assign( { real_c( 1 ), real_c( 1 ) }, { function, *tmpProlongate }, sourceLevel + 1, flag );
            }
            else
            {
               function.assign( { real_c( 1 ), real_c( 1 ) }, { function, *tmpProlongate }, sourceLevel + 1, flag );
               projection_->project( function, sourceLevel + 1, projectionFlag_ );
            }
         }
      }
      else
      {
         P1VectorFunctionLinearProlongation::prolongateAndAdd( function, sourceLevel, flag );

         // postproject
         if constexpr ( ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( postProject ) )
         {
            projection_->project( function, sourceLevel + 1, projectionFlag_ );
         }
      }
   }

   static constexpr bool useTmpSrc = ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( preProject ) &&
                                     ( !allowPreProjectionToChangeSrc );
   static constexpr bool useTmpDst = ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( postProject ) &&
                                     ( !allowPostProjectionToChangeDst );

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< ProjectionOperatorType >     projection_;
   DoFType                                       projectionFlag_;
   std::shared_ptr< P1VectorFunction< real_t > > tmp_;

   bool lowMemoryMode_;
};

} // namespace hyteg