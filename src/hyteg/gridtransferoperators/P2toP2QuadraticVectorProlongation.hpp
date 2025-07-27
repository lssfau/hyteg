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
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

class P2toP2QuadraticVectorProlongation : public ProlongationOperator< P2VectorFunction< real_t > >
{
 public:
   P2toP2QuadraticVectorProlongation() {}

   void prolongate( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongate( function[k], sourceLevel, flag );
      }
   }

   void prolongateAndAdd( const P2VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongateAndAdd( function[k], sourceLevel, flag );
      }
   }

 private:
   P2toP2QuadraticProlongation quadraticProlongationOperator_;
};

template < typename ProjectionOperatorType     = hyteg::NoOperator,
           bool preProject                     = false,
           bool postProject                    = true,
           bool allowPreProjectionToChangeSrc  = true,
           bool allowPostProjectionToChangeDst = true >
class P2toP2QuadraticVectorProlongationWithProjection : public P2toP2QuadraticVectorProlongation
{
 public:
   P2toP2QuadraticVectorProlongationWithProjection( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                                    uint_t                                            minLevel,
                                                    uint_t                                            maxLevel,
                                                    std::shared_ptr< ProjectionOperatorType >         projection = nullptr,
                                                    DoFType projectionFlag = hyteg::FreeslipBoundary,
                                                    bool    lowMemoryMode  = false )
   : P2toP2QuadraticVectorProlongation()
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
            tmp_ = std::make_shared< P2VectorFunction< real_t > >(
                "P2toP2QuadraticVectorRestrictionWithProjection tmp", storage, minLevel, maxLevel );
         }
      }
   }

   void prolongate( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      std::shared_ptr< P2VectorFunction< real_t > > tmpProlongate;

      if constexpr ( useTmpSrc )
      {
         if ( !lowMemoryMode_ )
         {
            tmpProlongate = tmp_;
         }
         else
         {
            tmpProlongate = getTemporaryFunction< P2VectorFunction< real_t > >( storage_, minLevel_, maxLevel_ );
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
         P2toP2QuadraticVectorProlongation::prolongate( *tmpProlongate, sourceLevel, flag );
         function.assign( { real_c( 1 ) }, { *tmpProlongate }, sourceLevel + 1, flag );
      }
      else
      {
         P2toP2QuadraticVectorProlongation::prolongate( function, sourceLevel, flag );
      }

      // postproject
      if constexpr ( ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value ) && ( postProject ) )
      {
         projection_->project( function, sourceLevel + 1, projectionFlag_ );
      }
   }

   void prolongateAndAdd( const P2VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      std::shared_ptr< P2VectorFunction< real_t > > tmpProlongate;
      
      if constexpr ( useTmpSrc || useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmpProlongate = tmp_;
         }
         else
         {
            tmpProlongate = getTemporaryFunction< P2VectorFunction< real_t > >( storage_, minLevel_, maxLevel_ );
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
         P2toP2QuadraticVectorProlongation::prolongateAndAdd( *tmpProlongate, sourceLevel, flag );

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
         P2toP2QuadraticVectorProlongation::prolongateAndAdd( function, sourceLevel, flag );

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
   std::shared_ptr< P2VectorFunction< real_t > > tmp_;

   bool lowMemoryMode_;
};

/***************************************************************************
NOTE: This prolongates the FE function and calls the project function on it 
      so that the normal components are set to zero on the FreeslipBoundary
***************************************************************************/
class P2toP2QuadraticVectorProlongationWithFreeSlipProjection : public P2toP2QuadraticVectorProlongation
{
 public:
   P2toP2QuadraticVectorProlongationWithFreeSlipProjection( std::shared_ptr< P2VectorFunction< real_t > > temp,
                                                            std::shared_ptr< P2ProjectNormalOperator >    projection )
   : P2toP2QuadraticVectorProlongation()
   , temp_( temp )
   , projection_( projection )
   {}

   void prolongate( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      P2toP2QuadraticVectorProlongation::prolongate( function, sourceLevel, flag );
      projection_->project( function, sourceLevel + 1, FreeslipBoundary );
   }

   // prolongateAndAdd has a different implementation than prolongate and seemingly needs internal projections
   // as a quick fix we are using prolongate on a temporary function and add it manually
   void prolongateAndAdd( const P2VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      temp_->assign( { 1.0 }, { function }, sourceLevel, All );
      prolongate( *temp_, sourceLevel, flag );
      function.assign( { 1.0, 1.0 }, { function, *temp_ }, sourceLevel + 1, flag );
   }

 private:
   std::shared_ptr< P2VectorFunction< real_t > > temp_;
   std::shared_ptr< P2ProjectNormalOperator >    projection_;
};

} // namespace hyteg
