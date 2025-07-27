/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Nils Kohl, Marcus Mohr, Andreas Burkhart.
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

#include "hyteg/functions/PressureMeanProjection.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

template < typename ValueType = real_t >
class P1toP1LinearProlongation : public ProlongationOperator< P1Function< ValueType > >
{
 public:
   void prolongate( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      if ( function.getStorage()->hasGlobalCells() )
      {
         prolongate3DAdditively( function, sourceLevel, flag, Replace );
      }
      else
      {
         prolongate2DAdditively( function, sourceLevel, flag, Replace );
      }
   }

   void prolongateAndAdd( const P1Function< ValueType >& function,
                          const walberla::uint_t&        sourceLevel,
                          const DoFType&                 flag ) const override
   {
      if ( function.getStorage()->hasGlobalCells() )
      {
         prolongate3DAdditively( function, sourceLevel, flag, Add );
      }
      else
      {
         prolongate2DAdditively( function, sourceLevel, flag, Add );
      }
   }

 private:
   void prolongate2D( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
   void prolongate2DAdditively( const P1Function< ValueType >& function,
                                const uint_t&                  sourceLevel,
                                const DoFType&                 flag,
                                const UpdateType&              updateType ) const;

   void prolongate3DAdditively( const P1Function< ValueType >& function,
                                const uint_t&                  sourceLevel,
                                const DoFType&                 flag,
                                const UpdateType&              updateType ) const;

   void prolongateMacroVertex2D( const ValueType* src, ValueType* dst, const uint_t& sourceLevel ) const;

   void prolongateMacroEdge2D( const ValueType* src, ValueType* dst, const uint_t& sourceLevel ) const;

   void prolongateMacroFace2D( const ValueType* src, ValueType* dst, const uint_t& sourceLevel ) const;
};

template < typename ValueType                  = real_t,
           bool preProject                     = false,
           bool postProject                    = true,
           bool allowPreProjectionToChangeSrc  = true,
           bool allowPostProjectionToChangeDst = true >
class P1toP1LinearProlongationWithProjection : public P1toP1LinearProlongation< ValueType >
{
 public:
   P1toP1LinearProlongationWithProjection( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                           uint_t                                            minLevel,
                                           uint_t                                            maxLevel,
                                           bool                                              lowMemoryMode = false )
   : P1toP1LinearProlongation< ValueType >()
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , lowMemoryMode_( lowMemoryMode )
   {
      if constexpr ( useTmpSrc || useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmp_ = std::make_shared< P1Function< ValueType > >(
                "P2toP2QuadraticVectorRestrictionWithProjection tmp", storage, minLevel, maxLevel );
         }
      }
   }

   void prolongate( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      std::shared_ptr< P1Function< ValueType > > tmpProlongate;

      if constexpr ( useTmpSrc )
      {
         if ( !lowMemoryMode_ )
         {
            tmpProlongate = tmp_;
         }
         else
         {
            tmpProlongate = getTemporaryFunction< P1Function< ValueType > >( storage_, minLevel_, maxLevel_ );
         }
      }

      // init tmp
      if constexpr ( useTmpSrc )
      {
         tmpProlongate->copyBoundaryConditionFromFunction( function );
         tmpProlongate->assign( { real_c( 1 ) }, { function }, sourceLevel, All );
         // prolongate is precommunicating and project does not rely on halos -> no communication necessary here
      }

      // preproject
      if constexpr ( preProject )
      {
         if constexpr ( !allowPreProjectionToChangeSrc )
         {
            hyteg::projectPressureMean( *tmpProlongate, sourceLevel );
         }
         else
         {
            hyteg::projectPressureMean( function, sourceLevel );
         }
      }

      // prolongate
      if constexpr ( useTmpSrc )
      {
         P1toP1LinearProlongation< ValueType >::prolongate( *tmpProlongate, sourceLevel, flag );
         function.assign( { real_c( 1 ) }, { *tmpProlongate }, sourceLevel + 1, flag );
      }
      else
      {
         P1toP1LinearProlongation< ValueType >::prolongate( function, sourceLevel, flag );
      }

      // postproject
      if constexpr ( postProject )
      {
         hyteg::projectPressureMean( function, sourceLevel + 1 );
      }
   }

   void prolongateAndAdd( const P1Function< ValueType >& function,
                          const walberla::uint_t&        sourceLevel,
                          const DoFType&                 flag ) const override
   {
      std::shared_ptr< P1Function< ValueType > > tmpProlongate;

      if constexpr ( useTmpSrc || useTmpDst )
      {
         if ( !lowMemoryMode_ )
         {
            tmpProlongate = tmp_;
         }
         else
         {
            tmpProlongate = getTemporaryFunction< P1Function< ValueType > >( storage_, minLevel_, maxLevel_ );
         }
      }

      // init tmp
      if constexpr ( useTmpSrc || useTmpDst )
      {
         tmpProlongate->copyBoundaryConditionFromFunction( function );
         tmpProlongate->assign( { real_c( 1 ) }, { function }, sourceLevel, All );
         tmpProlongate->assign( { real_c( 1 ) }, { function }, sourceLevel + 1, All );
         // prolongate is precommunicating and project does not rely on halos -> no communication necessary here
      }

      // preproject
      if constexpr ( preProject )
      {
         if constexpr ( useTmpSrc || useTmpDst )
         {
            hyteg::projectPressureMean( *tmpProlongate, sourceLevel );
         }
         else
         {
            hyteg::projectPressureMean( function, sourceLevel );
         }
      }

      // prolongateAndAdd
      if constexpr ( useTmpSrc || useTmpDst )
      {
         P1toP1LinearProlongation< ValueType >::prolongateAndAdd( *tmpProlongate, sourceLevel, flag );

         // postproject
         if constexpr ( postProject )
         {
            if ( !allowPostProjectionToChangeDst )
            {
               hyteg::projectPressureMean( *tmpProlongate, sourceLevel + 1 );
               function.assign( { real_c( 1 ), real_c( 1 ) }, { function, *tmpProlongate }, sourceLevel + 1, flag );
            }
            else
            {
               function.assign( { real_c( 1 ), real_c( 1 ) }, { function, *tmpProlongate }, sourceLevel + 1, flag );
               hyteg::projectPressureMean( function, sourceLevel + 1 );
            }
         }
      }
      else
      {
         P1toP1LinearProlongation< ValueType >::prolongateAndAdd( function, sourceLevel, flag );

         // postproject
         if constexpr ( postProject )
         {
            hyteg::projectPressureMean( function, sourceLevel + 1 );
         }
      }
   }

   static constexpr bool useTmpSrc = ( preProject ) && ( !allowPreProjectionToChangeSrc );
   static constexpr bool useTmpDst = ( postProject ) && ( !allowPostProjectionToChangeDst );

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< P1Function< ValueType > > tmp_;

   bool lowMemoryMode_;
};

} // namespace hyteg
