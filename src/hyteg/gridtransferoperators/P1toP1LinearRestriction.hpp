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
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

template < typename ValueType = real_t >
class P1toP1LinearRestriction : public RestrictionOperator< P1Function< ValueType > >
{
 public:
   void restrict( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      if ( function.getStorage()->hasGlobalCells() )
      {
         restrict3D( function, sourceLevel, flag );
      }
      else
      {
         restrict2DAdditively( function, sourceLevel, flag );
      }
   }

 private:
   void restrict2D( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
   void restrict2DAdditively( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const;

   void restrict3D( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const;

   void restrictMacroVertex( const ValueType* src,
                             ValueType*       dst,
                             const uint_t&    sourceLevel,
                             const uint_t&    numNeighborEdges ) const;

   void
       restrictMacroEdge( const ValueType* src, ValueType* dst, const uint_t& sourceLevel, const uint_t& numNeighborFaces ) const;

   void
       restrictMacroFace( const ValueType* src, ValueType* dst, const uint_t& sourceLevel, const uint_t& numNeighborCells ) const;
};

template < typename ValueType                 = real_t,
           bool preProject                    = false,
           bool postProject                   = true,
           bool allowPreProjectionToChangeSrc = true >
class P1toP1LinearRestrictionWithProjection : public P1toP1LinearRestriction< ValueType >
{
 public:
   P1toP1LinearRestrictionWithProjection( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                          uint_t                                            minLevel,
                                          uint_t                                            maxLevel,
                                          bool                                              lowMemoryMode = false )
   : P1toP1LinearRestriction< ValueType >()
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , lowMemoryMode_( lowMemoryMode )
   {
      if constexpr ( useTmp )
      {
         if ( !lowMemoryMode_ )
         {
            tmp_ = std::make_shared< P1Function< ValueType > >(
                "P1toP1LinearRestrictionWithProjection tmp", storage, minLevel, maxLevel );
         }
      }
   }
   void restrict( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      std::shared_ptr< P1Function< ValueType > > tmpRestrict;

      if constexpr ( useTmp )
      {
         if ( !lowMemoryMode_ )
         {
            tmpRestrict = tmp_;
         }
         else
         {
            tmpRestrict = getTemporaryFunction< P1Function< ValueType > >( storage_, minLevel_, maxLevel_ );
         }
      }

      // init tmp
      if constexpr ( useTmp )
      {
         tmpRestrict->copyBoundaryConditionFromFunction( function );

         tmpRestrict->assign( { real_c( 1 ) }, { function }, sourceLevel, All );
         // restrict is precommunicating and project does not rely on halos -> no communication necessary here
      }

      // preproject
      if constexpr ( preProject )
      {
         if constexpr ( !allowPreProjectionToChangeSrc )
         {
            hyteg::projectPressureMean( *tmpRestrict, sourceLevel );
         }
         else
         {
            hyteg::projectPressureMean( function, sourceLevel );
         }
      }

      // restrict
      if constexpr ( useTmp )
      {
         P1toP1LinearRestriction< ValueType >::restrict( *tmpRestrict, sourceLevel, flag );
         function.assign( { real_c( 1 ) }, { *tmpRestrict }, sourceLevel - 1, flag );
      }
      else
      {
         P1toP1LinearRestriction< ValueType >::restrict( function, sourceLevel, flag );
      }

      // postproject
      if constexpr ( postProject )
      {
         hyteg::projectPressureMean( function, sourceLevel - 1 );
      }
   }

   static constexpr bool useTmp = ( preProject ) && ( !allowPreProjectionToChangeSrc );

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< P1Function< ValueType > > tmp_;

   bool lowMemoryMode_;
};

} // namespace hyteg
