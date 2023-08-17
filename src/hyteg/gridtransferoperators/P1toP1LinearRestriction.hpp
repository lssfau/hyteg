/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

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

} // namespace hyteg
