/*
* Copyright (c) 2024 Andreas Burkhart.
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

#include "hyteg/operators/Operator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

template < class SrcType, class DstType >
class ZeroOperator : public Operator< SrcType, DstType >
{
 public:
   ZeroOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< SrcType, DstType >( storage, minLevel, maxLevel )
   {}

   void apply( const SrcType& src,
               const DstType& dst,
               const uint_t   level,
               const DoFType  flag,
               UpdateType     updateType = Replace ) const
   {
      if ( updateType == Replace )
      {
         dst.interpolate( real_c( 0.0 ), level, flag );
      }
   }
};

} // namespace hyteg