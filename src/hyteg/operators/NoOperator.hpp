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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

class NoOperator
{
 public:
   NoOperator() {}

   NoOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel ) {}

   template < class SrcType, class DstType >
   void apply( const SrcType&       src,
               const DstType&       dst,
               const uint_t         level,
               const hyteg::DoFType flag,
               hyteg::UpdateType    updateType = Replace ) const
   {}

   template < class SrcType, class DstType >
   void applyScaled( const SrcType&       src,
                     const DstType&       dst,
                     real_t               OperatorScaling,
                     const uint_t         level,
                     const hyteg::DoFType flag,
                     hyteg::UpdateType    updateType = Replace ) const
   {}

   template < class SrcType, class DstType >
   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                  const SrcType&                                     src,
                  const DstType&                                     dst,
                  size_t                                             level,
                  hyteg::DoFType                                     flag ) const
   {}

   template < class SrcType >
   void project( const SrcType& dst, size_t level, hyteg::DoFType flag ) const
   {}

   template < class OperatorType >
   void reassembleMatrix( const OperatorType& A )
   {}
};

} // namespace hyteg