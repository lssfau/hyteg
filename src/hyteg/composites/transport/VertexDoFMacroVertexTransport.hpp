/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
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
#include "core/debug/all.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {
namespace vertexdof {
namespace transport {
namespace macrovertex {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using indexing::Index;

template < typename ValueType, bool AlgebraicUpwind >
inline void apply( const uint_t&                                                 level,
                   Vertex&                                                       vertex,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dstId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& uxId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& uyId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& uzId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Vertex >&  xOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Vertex >&  yOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Vertex >&  zOprId )
{
   const ValueType* src = vertex.getData( srcId )->getPointer( level );
   ValueType*       dst = vertex.getData( dstId )->getPointer( level );

   const ValueType* ux = vertex.getData( uxId )->getPointer( level );
   const ValueType* uy = vertex.getData( uyId )->getPointer( level );
   const ValueType* uz = vertex.getData( uzId )->getPointer( level );

   const ValueType* xOperatorData = vertex.getData( xOprId )->getPointer( level );
   const ValueType* yOperatorData = vertex.getData( yOprId )->getPointer( level );
   const ValueType* zOperatorData = vertex.getData( zOprId )->getPointer( level );

   std::vector< ValueType > stencil( vertex.getData( xOprId )->getSize( level ) );
   real_t                   dTmp;

   // fill stencil
   for( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      stencil[i + 1] = 0.5 * ( ux[0] + ux[i + 1] ) * xOperatorData[i + 1];
      stencil[i + 1] += 0.5 * ( uy[0] + uy[i + 1] ) * yOperatorData[i + 1];
      stencil[i + 1] += 0.5 * ( uz[0] + uz[i + 1] ) * zOperatorData[i + 1];
      stencil[0] -= stencil[i + 1];

      // algebraic upwind
      if( AlgebraicUpwind )
      {
         dTmp = std::abs( stencil[i + 1] );
         stencil[0] += dTmp;
         stencil[i + 1] -= dTmp;
      }
   }

   // apply stencil
   dst[0] = stencil[0] * src[0];

   for( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      dst[0] += stencil[i + 1] * src[i + 1];
   }
}

} // namespace macrovertex
} // namespace transport
} // namespace vertexdof
} // namespace hyteg
