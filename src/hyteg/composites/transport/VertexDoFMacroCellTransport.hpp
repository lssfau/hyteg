/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {
namespace vertexdof {
namespace transport {
namespace macrocell {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using indexing::Index;

template < typename ValueType, bool AlgebraicUpwind >
inline void apply( const uint_t&                                               level,
                   Cell&                                                       cell,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dstId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& uxId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& uyId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& uzId,
                   const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >&  xOprId,
                   const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >&  yOprId,
                   const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >&  zOprId )
{
   typedef stencilDirection sd;

   const ValueType* src = cell.getData( srcId )->getPointer( level );
   ValueType*       dst = cell.getData( dstId )->getPointer( level );

   const ValueType* ux = cell.getData( uxId )->getPointer( level );
   const ValueType* uy = cell.getData( uyId )->getPointer( level );
   const ValueType* uz = cell.getData( uzId )->getPointer( level );

   auto xOperatorData = cell.getData( xOprId )->getData( level );
   auto yOperatorData = cell.getData( yOprId )->getData( level );
   auto zOperatorData = cell.getData( zOprId )->getData( level );

   vertexdof::macrocell::StencilMap_T stencil;
   real_t                             dTmp;

   ValueType tmp;
   for( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const uint_t x = it.x();
      const uint_t y = it.y();
      const uint_t z = it.z();

      const uint_t centerIdx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );
      auto centerStencilIdx = vertexdof::logicalIndexOffsetFromVertex( stencilDirection::VERTEX_C );

      // fill stencil
      stencil[centerStencilIdx] = 0.0;

      for( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         auto stencilIdx = vertexdof::logicalIndexOffsetFromVertex( neighbor );
         const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
         stencil[stencilIdx]     = 0.5 * ( ux[centerIdx] + ux[idx] ) * xOperatorData[stencilIdx];
         stencil[stencilIdx] += 0.5 * ( uy[centerIdx] + uy[idx] ) * yOperatorData[stencilIdx];
         stencil[stencilIdx] += 0.5 * ( uz[centerIdx] + uz[idx] ) * zOperatorData[stencilIdx];

         stencil[centerStencilIdx] -= stencil[stencilIdx];
      }

      // algebraic upwind
      if( AlgebraicUpwind )
      {
         for( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
         {
            auto stencilIdx = vertexdof::logicalIndexOffsetFromVertex( neighbor );

            dTmp = std::abs( stencil[stencilIdx] );
            stencil[centerStencilIdx] += dTmp;
            stencil[stencilIdx] -= dTmp;
         }
      }

      tmp = stencil[centerStencilIdx] * src[centerIdx];

      for( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         auto stencilIdx = vertexdof::logicalIndexOffsetFromVertex( neighbor );
         const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
         tmp += stencil[stencilIdx] * src[idx];
      }

      dst[centerIdx] = tmp;
   }
}

} // namespace macrocell
} // namespace transport
} // namespace vertexdof
} // namespace hyteg
