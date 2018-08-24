#pragma once

#include "core/DataTypes.h"
#include "core/debug/all.h"

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/primitives/Cell.hpp"
#include "tinyhhg_core/types/flags.hpp"

namespace hhg {
namespace vertexdof {
namespace transport {
namespace macrocell {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using indexing::Index;

template < typename ValueType >
inline void apply( const uint_t&                                               level,
                   Cell&                                                       cell,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dstId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& uxId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& uyId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& uzId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Cell >&  xOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Cell >&  yOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Cell >&  zOprId )
{
   typedef stencilDirection sd;

   const ValueType* src = cell.getData( srcId )->getPointer( level );
   ValueType*       dst = cell.getData( dstId )->getPointer( level );

   const ValueType* ux = cell.getData( uxId )->getPointer( level );
   const ValueType* uy = cell.getData( uyId )->getPointer( level );
   const ValueType* uz = cell.getData( uzId )->getPointer( level );

   const ValueType* xOperatorData = cell.getData( xOprId )->getPointer( level );
   const ValueType* yOperatorData = cell.getData( yOprId )->getPointer( level );
   const ValueType* zOperatorData = cell.getData( zOprId )->getPointer( level );

   std::vector< ValueType > stencil( 27 );
   real_t                   dTmp;

   ValueType tmp;
   for( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const uint_t x = it.x();
      const uint_t y = it.y();
      const uint_t z = it.z();

      const uint_t centerIdx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );
      const uint_t centerStencilIdx = vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C );

      // fill stencil
      stencil[centerStencilIdx] = 0.0;

      for( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( neighbor );
         const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
         stencil[stencilIdx]     = 0.5 * ( ux[centerIdx] + ux[idx] ) * xOperatorData[stencilIdx];
         stencil[stencilIdx] += 0.5 * ( uy[centerIdx] + uy[idx] ) * yOperatorData[stencilIdx];
         stencil[stencilIdx] += 0.5 * ( uz[centerIdx] + uz[idx] ) * zOperatorData[stencilIdx];

         stencil[centerStencilIdx] -= stencil[stencilIdx];
      }

      // algebraic upwind
      for( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( neighbor );

         dTmp = std::abs( stencil[stencilIdx] );
         stencil[centerStencilIdx] += dTmp;
         stencil[stencilIdx] -= dTmp;
      }

      tmp = stencil[centerStencilIdx] * src[centerIdx];

      for( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( neighbor );
         const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
         tmp += stencil[stencilIdx] * src[idx];
      }

      dst[centerIdx] = tmp;
   }
}

} // namespace macrocell
} // namespace transport
} // namespace vertexdof
} // namespace hhg