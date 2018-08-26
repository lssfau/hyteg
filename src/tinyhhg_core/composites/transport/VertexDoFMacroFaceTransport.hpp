#pragma once

#include "core/DataTypes.h"
#include "core/debug/all.h"

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/types/flags.hpp"

namespace hhg {
namespace vertexdof {
namespace transport {
namespace macroface {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using indexing::Index;

template < typename ValueType >
inline void apply( const uint_t&                                               Level,
                   Face&                                                       face,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& uxId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& uyId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& uzId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Face >&  xOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Face >&  yOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Face >&  zOprId )
{
  uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
  uint_t inner_rowsize = rowsize;

  const ValueType* src = face.getData( srcId )->getPointer( Level );
  ValueType*       dst = face.getData( dstId )->getPointer( Level );

  const ValueType* ux = face.getData( uxId )->getPointer( Level );
  const ValueType* uy = face.getData( uyId )->getPointer( Level );
  const ValueType* uz = face.getData( uzId )->getPointer( Level );

  const ValueType* xOperatorData = face.getData( xOprId )->getPointer( Level );
  const ValueType* yOperatorData = face.getData( yOprId )->getPointer( Level );
  const ValueType* zOperatorData = face.getData( zOprId )->getPointer( Level );

  std::vector< ValueType > stencil( 27 );
  real_t                   dTmp;

  ValueType tmp;

  if( face.getNumNeighborCells() == 0 )
  {
    WALBERLA_ABORT("Not implemented")
  }

  for( uint_t j = 1; j < rowsize - 2; ++j )
  {
    for( uint_t i = 1; i < inner_rowsize - 2; ++i )
    {
      // fill stencil
      if( face.getNumNeighborCells() == 1 )
      {
        auto centerIdx = vertexdof::macroface::indexFromVertex( Level, i, j, hhg::stencilDirection::VERTEX_C );
        auto centerStencilIdx = vertexdof::stencilIndexFromVertex( hhg::stencilDirection::VERTEX_C );

        stencil[centerStencilIdx] = 0.0;
        for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
        {
          const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( direction );
          const uint_t idx        = vertexdof::macroface::indexFromVertex( Level, i, j, direction );
          stencil[stencilIdx]     = 0.5 * ( ux[centerIdx] + ux[idx] ) * xOperatorData[stencilIdx];
          stencil[stencilIdx] += 0.5 * ( uy[centerIdx] + uy[idx] ) * yOperatorData[stencilIdx];
          stencil[stencilIdx] += 0.5 * ( uz[centerIdx] + uz[idx] ) * zOperatorData[stencilIdx];
          stencil[centerStencilIdx] -= stencil[stencilIdx];
        }

        // algebraic upwind
        for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
        {
          const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( direction );

          dTmp = std::abs( stencil[stencilIdx] );
          stencil[centerStencilIdx] += dTmp;
          stencil[stencilIdx] -= dTmp;
        }

      } else if( face.getNumNeighborCells() == 2 )
      {
        auto centerIdx = vertexdof::macroface::indexFromVertex( Level, i, j, hhg::stencilDirection::VERTEX_C );
        auto centerStencilIdx = vertexdof::stencilIndexFromVertex( hhg::stencilDirection::VERTEX_C );

        stencil[centerStencilIdx] = 0.0;
        for( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
        {
          const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( direction );
          const uint_t idx        = vertexdof::macroface::indexFromVertex( Level, i, j, direction );
          stencil[stencilIdx]     = 0.5 * ( ux[centerIdx] + ux[idx] ) * xOperatorData[stencilIdx];
          stencil[stencilIdx] += 0.5 * ( uy[centerIdx] + uy[idx] ) * yOperatorData[stencilIdx];
          stencil[stencilIdx] += 0.5 * ( uz[centerIdx] + uz[idx] ) * zOperatorData[stencilIdx];
          stencil[centerStencilIdx] -= stencil[stencilIdx];
        }

        // algebraic upwind
        for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
        {
          const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( direction );

          dTmp = std::abs( stencil[stencilIdx] );
          stencil[centerStencilIdx] += dTmp;
          stencil[stencilIdx] -= dTmp;
        }
      }

      // apply stencil
      if( face.getNumNeighborCells() == 1 )
      {
        tmp = real_c( 0 );
        for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter )
        {
          tmp += stencil[vertexdof::stencilIndexFromVertex( direction )] *
              src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
        }
      } else if( face.getNumNeighborCells() == 2 )
      {
        tmp = real_c( 0 );
        for( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter )
        {
          tmp += stencil[vertexdof::stencilIndexFromVertex( direction )] *
              src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
        }
      }

      WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

      dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
    }
    --inner_rowsize;
  }
}

} // namespace macroface
} // namespace transport
} // namespace vertexdof
} // namespace hhg