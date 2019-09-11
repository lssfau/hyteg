#pragma once

#include "core/DataTypes.h"
#include "core/debug/all.h"

#include "hyteg/FunctionMemory.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/StencilMemory.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/types/flags.hpp"

namespace hyteg {
namespace vertexdof {
namespace transport {
namespace macroface {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using indexing::Index;

template < typename ValueType, bool AlgebraicUpwind >
inline void apply( const uint_t&                                                                         Level,
                   Face&                                                                                 face,
                   const PrimitiveStorage&                                                               storage,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                           srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                           dstId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                           uxId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                           uyId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                           uzId,
                   const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face >& xOprId,
                   const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face >& yOprId,
                   const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face >& zOprId )
{
   const ValueType* src = face.getData( srcId )->getPointer( Level );
   ValueType*       dst = face.getData( dstId )->getPointer( Level );

   const ValueType* ux = face.getData( uxId )->getPointer( Level );
   const ValueType* uy = face.getData( uyId )->getPointer( Level );
   const ValueType* uz = face.getData( uzId )->getPointer( Level );

   auto xOperatorData = face.getData( xOprId )->getData( Level );
   auto yOperatorData = face.getData( yOprId )->getData( Level );
   auto zOperatorData = face.getData( zOprId )->getData( Level );

   std::map< uint_t, ValueType > stencil;

   if ( face.getNumNeighborCells() == 0 )
   {
      WALBERLA_ABORT( "Not implemented" )
   }

   for ( const auto& idxIt : vertexdof::macroface::Iterator( Level, 1 ) )
   {
      const auto centerArrayIndexOnFace = vertexdof::macroface::index( Level, idxIt.x(), idxIt.y() );

      stencil.clear();

      stencil[centerArrayIndexOnFace] = 0;

      // fill stencil
      for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
      {
         auto neighborCell = storage.getCell( face.neighborCells().at( neighborCellIdx ) );
         auto centerIndexInCell =
             vertexdof::macroface::getIndexInNeighboringMacroCell( idxIt, face, neighborCellIdx, storage, Level );

         for ( const auto& stencilIt : xOperatorData[neighborCellIdx] )
         {
            auto direction = stencilIt.first;

            if ( direction == indexing::IndexIncrement( {0, 0, 0} ) )
               continue;

            auto weightX = xOperatorData[neighborCellIdx][direction];
            auto weightY = yOperatorData[neighborCellIdx][direction];
            auto weightZ = zOperatorData[neighborCellIdx][direction];

            auto leafIndexInMacroCell = centerIndexInCell + direction;
            auto leafIndexInMacroFace = vertexdof::macrocell::getIndexInNeighboringMacroFace(
                leafIndexInMacroCell, *neighborCell, neighborCell->getLocalFaceID( face.getID() ), storage, Level );

            uint_t leafArrayIndexInMacroFace;
            if ( leafIndexInMacroFace.z() == 0 )
            {
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y() );
            }
            else
            {
               WALBERLA_ASSERT_EQUAL( leafIndexInMacroFace.z(), 1 );
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y(), neighborCellIdx );
            }

            if ( stencil.count( leafArrayIndexInMacroFace ) == 0 )
            {
               stencil[leafArrayIndexInMacroFace] = 0;
            }
            stencil[leafArrayIndexInMacroFace] += 0.5 * ( ux[centerArrayIndexOnFace] + ux[leafArrayIndexInMacroFace] ) * weightX;
            stencil[leafArrayIndexInMacroFace] += 0.5 * ( uy[centerArrayIndexOnFace] + uy[leafArrayIndexInMacroFace] ) * weightY;
            stencil[leafArrayIndexInMacroFace] += 0.5 * ( uz[centerArrayIndexOnFace] + uz[leafArrayIndexInMacroFace] ) * weightZ;
            stencil[centerArrayIndexOnFace] -= 0.5 * ( ux[centerArrayIndexOnFace] + ux[leafArrayIndexInMacroFace] ) * weightX;
            stencil[centerArrayIndexOnFace] -= 0.5 * ( uy[centerArrayIndexOnFace] + uy[leafArrayIndexInMacroFace] ) * weightY;
            stencil[centerArrayIndexOnFace] -= 0.5 * ( uz[centerArrayIndexOnFace] + uz[leafArrayIndexInMacroFace] ) * weightZ;
         }
      }

      if ( AlgebraicUpwind )
      {
         for ( auto& it : stencil )
         {
            if ( it.first == centerArrayIndexOnFace )
               continue;

            const auto dTmp = std::abs( it.second );
            stencil[centerArrayIndexOnFace] += dTmp;
            stencil[it.first] -= dTmp;
         }
      }

      real_t tmp = 0;

      // apply stencil
      for ( const auto& it : stencil )
      {
         tmp += src[it.first] * it.second;
      }

      dst[centerArrayIndexOnFace] = tmp;
   }
}

} // namespace macroface
} // namespace transport
} // namespace vertexdof
} // namespace hyteg