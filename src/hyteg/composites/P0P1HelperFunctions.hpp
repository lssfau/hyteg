/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/types/PointND.hpp"

using walberla::real_c;
using walberla::real_t;

namespace hyteg {

void integrateP0P1ToP1( P0Function< real_t >& srcP0,
                        P1Function< real_t >& srcP1,
                        P1Function< real_t >& dstP1,
                        uint_t                level,
                        DoFType               flag )
{
   using namespace vertexdof::macroface;

   // We assume that we always want to treat nodes interior to a macro-face
   WALBERLA_ASSERT( flag == All || flag == Inner );

   // src functions need to be up-to-date
   communication::syncFunctionBetweenPrimitives( srcP1, level );

   // all three functions are assumed to have the same storage or that at least
   // the ones for srcP0 and srcP1 are supersets of that of dstP1
   WALBERLA_ASSERT( srcP0.getStorage() == srcP1.getStorage() );
   WALBERLA_ASSERT( srcP0.getStorage() == dstP1.getStorage() );
   auto storage = dstP1.getStorage();

   // for additive assembly the target must first be zeroed
   dstP1.setToZero( level );

   // loop over all macro-faces to locally integrate to P1
   for ( auto& it : storage->getFaces() )
   {
      Face&      face   = *it.second;
      const auto faceID = face.getID();

      // obtain memory access
      auto srcP1dofs = face.getData( srcP1.getFaceDataID() )->getPointer( level );
      auto dstP1dofs = face.getData( dstP1.getFaceDataID() )->getPointer( level );

      const auto srcP0MemLayout = srcP0.getDGFunction()->volumeDoFFunction()->memoryLayout();
      auto       srcP0dofs      = srcP0.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );

      real_t weightedFaceArea = std::pow( real_c( 4 ), -real_c( level ) ) * face.getArea() / real_c( 3 );

      // got blue and gray triangles
      for ( uint_t microVolType = 0; microVolType < 2; microVolType++ )
      {
         auto faceType = facedof::allFaceTypes[microVolType];

         // loop of all micro-faces of current type belonging to current macro-face
         for ( auto itFace = facedof::macroface::Iterator( level, faceType, 0 ).begin(); itFace != itFace.end(); ++itFace )
         {
            indexing::Index elementIdx   = *itFace;
            auto            neighborInfo = volumedofspace::indexing::ElementNeighborInfo(
                elementIdx, faceType, level, srcP0.getBoundaryCondition(), faceID, storage );

            // loop over the three micro-edges
            for ( uint_t k = 0; k < 3; ++k )
            {
               // if edge is on macro boundary we need to check the DoFType flag
               if ( neighborInfo.atMacroBoundary( k ) && !testFlag( neighborInfo.neighborBoundaryType( k ), flag ) )
               {
                  continue;
               }
               const std::vector< indexing::Index > vertexID = neighborInfo.interfaceVertexIndices( k );

               // compute contribution to final value at both endpoints
               real_t aux = real_c( 0.25 ) *
                            ( srcP1dofs[indexFromVertex( level, vertexID[0].x(), vertexID[0].y(), stencilDirection::VERTEX_C )] +
                              srcP1dofs[indexFromVertex( level, vertexID[1].x(), vertexID[1].y(), stencilDirection::VERTEX_C )] );

               aux *= weightedFaceArea * srcP0dofs[volumedofspace::indexing::index(
                                             elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, srcP0MemLayout )];

               // add value to endpoints
               dstP1dofs[indexFromVertex( level, vertexID[0].x(), vertexID[0].y(), stencilDirection::VERTEX_C )] += aux;
               dstP1dofs[indexFromVertex( level, vertexID[1].x(), vertexID[1].y(), stencilDirection::VERTEX_C )] += aux;
            }
         }
      }

      // need to additively communicate the resuls
      dstP1.communicateAdditively< Face, Edge >( level, true );
      dstP1.communicateAdditively< Face, Vertex >( level, true );
   }
}

} // namespace hyteg
