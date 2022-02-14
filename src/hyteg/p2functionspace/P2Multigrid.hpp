/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include <core/DataTypes.h>

#include "hyteg/primitives/all.hpp"
#include "hyteg/StencilDirections.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/memory/FunctionMemory.hpp"

namespace hyteg {
namespace P2 {

using walberla::uint_t;

namespace macroface {

template< typename ValueType >
void prolongate(const uint_t sourceLevel,
                    const Face &face,
                    const PrimitiveDataID < FunctionMemory < ValueType >, Face> & vertexDoFMemoryID,
                    const PrimitiveDataID < FunctionMemory < ValueType >, Face> & edgeDoFMemoryID   ){

  ValueType* vertexDofFineData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel + 1 );
  ValueType* edgeDofFineData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel + 1);
  ValueType* vertexDofCoarseData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel );
  ValueType* edgeDofCoarseData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );

  typedef hyteg::stencilDirection sD;

  /// update vertexdofs from vertexdofs
  for( const auto & it : hyteg::vertexdof::macroface::Iterator( sourceLevel , 1)) {

    using hyteg::vertexdof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    vertexDofFineData[indexFromVertex(sourceLevel + 1,fineCol, fineRow, sD::VERTEX_C)] =
      vertexDofCoarseData[indexFromVertex(sourceLevel,it.col(), it.row(), sD::VERTEX_C)];
  }

  /// update vertexdofs from edgedofs
  for( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel , 0)) {

    using hyteg::vertexdof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    if(fineRow != 0) {
      vertexDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow, sD::VERTEX_C)] =
        edgeDofCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_HO_E)];
    }

    if(fineCol + 1 + fineRow != ( hyteg::levelinfo::num_microedges_per_edge( sourceLevel + 1 ) - 1)) {
      vertexDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow + 1, sD::VERTEX_C)] =
        edgeDofCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)];
    }

    if(fineCol != 0) {
      vertexDofFineData[indexFromVertex(sourceLevel + 1,fineCol, fineRow + 1, sD::VERTEX_C)] =
        edgeDofCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_VE_N)];
    }
  }

  /// update edgedofs
  for( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel , 0)) {
    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    using hyteg::edgedof::macroface::indexFromVertex;

    if(fineRow != 0) {
      /// lower left horizontal edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol, fineRow, sD::EDGE_HO_E)] =
        0.75 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_HO_E)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1, it.row(), sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row(), sD::VERTEX_C)];

      /// lower right horizontal edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow, sD::EDGE_HO_E)] =
        0.75 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_HO_E)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col()   , it.row(), sD::VERTEX_C)] +
        0.375 *  vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1,it.row(), sD::VERTEX_C)];
    }

    /// inner horizontal edge
    edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol , fineRow + 1, sD::EDGE_HO_E)] =
      0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_VE_N )] +
      0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)] +
      0.25 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_HO_E )] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col()   , it.row(), sD::VERTEX_C)] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1,it.row(), sD::VERTEX_C)];


    if(fineCol + 1 + fineRow != ( hyteg::levelinfo::num_microedges_per_edge( sourceLevel + 1 ) - 1)) {
      /// lower outer diagonal edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow, sD::EDGE_DI_NE)] =
        0.75 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row() + 1, sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1, it.row(), sD::VERTEX_C)];

      /// upper outer diagonal edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol , fineRow + 1, sD::EDGE_DI_NE)] =
        0.75 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1, it.row()    , sD::VERTEX_C)] +
        0.375  * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col()    , it.row() + 1, sD::VERTEX_C)];
    }

    /// inner diagonal edge
    edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol , fineRow , sD::EDGE_DI_NE)] =
      0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_VE_N )] +
      0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_HO_E )] +
      0.25 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col()    ,it.row() + 1, sD::VERTEX_C)] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1,it.row()    , sD::VERTEX_C)];

    if(fineCol != 0){
      /// lower vertical edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol, fineRow, sD::EDGE_VE_N)] =
        0.75 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_VE_N)] +
        -0.125 *vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row() + 1, sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row()    , sD::VERTEX_C)];

      /// upper vertical edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol, fineRow + 1, sD::EDGE_VE_N)] =
        0.75 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_VE_N)] +
        -0.125 *vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row()    , sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row() + 1, sD::VERTEX_C)];
    }

    /// inner vertical edge
    edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow , sD::EDGE_VE_N )] =
      0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)] +
      0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_HO_E )] +
      0.25 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_VE_N )] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col()    ,it.row()    , sD::VERTEX_C)] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col()    ,it.row() + 1, sD::VERTEX_C)];

    /// we have to update some edge dof which are contained in the upside down triangles
    if(it.col() + 1 + it.row() != hyteg::levelinfo::num_microvertices_per_edge(sourceLevel) - 1) {
      /// horzitonal edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow + 1, sD::EDGE_HO_E)] =
        0.5 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col() + 1, it.row(), sD::EDGE_VE_N)] +
        0.5 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)] +
        0.25 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row() + 1, sD::EDGE_HO_E)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row() + 1, sD::VERTEX_C)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1, it.row() + 1, sD::VERTEX_C)];

      /// diagonal edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow + 1, sD::EDGE_DI_NE)] =
        0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col() + 1, it.row() + 1, sD::EDGE_VE_S )] +
        0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col() + 1, it.row() + 1, sD::EDGE_HO_W )] +
        0.25 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_DI_NE)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col()    ,it.row() + 1, sD::VERTEX_C)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1,it.row()    , sD::VERTEX_C)];

      /// vertical edge
      edgeDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, fineRow + 1, sD::EDGE_VE_N )] =
        0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col()    , it.row()    , sD::EDGE_DI_NE)] +
        0.5  * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col()    , it.row() + 1, sD::EDGE_HO_E )] +
        0.25 * edgeDofCoarseData[indexFromVertex(sourceLevel ,it.col() + 1, it.row(), sD::EDGE_VE_N )] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1,it.row()    , sD::VERTEX_C)] +
        -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel ,it.col() + 1,it.row() + 1, sD::VERTEX_C)];

    }
  }

}


template< typename ValueType >
void restrict(const uint_t sourceLevel,
                  const Face & face,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & vertexDoFMemoryID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & edgeDoFMemoryID){

  ValueType* vertexDofFineData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel );
  ValueType* edgeDofFineData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );
  ValueType* vertexDofCoarseData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofCoarseData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel - 1);

  typedef hyteg::stencilDirection sD;
  sD targetDirection;

  indexing::FaceBoundaryDirection firstFaceBorderDirection  = indexing::getFaceBoundaryDirection( 0, face.edge_orientation[0] );
  indexing::FaceBoundaryDirection secondFaceBorderDirection = indexing::getFaceBoundaryDirection( 1, face.edge_orientation[1] );
  indexing::FaceBoundaryDirection thirdFaceBorderDirection  = indexing::getFaceBoundaryDirection( 2, face.edge_orientation[2] );

  real_t tmp;

  /// update vertex dof entries
  for( const auto & it : hyteg::vertexdof::macroface::Iterator( sourceLevel -1, 1)){
    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    using hyteg::edgedof::macroface::indexFromVertex;

    tmp = 0;
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 2,fineRow + 1,sD::EDGE_HO_E )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 2,fineRow + 1,sD::EDGE_DI_NE)];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 1,fineRow + 1,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 1,fineRow + 1,sD::EDGE_DI_NE)];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow + 1,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow + 1,sD::EDGE_HO_E )];

    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 2,fineRow    ,sD::EDGE_HO_E )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 2,fineRow    ,sD::EDGE_DI_NE)];
    tmp += edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 1,fineRow    ,sD::EDGE_HO_E )] * 3;
    tmp += edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 1,fineRow    ,sD::EDGE_DI_NE)] * 3;
    tmp += edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow    ,sD::EDGE_VE_N )] * 3;
    tmp += edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow    ,sD::EDGE_HO_E )] * 3;
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1,fineRow    ,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1,fineRow    ,sD::EDGE_HO_E )];

    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 1,fineRow - 1,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol - 1,fineRow - 1,sD::EDGE_HO_E )];
    tmp += edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow - 1,sD::EDGE_VE_N )] * 3;
    tmp += edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow - 1,sD::EDGE_DI_NE)] * 3;
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1,fineRow - 1,sD::EDGE_HO_E )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1,fineRow - 1,sD::EDGE_DI_NE)];

    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow - 2,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol    ,fineRow - 2,sD::EDGE_DI_NE)];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1,fineRow - 2,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1,fineRow - 2,sD::EDGE_DI_NE)];

    tmp *= 0.125;
    tmp += vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol, fineRow, sD::VERTEX_C)];

    vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::VERTEX_C)] = tmp;

  }

  /// update edge dof entries
  for ( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel - 1, 0 ) )
  {
    using hyteg::edgedof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;
    /// horizontal
    if( it.row() != 0) {
      tmp  = vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol +1, fineRow, sD::VERTEX_C)];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_DI_NW )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_VE_N  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_HO_NW )];

      tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_HO_W  )];
      tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_HO_E  )];

      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_DI_SE )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_VE_S  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow,sD::EDGE_HO_SE )];

      edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::EDGE_HO_E)] = tmp;
    }
    /// diagonal
    if( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( sourceLevel - 1 ) - 1)) {
      tmp  = vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1, sD::VERTEX_C)];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_HO_W  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_VE_S  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_DI_SW )];

      tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_DI_NW )];
      tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_DI_SE )];

      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_HO_E  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_VE_N  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1,sD::EDGE_DI_NW )];

      edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::EDGE_DI_NE)] = tmp;
    }
    /// vertical
    if( it.col() != 0) {
      tmp  = vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol, fineRow + 1, sD::VERTEX_C)];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_HO_W  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_DI_NW )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_VE_NW )];

      tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_VE_N  )];
      tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_VE_S  )];

      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_HO_E  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_DI_SE )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel,fineCol, fineRow + 1,sD::EDGE_VE_SE )];

      edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::EDGE_VE_N)] = tmp;
    }
  }

  ///this need to be done in two steps since we would overwrite already updated entries otherwise
  ///put edge dofs which cant be reached from edge onto vertexdof
  ///first edge
  ///use ghost layer as temp storage
  for( const auto & it : hyteg::vertexdof::macroface::BoundaryIterator( sourceLevel, firstFaceBorderDirection, 0, 1) ) {
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), sD::VERTEX_C);

    vertexDofFineData[vertexIndex] = edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col()    ,it.row() + 1,sD::EDGE_VE_N)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col()    ,it.row() + 1,sD::EDGE_DI_NW)]+
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() - 1,it.row() + 1,sD::EDGE_VE_N)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() - 1,it.row() + 1,sD::EDGE_DI_NW)];
  }

  ///second edge
  for( const auto & it : hyteg::vertexdof::macroface::BoundaryIterator( sourceLevel, secondFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), sD::VERTEX_C);

    vertexDofFineData[vertexIndex] = edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col()    ,it.row() - 1,sD::EDGE_VE_S)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col()    ,it.row() - 1,sD::EDGE_HO_W)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() - 1,it.row()    ,sD::EDGE_VE_S)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() - 1,it.row()    ,sD::EDGE_HO_W)];
  }

  ///third edge
  for( const auto & it : hyteg::vertexdof::macroface::BoundaryIterator( sourceLevel, thirdFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), sD::VERTEX_C);

    vertexDofFineData[vertexIndex] = edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() + 1,it.row() - 1,sD::EDGE_DI_SE)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() + 1,it.row() - 1,sD::EDGE_HO_E)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() + 1,it.row()    ,sD::EDGE_DI_SE)] +
                                 edgeDofFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col() + 1,it.row()    ,sD::EDGE_HO_E)];
  }

  ///write data to edgedofs
  if(face.edge_orientation[0] == 1){
    targetDirection = sD::EDGE_VE_N;
  } else {
    targetDirection = sD::EDGE_DI_NW;
  }
  for( const auto & it : hyteg::vertexdof::macroface::BoundaryIterator( sourceLevel, firstFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t targetIndex = hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] *= 3.;
    edgeDofFineData[targetIndex] -= vertexDofFineData[vertexIndex];
  }

  if(face.edge_orientation[1] == 1){
    targetDirection = sD::EDGE_HO_W;
  } else {
    targetDirection = sD::EDGE_VE_S;
  }
  for( const auto & it : hyteg::vertexdof::macroface::BoundaryIterator( sourceLevel, secondFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t targetIndex = hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] *= 3.;
    edgeDofFineData[targetIndex] -= vertexDofFineData[vertexIndex];
  }

  if(face.edge_orientation[2] == 1){
    targetDirection = sD::EDGE_DI_SE;
  } else {
    targetDirection = sD::EDGE_HO_E;
  }
  for( const auto & it : hyteg::vertexdof::macroface::BoundaryIterator( sourceLevel, thirdFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.row()%2 == 1){
      continue;
    }
    uint_t targetIndex = hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] *= 3.;
    edgeDofFineData[targetIndex] -= vertexDofFineData[vertexIndex];
  }
}


template< typename ValueType >
void postRestrict(const uint_t sourceLevel,
                      const Face & face,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face > & vertexDoFMemoryID,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face > & edgeDoFMemoryID) {
  typedef hyteg::stencilDirection sD;
  sD targetDirection;

  indexing::FaceBoundaryDirection firstFaceBorderDirection = indexing::getFaceBoundaryDirection( 0, face.edge_orientation[0] );
  indexing::FaceBoundaryDirection secondFaceBorderDirection = indexing::getFaceBoundaryDirection( 1, face.edge_orientation[1] );
  indexing::FaceBoundaryDirection thirdFaceBorderDirection = indexing::getFaceBoundaryDirection( 2, face.edge_orientation[2] );

  ValueType *vertexDofFineData = face.getData(vertexDoFMemoryID)->getPointer(sourceLevel);
  ValueType *edgeDofFineData = face.getData(edgeDoFMemoryID)->getPointer(sourceLevel);

  if (face.edge_orientation[0] == 1) {
    targetDirection = sD::EDGE_VE_N;
  } else {
    targetDirection = sD::EDGE_DI_NW;
  }
  for (const auto &it : hyteg::vertexdof::macroface::BoundaryIterator(sourceLevel, firstFaceBorderDirection, 0, 1)) {
    ///ignore every second entry
    if (it.col() % 2 == 1) {
      continue;
    }
    uint_t targetIndex = hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col(), it.row(), targetDirection);
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(), it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] += vertexDofFineData[vertexIndex];
    edgeDofFineData[targetIndex] /= 3.;
  }

  if (face.edge_orientation[1] == 1) {
    targetDirection = sD::EDGE_HO_W;
  } else {
    targetDirection = sD::EDGE_VE_S;
  }
  for (const auto &it : hyteg::vertexdof::macroface::BoundaryIterator(sourceLevel, secondFaceBorderDirection, 0, 1)) {
    ///ignore every second entry
    if (it.col() % 2 == 1) {
      continue;
    }
    uint_t targetIndex = hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col(), it.row(), targetDirection);
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(), it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] += vertexDofFineData[vertexIndex];
    edgeDofFineData[targetIndex] /= 3.;

  }

  if (face.edge_orientation[2] == 1) {
    targetDirection = sD::EDGE_DI_SE;
  } else {
    targetDirection = sD::EDGE_HO_E;
  }
  for (const auto &it : hyteg::vertexdof::macroface::BoundaryIterator(sourceLevel, thirdFaceBorderDirection, 0, 1)) {
    ///ignore every second entry
    if (it.row() % 2 == 1) {
      continue;
    }
    uint_t targetIndex = hyteg::edgedof::macroface::indexFromVertex(sourceLevel,it.col(), it.row(), targetDirection);
    uint_t vertexIndex = hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,it.col(), it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] += vertexDofFineData[vertexIndex];
    edgeDofFineData[targetIndex] /= 3.;
  }
}

template< typename ValueType >
void restrictInjection(const uint_t sourceLevel,
                           const Face & face,
                           const PrimitiveDataID< FunctionMemory< ValueType >, Face > & vertexDoFMemoryID,
                           const PrimitiveDataID< FunctionMemory< ValueType >, Face > & edgeDoFMemoryID){

  ValueType* vertexDofFineData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel );
  ValueType* vertexDofCoarseData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofCoarseData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel - 1);

  typedef hyteg::stencilDirection sD;

  real_t tmp;

  /// update vertex dof entries
  for( const auto & it : hyteg::vertexdof::macroface::Iterator( sourceLevel -1, 1)){
    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    using hyteg::edgedof::macroface::indexFromVertex;

    tmp = vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol, fineRow, sD::VERTEX_C)];
    vertexDofCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::VERTEX_C)] = ValueType( tmp );
  }

  /// update edge dof entries
  for ( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel - 1, 0 ) )
  {
    using hyteg::edgedof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;
    /// horizontal
    if( it.row() != 0) {
      tmp  = vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol +1, fineRow, sD::VERTEX_C)];
      edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::EDGE_HO_E)] = ValueType( tmp );
    }
    /// diagonal
    if( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( sourceLevel - 1 ) - 1)) {
      tmp  = vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol + 1, fineRow + 1, sD::VERTEX_C)];
      edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::EDGE_DI_NE)] = ValueType( tmp );
    }
    /// vertical
    if( it.col() != 0) {
      tmp  = vertexDofFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,fineCol, fineRow + 1, sD::VERTEX_C)];
      edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),it.row(),sD::EDGE_VE_N)] = ValueType( tmp );
    }
  }

}


}// namespace macroface

namespace macroedge {

template< typename ValueType >
void prolongate(const uint_t sourceLevel,
                    const Edge & edge,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & vertexDoFMemoryID,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeDoFMemoryID) {
  ValueType *vertexDofFineData = edge.getData(vertexDoFMemoryID)->getPointer(sourceLevel + 1);
  ValueType *edgeDofFineData = edge.getData(edgeDoFMemoryID)->getPointer(sourceLevel + 1);
  ValueType *vertexDofCoarseData = edge.getData(vertexDoFMemoryID)->getPointer(sourceLevel);
  ValueType *edgeDofCoarseData = edge.getData(edgeDoFMemoryID)->getPointer(sourceLevel);

  typedef hyteg::stencilDirection sD;

/// update vertexdofs from vertexdofs
  for (const auto &it : hyteg::vertexdof::macroedge::Iterator(sourceLevel, 1)) {

    using hyteg::vertexdof::macroedge::indexFromVertex;

    uint_t fineCol = it.col() * 2;

    vertexDofFineData[indexFromVertex(sourceLevel + 1,fineCol, sD::VERTEX_C)] =
      vertexDofCoarseData[indexFromVertex(sourceLevel ,it.col(), sD::VERTEX_C)];

    /// since the iterator does not include the fist and last vertex we need to handle the first edge separately
    if(it.col() == 1){
      vertexDofFineData[indexFromVertex(sourceLevel + 1,fineCol - 1, sD::VERTEX_C)] =
        edgeDofCoarseData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel ,it.col(), sD::EDGE_HO_W)];
    }

    vertexDofFineData[indexFromVertex(sourceLevel + 1,fineCol + 1, sD::VERTEX_C)] =
      edgeDofCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel ,it.col(), it.row(), sD::EDGE_HO_E)];
  }

  for (const auto &it : hyteg::edgedof::macroedge::Iterator(sourceLevel, 0)) {
    uint_t fineCol = it.col() * 2;

    /// left horizontal edge
    edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel + 1,fineCol, sD::EDGE_HO_E)] =
      0.75 * edgeDofCoarseData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel ,it.col(), sD::EDGE_HO_E)] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel ,it.col() + 1, sD::VERTEX_C)] +
      0.375 * vertexDofCoarseData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel ,it.col(), sD::VERTEX_C)];

    /// right horizontal edge
    edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel + 1,fineCol + 1, sD::EDGE_HO_E)] =
      0.75 * edgeDofCoarseData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel ,it.col(), sD::EDGE_HO_E)] +
      -0.125 * vertexDofCoarseData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel ,it.col(), sD::VERTEX_C)] +
      0.375 * vertexDofCoarseData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel ,it.col() + 1, sD::VERTEX_C)];
  }

}

template< typename ValueType >
void restrict(const uint_t sourceLevel,
                  const Edge & edge,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & vertexDoFMemoryID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeDoFMemoryID){

  ValueType* vertexDofCoarseData = edge.getData( vertexDoFMemoryID )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofCoarseData   = edge.getData( edgeDoFMemoryID   )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofFineData     = edge.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );
  ValueType* vertexDofFineData   = edge.getData( vertexDoFMemoryID )->getPointer( sourceLevel );

  typedef hyteg::stencilDirection sD;
  real_t tmp;

  for( const auto& it : hyteg::vertexdof::macroedge::Iterator(sourceLevel -1,1 )){
    uint_t targetIndex = hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel - 1,it.col(), sD::VERTEX_C);
    tmp = vertexDofFineData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel,it.col() * 2, sD::VERTEX_C)];
    ///south face
    tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 - 1, sD::EDGE_HO_SE)];
    tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 - 1, sD::EDGE_VE_S )];

    tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 + 1, sD::EDGE_HO_SE)];
    tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 + 1, sD::EDGE_DI_SE)];

    tmp +=  3./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2    , sD::EDGE_VE_S )];
    ///this weight is adjusted to to precomputation from the edge
    tmp +=  1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2    , sD::EDGE_DI_SE)];

    ///on edge
    tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 - 1, sD::EDGE_HO_W )];
    tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 + 1, sD::EDGE_HO_E )];
    tmp +=  3./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2    , sD::EDGE_HO_W )];
    tmp +=  3./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2    , sD::EDGE_HO_E )];


    if(edge.getNumNeighborFaces() == 2){
      tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 - 1, sD::EDGE_HO_NW)];
      tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 - 1, sD::EDGE_DI_NW)];

      tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 + 1, sD::EDGE_HO_NW)];
      tmp += -1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2 + 1, sD::EDGE_VE_N )];

      tmp +=  3./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2    , sD::EDGE_DI_NW)];
      ///this weight is adjusted to to precomputation from the edge
      tmp +=  1./8. * edgeDofFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel,it.col() * 2    , sD::EDGE_VE_N )];

    }
    vertexDofCoarseData[targetIndex] = tmp;
  }

  for( const auto& it : hyteg::edgedof::macroedge::Iterator(sourceLevel -1,0 )){
    using hyteg::edgedof::macroedge::indexFromVertex;

    tmp  = vertexDofFineData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel ,it.col() * 2 + 1, sD::VERTEX_C)];
    tmp += 0.5 * edgeDofFineData[indexFromVertex(sourceLevel ,it.col() * 2 + 1, sD::EDGE_DI_SE)];
    tmp += 0.5 * edgeDofFineData[indexFromVertex(sourceLevel ,it.col() * 2 + 1, sD::EDGE_VE_S)];
    tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel ,it.col() * 2 + 1, sD::EDGE_HO_SE)];

    tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,it.col() * 2 + 1,sD::EDGE_HO_W  )];
    tmp += 0.75 * edgeDofFineData[indexFromVertex(sourceLevel,it.col() * 2 + 1,sD::EDGE_HO_E  )];

    if(edge.getNumNeighborFaces() == 2) {
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,it.col() * 2 + 1,sD::EDGE_DI_NW )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex(sourceLevel,it.col() * 2 + 1,sD::EDGE_VE_N  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex(sourceLevel,it.col() * 2 + 1,sD::EDGE_HO_NW )];
    }

    edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),sD::EDGE_HO_E)] = tmp;
  }
}

template< typename ValueType >
void restrictInjection(const uint_t sourceLevel,
                           const Edge & edge,
                           const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & vertexDoFMemoryID,
                           const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeDoFMemoryID){

  ValueType* vertexDofCoarseData = edge.getData( vertexDoFMemoryID )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofCoarseData   = edge.getData( edgeDoFMemoryID   )->getPointer( sourceLevel - 1 );
  ValueType* vertexDofFineData   = edge.getData( vertexDoFMemoryID )->getPointer( sourceLevel );

  typedef hyteg::stencilDirection sD;
  real_t tmp;

  for( const auto& it : hyteg::vertexdof::macroedge::Iterator(sourceLevel -1,1 )){
    uint_t targetIndex = hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel - 1,it.col(), sD::VERTEX_C);
    tmp = vertexDofFineData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel,it.col() * 2, sD::VERTEX_C)];

    vertexDofCoarseData[targetIndex] = ValueType( tmp );
  }

  for( const auto& it : hyteg::edgedof::macroedge::Iterator(sourceLevel -1,0 )){
    using hyteg::edgedof::macroedge::indexFromVertex;

    tmp  = vertexDofFineData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel ,it.col() * 2 + 1, sD::VERTEX_C)];

    edgeDofCoarseData[indexFromVertex(sourceLevel - 1,it.col(),sD::EDGE_HO_E)] = ValueType( tmp );
  }
}


}// namespace macroedge

namespace macrovertex {

template<typename ValueType >
void prolongate(const uint_t sourceLevel,
                    const Vertex &vertex,
                    const PrimitiveDataID <FunctionMemory<ValueType>, Vertex> &vertexDoFMemoryID,
                    const PrimitiveDataID <FunctionMemory<ValueType>, Vertex> &edgeDoFMemoryID) {
  ValueType *vertexDofFineData = vertex.getData(vertexDoFMemoryID)->getPointer(sourceLevel + 1);
  ValueType *vertexDofCoarseData = vertex.getData(vertexDoFMemoryID)->getPointer(sourceLevel);

  vertexDofFineData[0] = vertexDofCoarseData[0];

}


template<typename ValueType >
void restrictInjection(const uint_t sourceLevel,
                           const Vertex &vertex,
                           const PrimitiveDataID <FunctionMemory<ValueType>, Vertex> &vertexDoFMemoryID,
                           const PrimitiveDataID <FunctionMemory<ValueType>, Vertex> &edgeDoFMemoryID) {
  ValueType *vertexDofFineData   = vertex.getData(vertexDoFMemoryID)->getPointer(sourceLevel     );
  ValueType *vertexDofCoarseData = vertex.getData(vertexDoFMemoryID)->getPointer(sourceLevel - 1 );

  vertexDofCoarseData[0] = vertexDofFineData[0];

}

}// namespace macrovertex


}// namespace P2
}// namespace hyteg
