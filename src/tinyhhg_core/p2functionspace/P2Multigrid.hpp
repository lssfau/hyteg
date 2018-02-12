#pragma once

#include <core/DataTypes.h>

#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/macros.hpp"

namespace hhg {
namespace P2 {

using walberla::uint_t;

namespace macroface {

template< typename ValueType, uint_t sourceLevel >
void prolongateTmpl(const Face &face,
                    const PrimitiveDataID < FunctionMemory < ValueType >, Face> & vertexDoFMemoryID,
                    const PrimitiveDataID < FunctionMemory < ValueType >, Face> & edgeDoFMemoryID   ){

  ValueType* vertexDofFineData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel + 1 );
  ValueType* edgeDofFineData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel + 1);
  ValueType* vertexDofCoarseData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel );
  ValueType* edgeDofCoarseData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );

  typedef hhg::stencilDirection sD;

  /// update vertexdofs from vertexdofs
  for( const auto & it : hhg::vertexdof::macroface::Iterator( sourceLevel , 1)) {

    using hhg::vertexdof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    vertexDofFineData[indexFromVertex< sourceLevel +1 >(fineCol, fineRow, sD::VERTEX_C)] =
      vertexDofCoarseData[indexFromVertex< sourceLevel >(it.col(), it.row(), sD::VERTEX_C)];
  }

  /// update vertexdofs from edgedofs
  for( const auto & it : hhg::edgedof::macroface::Iterator( sourceLevel , 0)) {

    using hhg::vertexdof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    if(fineRow != 0) {
      vertexDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow, sD::VERTEX_C)] =
        edgeDofCoarseData[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_HO_E)];
    }

    if(fineCol + 1 + fineRow != (hhg::levelinfo::num_microedges_per_edge( sourceLevel + 1 ) - 1)) {
      vertexDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow + 1, sD::VERTEX_C)] =
        edgeDofCoarseData[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)];
    }

    if(fineCol != 0) {
      vertexDofFineData[indexFromVertex<sourceLevel + 1>(fineCol, fineRow + 1, sD::VERTEX_C)] =
        edgeDofCoarseData[hhg::edgedof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_VE_N)];
    }
  }

  /// update edgedofs
  for( const auto & it : hhg::edgedof::macroface::Iterator( sourceLevel , 0)) {
    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    using hhg::edgedof::macroface::indexFromVertex;

    if(fineRow != 0) {
      /// lower left horizontal edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol, fineRow, sD::EDGE_HO_E)] =
        0.75 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_HO_E)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1, it.row(), sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row(), sD::VERTEX_C)];

      /// lower right horizontal edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow, sD::EDGE_HO_E)] =
        0.75 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_HO_E)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col()   , it.row(), sD::VERTEX_C)] +
        0.375 *  vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1,it.row(), sD::VERTEX_C)];
    }

    /// inner horizontal edge
    edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol , fineRow + 1, sD::EDGE_HO_E)] =
      0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_VE_N )] +
      0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)] +
      0.25 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_HO_E )] +
      -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col()   , it.row(), sD::VERTEX_C)] +
      -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1,it.row(), sD::VERTEX_C)];


    if(fineCol + 1 + fineRow != (hhg::levelinfo::num_microedges_per_edge( sourceLevel + 1 ) - 1)) {
      /// lower outer diagonal edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow, sD::EDGE_DI_NE)] =
        0.75 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row() + 1, sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1, it.row(), sD::VERTEX_C)];

      /// upper outer diagonal edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol , fineRow + 1, sD::EDGE_DI_NE)] =
        0.75 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1, it.row()    , sD::VERTEX_C)] +
        0.375  * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col()    , it.row() + 1, sD::VERTEX_C)];
    }

    /// inner diagonal edge
    edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol , fineRow , sD::EDGE_DI_NE)] =
      0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_VE_N )] +
      0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_HO_E )] +
      0.25 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)] +
      -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col()    ,it.row() + 1, sD::VERTEX_C)] +
      -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1,it.row()    , sD::VERTEX_C)];

    if(fineCol != 0){
      /// lower vertical edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol, fineRow, sD::EDGE_VE_N)] =
        0.75 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_VE_N)] +
        -0.125 *vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row() + 1, sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row()    , sD::VERTEX_C)];

      /// upper vertical edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol, fineRow + 1, sD::EDGE_VE_N)] =
        0.75 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_VE_N)] +
        -0.125 *vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row()    , sD::VERTEX_C)] +
        0.375 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row() + 1, sD::VERTEX_C)];
    }

    /// inner vertical edge
    edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow , sD::EDGE_VE_N )] =
      0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)] +
      0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_HO_E )] +
      0.25 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_VE_N )] +
      -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col()    ,it.row()    , sD::VERTEX_C)] +
      -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col()    ,it.row() + 1, sD::VERTEX_C)];

    /// we have to update some edge dof which are contained in the upside down triangles
    if(it.col() + 1 + it.row() != hhg::levelinfo::num_microvertices_per_edge(sourceLevel) - 1) {
      /// horzitonal edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow + 1, sD::EDGE_HO_E)] =
        0.5 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col() + 1, it.row(), sD::EDGE_VE_N)] +
        0.5 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)] +
        0.25 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row() + 1, sD::EDGE_HO_E)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col(), it.row() + 1, sD::VERTEX_C)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1, it.row() + 1, sD::VERTEX_C)];

      /// diagonal edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow + 1, sD::EDGE_DI_NE)] =
        0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col() + 1, it.row() + 1, sD::EDGE_VE_S )] +
        0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col() + 1, it.row() + 1, sD::EDGE_HO_W )] +
        0.25 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col(), it.row(), sD::EDGE_DI_NE)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col()    ,it.row() + 1, sD::VERTEX_C)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1,it.row()    , sD::VERTEX_C)];

      /// vertical edge
      edgeDofFineData[indexFromVertex<sourceLevel + 1>(fineCol + 1, fineRow + 1, sD::EDGE_VE_N )] =
        0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col()    , it.row()    , sD::EDGE_DI_NE)] +
        0.5  * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col()    , it.row() + 1, sD::EDGE_HO_E )] +
        0.25 * edgeDofCoarseData[indexFromVertex<sourceLevel>(it.col() + 1, it.row(), sD::EDGE_VE_N )] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1,it.row()    , sD::VERTEX_C)] +
        -0.125 * vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex<sourceLevel>(it.col() + 1,it.row() + 1, sD::VERTEX_C)];

    }
  }

}

SPECIALIZE_WITH_VALUETYPE(void, prolongateTmpl, prolongate)

template< typename ValueType, uint_t sourceLevel >
void restrictTmpl(const Face & face,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & vertexDoFMemoryID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & edgeDoFMemoryID){

  ValueType* vertexDofFineData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel );
  ValueType* edgeDofFineData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );
  ValueType* vertexDofCoarseData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofCoarseData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel - 1);

  typedef hhg::stencilDirection sD;
  sD targetDirection;

  indexing::FaceBorderDirection firstFaceBorderDirection  = indexing::getFaceBorderDirection( 0, face.edge_orientation[0] );
  indexing::FaceBorderDirection secondFaceBorderDirection = indexing::getFaceBorderDirection( 1, face.edge_orientation[1] );
  indexing::FaceBorderDirection thirdFaceBorderDirection  = indexing::getFaceBorderDirection( 2, face.edge_orientation[2] );

  real_t tmp;

  /// update vertex dof entries
  for( const auto & it : hhg::vertexdof::macroface::Iterator( sourceLevel -1, 1)){
    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;

    using hhg::edgedof::macroface::indexFromVertex;

    tmp = 0;
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 2,fineRow + 1,sD::EDGE_HO_E )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 2,fineRow + 1,sD::EDGE_DI_NE)];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 1,fineRow + 1,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 1,fineRow + 1,sD::EDGE_DI_NE)];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow + 1,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow + 1,sD::EDGE_HO_E )];

    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 2,fineRow    ,sD::EDGE_HO_E )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 2,fineRow    ,sD::EDGE_DI_NE)];
    tmp += edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 1,fineRow    ,sD::EDGE_HO_E )] * 3;
    tmp += edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 1,fineRow    ,sD::EDGE_DI_NE)] * 3;
    tmp += edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow    ,sD::EDGE_VE_N )] * 3;
    tmp += edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow    ,sD::EDGE_HO_E )] * 3;
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow    ,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow    ,sD::EDGE_HO_E )];

    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 1,fineRow - 1,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol - 1,fineRow - 1,sD::EDGE_HO_E )];
    tmp += edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow - 1,sD::EDGE_VE_N )] * 3;
    tmp += edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow - 1,sD::EDGE_DI_NE)] * 3;
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow - 1,sD::EDGE_HO_E )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow - 1,sD::EDGE_DI_NE)];

    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow - 2,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow - 2,sD::EDGE_DI_NE)];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow - 2,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow - 2,sD::EDGE_DI_NE)];

    tmp *= 0.125;
    tmp += vertexDofFineData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(fineCol, fineRow, sD::VERTEX_C)];

    vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel - 1>(it.col(),it.row(),sD::VERTEX_C)] = tmp;

  }

  /// update edge dof entries
  for ( const auto & it : hhg::edgedof::macroface::Iterator( sourceLevel - 1, 0 ) )
  {
    using hhg::edgedof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;
    /// horizontal
    if( it.row() != 0) {
      tmp  = vertexDofFineData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(fineCol +1, fineRow, sD::VERTEX_C)];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_DI_NW )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_VE_N  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_HO_NW )];

      tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_HO_W  )];
      tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_HO_E  )];

      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_DI_SE )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_VE_S  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_HO_SE )];

      edgeDofCoarseData[indexFromVertex< sourceLevel -1 >(it.col(),it.row(),sD::EDGE_HO_E)] = tmp;
    }
    /// diagonal
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( sourceLevel - 1 ) - 1)) {
      tmp  = vertexDofFineData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1, sD::VERTEX_C)];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_HO_W  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_VE_S  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_DI_SW )];

      tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_DI_NW )];
      tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_DI_SE )];

      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_HO_E  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_VE_N  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow + 1,sD::EDGE_DI_NW )];

      edgeDofCoarseData[indexFromVertex< sourceLevel -1 >(it.col(),it.row(),sD::EDGE_DI_NE)] = tmp;
    }
    /// vertical
    if( it.col() != 0) {
      tmp  = vertexDofFineData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(fineCol, fineRow + 1, sD::VERTEX_C)];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_HO_W  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_DI_NW )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_VE_NW )];

      tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_VE_N  )];
      tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_VE_S  )];

      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_HO_E  )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_DI_SE )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol, fineRow + 1,sD::EDGE_VE_SE )];

      edgeDofCoarseData[indexFromVertex< sourceLevel -1 >(it.col(),it.row(),sD::EDGE_VE_N)] = tmp;
    }
  }




  ///this need to be done in two steps since we would overwrite already updated entries otherwise
  ///put edge dofs which cant be reached from edge onto vertexdof
  ///first edge
  ///use ghost layer as temp storage
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, firstFaceBorderDirection, 0, 1) ) {
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    vertexDofFineData[vertexIndex] = edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col()    ,it.row() + 1,sD::EDGE_VE_N)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col()    ,it.row() + 1,sD::EDGE_DI_NW)]+
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() - 1,it.row() + 1,sD::EDGE_VE_N)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() - 1,it.row() + 1,sD::EDGE_DI_NW)];
  }

  ///second edge
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, secondFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    vertexDofFineData[vertexIndex] = edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col()    ,it.row() - 1,sD::EDGE_VE_S)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col()    ,it.row() - 1,sD::EDGE_HO_W)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() - 1,it.row()    ,sD::EDGE_VE_S)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() - 1,it.row()    ,sD::EDGE_HO_W)];
  }

  ///third edge
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, thirdFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    vertexDofFineData[vertexIndex] = edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() + 1,it.row() - 1,sD::EDGE_DI_SE)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() + 1,it.row() - 1,sD::EDGE_HO_E)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() + 1,it.row()    ,sD::EDGE_DI_SE)] +
                                 edgeDofFineData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() + 1,it.row()    ,sD::EDGE_HO_E)];
  }

  ///write data to edgedofs
  if(face.edge_orientation[0] == 1){
    targetDirection = sD::EDGE_VE_N;
  } else {
    targetDirection = sD::EDGE_DI_NW;
  }
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, firstFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t targetIndex = hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] *= 3.;
    edgeDofFineData[targetIndex] -= vertexDofFineData[vertexIndex];
  }

  if(face.edge_orientation[1] == 1){
    targetDirection = sD::EDGE_HO_W;
  } else {
    targetDirection = sD::EDGE_VE_S;
  }
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, secondFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t targetIndex = hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] *= 3.;
    edgeDofFineData[targetIndex] -= vertexDofFineData[vertexIndex];
  }

  if(face.edge_orientation[2] == 1){
    targetDirection = sD::EDGE_DI_SE;
  } else {
    targetDirection = sD::EDGE_HO_E;
  }
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, thirdFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.row()%2 == 1){
      continue;
    }
    uint_t targetIndex = hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] *= 3.;
    edgeDofFineData[targetIndex] -= vertexDofFineData[vertexIndex];
  }

}

SPECIALIZE_WITH_VALUETYPE(void, restrictTmpl, restrict)

template< typename ValueType, uint_t sourceLevel >
void postRestrictTmpl(const Face & face,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face > & vertexDoFMemoryID,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face > & edgeDoFMemoryID){

  typedef hhg::stencilDirection sD;
  sD targetDirection;

  indexing::FaceBorderDirection firstFaceBorderDirection  = indexing::getFaceBorderDirection( 0, face.edge_orientation[0] );
  indexing::FaceBorderDirection secondFaceBorderDirection = indexing::getFaceBorderDirection( 1, face.edge_orientation[1] );
  indexing::FaceBorderDirection thirdFaceBorderDirection  = indexing::getFaceBorderDirection( 2, face.edge_orientation[2] );

  ValueType* vertexDofFineData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel );
  ValueType* edgeDofFineData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );

  if(face.edge_orientation[0] == 1){
    targetDirection = sD::EDGE_VE_N;
  } else {
    targetDirection = sD::EDGE_DI_NW;
  }
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, firstFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t targetIndex = hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] += vertexDofFineData[vertexIndex];
    edgeDofFineData[targetIndex] /= 3.;
  }

  if(face.edge_orientation[1] == 1){
    targetDirection = sD::EDGE_HO_W;
  } else {
    targetDirection = sD::EDGE_VE_S;
  }
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, secondFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t targetIndex = hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] += vertexDofFineData[vertexIndex];
    edgeDofFineData[targetIndex] /= 3.;

  }

  if(face.edge_orientation[2] == 1){
    targetDirection = sD::EDGE_DI_SE;
  } else {
    targetDirection = sD::EDGE_HO_E;
  }
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, thirdFaceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.row()%2 == 1){
      continue;
    }
    uint_t targetIndex = hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), targetDirection);
    uint_t vertexIndex = hhg::vertexdof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), sD::VERTEX_C);

    edgeDofFineData[targetIndex] += vertexDofFineData[vertexIndex];
    edgeDofFineData[targetIndex] /= 3.;
  }

}

SPECIALIZE_WITH_VALUETYPE(void, postRestrictTmpl, postRestrict)


}/// namespace macroface

namespace macroedge {

template< typename ValueType, uint_t sourceLevel >
void restrictTmpl(const Edge & edge,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & vertexDoFMemoryID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeDoFMemoryID){

  ValueType* vertexDofCoarseData = edge.getData( vertexDoFMemoryID )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofCoarseData   = edge.getData( edgeDoFMemoryID   )->getPointer( sourceLevel - 1 );
  ValueType* edgeDofFineData     = edge.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );
  ValueType* vertexDofFineData   = edge.getData( vertexDoFMemoryID )->getPointer( sourceLevel );

  typedef hhg::stencilDirection sD;
  real_t tmp;

  for( const auto& it : hhg::vertexdof::macroedge::Iterator(sourceLevel -1,1 )){
    uint_t targetIndex = hhg::vertexdof::macroedge::indexFromVertex< sourceLevel - 1 >(it.col(), sD::VERTEX_C);
    tmp = vertexDofFineData[hhg::vertexdof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2, sD::VERTEX_C)];
    ///south face
    tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_HO_SE)];
    tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_VE_S )];

    tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_HO_SE)];
    tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_DI_SE)];

    tmp +=  3./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_VE_S )];
    ///this weight is adjusted to to precomputation from the edge
    tmp +=  1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_DI_SE)];

    ///on edge
    tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_HO_W )];
    tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_HO_E )];
    tmp +=  3./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_HO_W )];
    tmp +=  3./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_HO_E )];


    if(edge.getNumNeighborFaces() == 2){
      tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_HO_NW)];
      tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_DI_NW)];

      tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_HO_NW)];
      tmp += -1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_VE_N )];

      tmp +=  3./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_DI_NW)];
      ///this weight is adjusted to to precomputation from the edge
      tmp +=  1./8. * edgeDofFineData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_VE_N )];

    }
    vertexDofCoarseData[targetIndex] = tmp;
  }

  for( const auto& it : hhg::edgedof::macroedge::Iterator(sourceLevel -1,0 )){
    using hhg::edgedof::macroedge::indexFromVertex;

    tmp  = vertexDofFineData[hhg::vertexdof::macroedge::indexFromVertex<sourceLevel>(it.col() * 2 + 1, sD::VERTEX_C)];
    tmp += 0.5 * edgeDofFineData[indexFromVertex<sourceLevel>(it.col() * 2 + 1, sD::EDGE_DI_SE)];
    tmp += 0.5 * edgeDofFineData[indexFromVertex<sourceLevel>(it.col() * 2 + 1, sD::EDGE_VE_S)];
    tmp += 0.25 * edgeDofFineData[indexFromVertex<sourceLevel>(it.col() * 2 + 1, sD::EDGE_HO_SE)];

    tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(it.col() * 2 + 1,sD::EDGE_HO_W  )];
    tmp += 0.75 * edgeDofFineData[indexFromVertex< sourceLevel >(it.col() * 2 + 1,sD::EDGE_HO_E  )];

    if(edge.getNumNeighborFaces() == 2) {
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(it.col() * 2 + 1,sD::EDGE_DI_NW )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(it.col() * 2 + 1,sD::EDGE_VE_N  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(it.col() * 2 + 1,sD::EDGE_HO_NW )];
    }

    edgeDofCoarseData[indexFromVertex< sourceLevel -1 >(it.col(),sD::EDGE_HO_E)] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, restrictTmpl, restrict)


}/// namespace macroedge

}/// namespace P2
}/// namespace hhg