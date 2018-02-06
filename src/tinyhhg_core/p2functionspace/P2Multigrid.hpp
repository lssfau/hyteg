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
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol    ,fineRow - 2,sD::EDGE_HO_E )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow - 2,sD::EDGE_VE_N )];
    tmp -= edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1,fineRow - 2,sD::EDGE_HO_E )];

    tmp *= 0.125;

    vertexDofCoarseData[hhg::vertexdof::macroface::indexFromVertex< sourceLevel - 1>(it.col(),it.row(),sD::VERTEX_C)] += tmp;

  }

  /// update edge dof entries
  for ( const auto & it : hhg::edgedof::macroface::Iterator( sourceLevel - 1, 0 ) )
  {
    using hhg::edgedof::macroface::indexFromVertex;

    uint_t fineCol = it.col() * 2;
    uint_t fineRow = it.row() * 2;
    /// horizontal
    if( it.row() != 0) {
      tmp  = 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_DI_NW )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_VE_N  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_HO_NW )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_DI_SE )];
      tmp += 0.5  * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_VE_S  )];
      tmp += 0.25 * edgeDofFineData[indexFromVertex< sourceLevel >(fineCol + 1, fineRow,sD::EDGE_HO_SE )];

    }
    /// diagonal
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( sourceLevel ) - 1)) {

    }
    /// vertical
    if( it.col() != 0) {

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

}/// namespace macroface

namespace macroedge {

template< typename ValueType, uint_t sourceLevel >
void restrictTmpl(const Edge & edge,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & vertexDoFMemoryID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeDoFMemoryID){

  ValueType* dstVertexDofData = edge.getData( vertexDoFMemoryID )->getPointer( sourceLevel - 1 );
  ValueType* srcEdgeDofData    = edge.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );

  typedef hhg::stencilDirection sD;

  for( const auto& it : hhg::vertexdof::macroedge::Iterator(sourceLevel -1,1 )){
    uint_t targetIndex = hhg::vertexdof::macroedge::indexFromVertex< sourceLevel - 1 >(it.col(), sD::VERTEX_C);
    real_t tmp = 0;
    ///south face
    tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_HO_SE)];
    tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_VE_S )];

    tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_HO_SE)];
    tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_DI_SE)];

    tmp +=  3./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_VE_S )];
    ///this weight is adjusted to to precomputation from the edge
    tmp +=  1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_DI_SE)];

    ///on edge
    tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_HO_W )];
    tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_HO_E )];
    tmp +=  3./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_HO_W )];
    tmp +=  3./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_HO_E )];


    if(edge.getNumNeighborFaces() == 2){
      tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_HO_NW)];
      tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 - 1, sD::EDGE_DI_NW)];

      tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_HO_NW)];
      tmp += -1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2 + 1, sD::EDGE_VE_N )];

      tmp +=  3./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_DI_NW)];
      ///this weight is adjusted to to precomputation from the edge
      tmp +=  1./8. * srcEdgeDofData[hhg::edgedof::macroedge::indexFromVertex< sourceLevel >(it.col() * 2    , sD::EDGE_VE_N )];

    }


    dstVertexDofData[targetIndex] += tmp;


  }
}

SPECIALIZE_WITH_VALUETYPE(void, restrictTmpl, restrict)


}/// namespace macroedge

}/// namespace P2
}/// namespace hhg