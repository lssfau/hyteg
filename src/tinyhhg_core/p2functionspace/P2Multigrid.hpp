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

  ValueType* vertexDofData = face.getData( vertexDoFMemoryID )->getPointer( sourceLevel );
  ValueType* edgeDofData    = face.getData( edgeDoFMemoryID   )->getPointer( sourceLevel );

  typedef hhg::stencilDirection sD;
  sD targetDirection;

  ///put edge dofs which cant be reached from edge onto vertexdof
  ///first edge
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( 0, face.edge_orientation[0] );
  if(face.edge_orientation[0] == 1){
    targetDirection = sD::EDGE_VE_N;
  } else {
    targetDirection = sD::EDGE_DI_NW;
  }
  for( const auto & it : hhg::vertexdof::macroface::BorderIterator( sourceLevel, faceBorderDirection, 0, 1) ){
    ///ignore every second entry
    if(it.col()%2 == 1){
      continue;
    }
    uint_t targetIndex = hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col(),it.row(), targetDirection);

    edgeDofData[targetIndex] = edgeDofData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col()    ,it.row(),sD::EDGE_VE_N)] +
                               edgeDofData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col()    ,it.row(),sD::EDGE_DI_NW)] +
                               edgeDofData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() - 1,it.row(),sD::EDGE_VE_N)] +
                               edgeDofData[hhg::edgedof::macroface::indexFromVertex< sourceLevel >(it.col() - 1,it.row(),sD::EDGE_DI_NW)];

  }

}

SPECIALIZE_WITH_VALUETYPE(void, restrictTmpl, restrict)

}/// namespace macroface

namespace macroedge {

template< typename ValueType, uint_t Level >
void restrictTmpl(const Edge & edge,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & vertexDoFMemoryID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeDoFMemoryID){

}

SPECIALIZE_WITH_VALUETYPE(void, restrictTmpl, restrict)


}/// namespace macroedge

}/// namespace P2
}/// namespace hhg