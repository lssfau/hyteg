#pragma once

#include <core/DataTypes.h>

#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/macros.hpp"

namespace hhg {
namespace P2 {

using walberla::uint_t;

namespace macroface {

template< typename ValueType, uint_t Level >
void restrictTmpl(const Face & face,
         const PrimitiveDataID< FunctionMemory< ValueType >, Face > & vertexDoFMemoryID,
         const PrimitiveDataID< FunctionMemory< ValueType >, Face > & edgeDoFMemoryID){

  ///put edge dofs which cant be reached from edge onto vertexdof


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