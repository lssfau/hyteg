#pragma once

#include <tinyhhg_core/p1functionspace/VertexDoFFunction.hpp>

namespace hhg {

template< typename ValueType >
using P1Function = vertexdof::VertexDoFFunction< ValueType >;

}
