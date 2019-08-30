#pragma once

#include <tinyhhg_core/p1functionspace/VertexDoFFunction.hpp>

namespace hyteg {

template< typename ValueType >
using P1Function = vertexdof::VertexDoFFunction< ValueType >;

}
