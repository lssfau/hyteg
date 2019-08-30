#pragma once

#include <hyteg/p1functionspace/VertexDoFFunction.hpp>

namespace hyteg {

template< typename ValueType >
using P1Function = vertexdof::VertexDoFFunction< ValueType >;

}
