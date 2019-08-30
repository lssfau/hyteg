
#pragma once

#include "hyteg/FunctionMemory.hpp"

namespace hyteg {

template< typename ValueType >
using StencilMemory = FunctionMemory< ValueType >;

template< typename DataType, typename PrimitiveType >
using StencilMemoryDataHandling = FunctionMemoryDataHandling< DataType, PrimitiveType >;

}
