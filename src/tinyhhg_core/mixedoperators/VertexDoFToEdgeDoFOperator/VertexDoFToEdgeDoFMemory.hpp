#pragma once

#include "tinyhhg_core/StencilMemory.hpp"

namespace hhg{

////////////////////
// Stencil memory //
////////////////////

template< typename ValueType >
using VertexVertexDoFToEdgeDoFStencilMemory = StencilMemory< ValueType >;

template< typename ValueType >
using EdgeVertexDoFToEdgeDoFStencilMemory = StencilMemory< ValueType >;

template< typename ValueType >
using FaceVertexDoFToEdgeDoFStencilMemory = StencilMemory< ValueType >;

}