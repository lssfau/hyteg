#pragma once

#include "tinyhhg_core/StencilMemory.hpp"

namespace hhg{

////////////////////
// Stencil memory //
////////////////////

template< typename ValueType >
using EdgeVertexDoFToEdgeDoFStencilMemory = StencilMemory< ValueType >;

template< typename ValueType >
using FaceVertexDoFToEdgeDoFStencilMemory = StencilMemory< ValueType >;

///
/// \param level stencil size is independent of level
/// \param numDependencies number of adjacent faces of the edge
/// \return number of the stencil entries
inline uint_t EdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  return 2 + numDependencies;
}

inline uint_t FaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( numDependencies );
  return 4;
}

}