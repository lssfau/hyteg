#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFunction.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"

#include <string>

namespace hhg
{

  using walberla::uint_t;

  template< typename funcType >
  syncFunctionBetweenPrimitives( const funcType *function, uint_t level ) {
    function->getCommunicator( level )->template communicate< Vertex, Edge   >();
    function->getCommunicator( level )->template communicate< Edge,   Vertex >();
    function->getCommunicator( level )->template communicate< Edge,   Face   >();
    function->getCommunicator( level )->template communicate< Face,   Edge   >();
  }

} // namespace hhg
