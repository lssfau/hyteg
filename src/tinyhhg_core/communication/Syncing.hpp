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
  void syncFunctionBetweenPrimitives( const funcType *function, uint_t level ) {
    function->template communicate< Vertex, Edge   >( level );
    function->template communicate< Edge,   Vertex >( level );
    function->template communicate< Edge,   Face   >( level );
    function->template communicate< Face,   Edge   >( level );
  }

} // namespace hhg
