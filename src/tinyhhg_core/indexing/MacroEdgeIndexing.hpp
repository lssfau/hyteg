
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/indexing/Common.hpp"

namespace hhg {
namespace indexing {

using walberla::uint_t;

namespace layout {

/// Required memory for the linear macro edge layout
template< uint_t width >
constexpr uint_t linearMacroEdgeSize()
{
  return width;
}

/// General linear memory layout indexing function for macro edges
template< uint_t width >
constexpr uint_t linearMacroEdgeIndex( const uint_t & col )
{
  return col;
}

}

template< uint_t width >
constexpr uint_t macroEdgeSize()
{
  return layout::linearMacroEdgeSize< width >();
}

template< uint_t width >
constexpr uint_t macroEdgeIndex( const uint_t & col )
{
  return layout::linearMacroEdgeIndex< width >( col );
}

}
}
