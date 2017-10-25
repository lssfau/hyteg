
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/indexing/Common.hpp"

namespace hhg {
namespace indexing {

using walberla::uint_t;

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
}
