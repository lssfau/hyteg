
#pragma once

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg {
namespace indexing {

// ##################
// ### Macro Edge ###
// ##################

template< uint_t level >
constexpr uint_t HorizontalEdgeDoFOnMacroEdgeIndex( const uint_t & col, const uint_t & neighbor )
{
  WALBERLA_ASSERT( false );
  return 0;
};

// ##################
// ### Macro Face ###
// ##################

/// Direct access functions

template< uint_t level >
constexpr uint_t HorizontalEdgeDoFOnMacroFaceIndex( const uint_t & col, const uint_t & row )
{
  return macroFaceIndex< levelToWidthAnyEdgeDoF< level > >( col, row );
};

// Stencil access functions



template< uint_t level >
constexpr uint_t EdgeDoFFaceIndexFromVertex( const uint_t & col, const uint_t & row,
                                             const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::EDGE_HO_C:
    return HorizontalEdgeDoFOnMacroFaceIndex< level >( col, row );

    // ...

  default:
    return 0;
    break;
  }
}



}
}
