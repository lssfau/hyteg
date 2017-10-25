
#pragma once

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg {
namespace indexing {

// TODO: All of this could also be pulled to the respective DoF space files / namespaces

// ##############
// ### Common ###
// ##############

// Conversion functions: level to width (i.e. number of DoFs in the longest 'line')

template< uint_t level >
constexpr uint_t levelToWidthVertexDoF = levelinfo::num_microvertices_per_edge( level );

template< uint_t level >
constexpr uint_t levelToWidthAnyEdgeDoF = levelinfo::num_microedges_per_edge( level );

template< uint_t level >
constexpr uint_t levelToFaceSizeAnyEdgeDoF = levelinfo::num_microedges_per_face( level ) / 3;

// ##################
// ### Macro Edge ###
// ##################

/// Index of a vertex DoF on a macro edge.
/// \param neighbor 0 to access the owned data, 1 to access first neighbor, ...
template< uint_t level >
constexpr uint_t VertexDoFOnMacroEdgeIndex( const uint_t & col, const uint_t & neighbor )
{
  if ( neighbor == 0 )
  {
    return linearMacroEdgeIndex< levelToWidthVertexDoF< level > >( col );
  }
  else
  {
    return                      linearMacroEdgeSize< levelToWidthVertexDoF< level >     >
           + ( neighbor - 1 ) * linearMacroEdgeSize< levelToWidthVertexDoF< level > - 1 >
           + linearMacroEdgeIndex< levelToWidthVertexDoF< level > - 1 >( col );
  }
};

// ##################
// ### Macro Face ###
// ##################

/// Direct access functions

template< uint_t level >
constexpr uint_t VertexDoFOnMacroFaceIndex( const uint_t & col, const uint_t & row )
{
  return linearMacroFaceIndex< levelToWidthVertexDoF< level > >( col, row );
};

// Stencil access functions

template< uint_t level >
constexpr uint_t VertexDoFOnMacroFaceIndexFromVertex( const uint_t & col, const uint_t & row,
                                                      const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_C:
    return VertexDoFOnMacroFaceIndex< level >( col, row );
  case sD::VERTEX_E:
    return VertexDoFOnMacroFaceIndex< level >( col + 1, row );

    // ...

  default:
    return 0;
    break;
  }
}

// Iterators

template< uint_t level >
using VertexDoFFaceBorderIterator = FaceBorderIterator< levelToWidthVertexDoF< level > >;

}
}
