
#pragma once

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg {
namespace indexing {
namespace vertexdof {

// ##############
// ### Common ###
// ##############

// Conversion functions: level to width (i.e. number of DoFs in the longest 'line')

template< uint_t level >
constexpr uint_t levelToWidth = levelinfo::num_microvertices_per_edge( level );


// ##################
// ### Macro Edge ###
// ##################

namespace macroedge {

/// Index of a vertex DoF on a macro edge (only access to owned DoFs, no ghost layers).
template< uint_t level >
constexpr uint_t index( const uint_t & col )
{
  return ::hhg::indexing::macroEdgeIndex< levelToWidth< level > >( col );
};

/// Index of a vertex DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor data, 1 to access second neighbor, ...
template< uint_t level >
constexpr uint_t index( const uint_t & col, const uint_t & neighbor )
{
  return                      macroEdgeSize< levelToWidth< level >     >()
         + ( neighbor - 1 ) * macroEdgeSize< levelToWidth< level > - 1 >()
         + macroEdgeIndex< levelToWidth< level > - 1 >( col );
};

}

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

/// Direct access functions

/// Index of a vertex DoF on a macro face (only access to owned DoFs, no ghost layers).
template< uint_t level >
constexpr uint_t index( const uint_t & col, const uint_t & row )
{
  return macroFaceIndex< levelToWidth< level > >( col, row );
};

/// Index of a vertex DoF on a ghost layer of a macro face.
/// \param neighbor 0 or 1 for the respective neighbor
template< uint_t level >
constexpr uint_t index( const uint_t & col, const uint_t & row, const uint_t & neighbor )
{
  WALBERLA_ASSERT( neighbor <= 1 );

  return                      macroFaceSize< levelToWidth< level >     >()
         + ( neighbor - 1 ) * macroFaceSize< levelToWidth< level > - 1 >()
         + macroFaceIndex< levelToWidth< level > - 1 >( col, row );

};

// Stencil access functions

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
template< uint_t level >
constexpr uint_t indexFromVertex( const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_C:
    return index< level >( col    , row     );
  case sD::VERTEX_E:
    return index< level >( col + 1, row     );
  case sD::VERTEX_W:
    return index< level >( col - 1, row     );
  case sD::VERTEX_N:
    return index< level >( col    , row + 1 );
  case sD::VERTEX_S:
    return index< level >( col    , row - 1 );
  case sD::VERTEX_NW:
    return index< level >( col - 1, row + 1 );
  case sD::VERTEX_SE:
    return index< level >( col + 1, row - 1 );
  default:
    return std::numeric_limits< uint_t >::max();
    break;
  }
}

// Iterators

/// Iterator over the border of a vertex DoF macro face.
/// See \ref FaceBorderIterator for more information.
template< uint_t level >
using BorderIterator = FaceBorderIterator< levelToWidth< level > >;

}
}
}
}
