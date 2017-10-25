
#pragma once

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg {
namespace indexing {
namespace vertexdof {

// TODO: All of this could also be pulled to the respective DoF space files / namespaces

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

/// Index of a vertex DoF on a macro edge.
/// \param neighbor 0 to access the owned data, 1 to access first neighbor, ...
template< uint_t level >
constexpr uint_t index( const uint_t & col, const uint_t & neighbor )
{
  if ( neighbor == 0 )
  {
    return hhg::indexing::macroEdgeIndex< levelToWidth< level > >( col );
  }
  else
  {
    return                      macroEdgeSize< levelToWidth< level >     >()
           + ( neighbor - 1 ) * macroEdgeSize< levelToWidth< level > - 1 >()
           + macroEdgeIndex< levelToWidth< level > - 1 >( col );
  }
};

}

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

/// Direct access functions

/// Index of a vertex DoF on a macro face.
/// \param neighbor 0 to access the owned data, 1 to access first neighbor, 2 to access second neighbor
template< uint_t level >
constexpr uint_t index( const uint_t & col, const uint_t & row, const uint_t & neighbor )
{
  WALBERLA_ASSERT( neighbor <= 2 );

  if ( neighbor == 0 )
  {
    return macroFaceIndex< levelToWidth< level > >( col, row );
  }
  else
  {
    return                      macroFaceSize< levelToWidth< level >     >()
           + ( neighbor - 1 ) * macroFaceSize< levelToWidth< level > - 1 >()
           + macroFaceIndex< levelToWidth< level > - 1 >( col, row );
  }
};

// Stencil access functions

template< uint_t level >
constexpr uint_t indexFromVertex( const uint_t & col, const uint_t & row,
                                                      const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_C:
    return index< level >( col, row, 0 );
  case sD::VERTEX_E:
    return index< level >( col + 1, row, 0 );

    // ...

  default:
    return 0;
    break;
  }
}

// Iterators

template< uint_t level >
using BorderIterator = FaceBorderIterator< levelToWidth< level > >;

}
}
}
}
