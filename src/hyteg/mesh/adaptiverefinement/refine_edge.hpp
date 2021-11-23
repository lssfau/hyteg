#pragma once

#include "simplex.hpp"

namespace hyteg {
namespace adaptiveRefinement {

/* apply edge bisection
      @param edge the edge to be bisected
      @param vtx index of the vertex on the edge midpoint
   */
inline void bisect_edge( std::shared_ptr< Simplex1 > edge, int64_t vtx )
{
   WALBERLA_ASSERT( not edge->has_children() );
   edge->set_midpoint_idx( vtx );
   edge->add_child( std::make_shared< Simplex1 >( edge->get_vertices()[0], vtx, edge ) );
   edge->add_child( std::make_shared< Simplex1 >( vtx, edge->get_vertices()[1], edge ) );
}

} // namespace adaptiveRefinement
} // namespace hyteg
