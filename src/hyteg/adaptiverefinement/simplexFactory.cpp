
#include "simplexFactory.hpp"

namespace hyteg {
namespace adaptiveRefinement {

std::shared_ptr< Simplex1 > SimplexFactory::make_edge( int a, int b, Color color )
{
   auto edge = std::make_shared< Simplex1 >( _vertices[a], _vertices[b], nullptr, color );
   this->add_edge( a, b, edge );
   return edge;
}

std::shared_ptr< Simplex2 > SimplexFactory::make_face( int a, int b, int c )
{
   std::array< int64_t, 3 >                     V{ _vertices[a], _vertices[b], _vertices[c] };
   std::array< std::shared_ptr< Simplex1 >, 3 > E{ _edges[{ a, b }], _edges[{ b, c }], _edges[{ c, a }] };

   auto face = std::make_shared< Simplex2 >( V, E );
   this->add_face( a, b, c, face );
   return face;
}

std::shared_ptr< Simplex3 > SimplexFactory::make_cell( int a, int b, int c, int d )
{
   std::array< int64_t, 4 > V{ _vertices[a], _vertices[b], _vertices[c], _vertices[d] };

   std::array< std::shared_ptr< Simplex1 >, 6 > E{
       _edges[{ a, b }], _edges[{ b, c }], _edges[{ c, a }], _edges[{ a, d }], _edges[{ b, d }], _edges[{ c, d }] };

   std::array< std::shared_ptr< Simplex2 >, 4 > F{
       _faces[{ a, b, c }], _faces[{ a, b, d }], _faces[{ b, c, d }], _faces[{ a, c, d }] };

   return std::make_shared< Simplex3 >( V, E, F, _parent );
}

} // namespace adaptiveRefinement
} // namespace hyteg