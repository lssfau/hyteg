#pragma once

#include "simplex.hpp"

namespace hyteg {
namespace adaptiveRefinement {

// Factory for constructing sub-simplices of a given tetrahedron
class SimplexFactory
{
 public:
   using FaceMap = std::map< Idx< 3 >, std::shared_ptr< Simplex2 > >;
   using EdgeMap = std::map< Idx< 2 >, std::shared_ptr< Simplex1 > >;

   /* Ctor
      @param parent     parent-tetrahedron
      @param vertices   set of vertices for sub-simplices (i.e. global_vertices \intersect parent),
                        given by their global vertex-id
      @param edges      set of edges for sub-simplices s.th.
                        edges(i,j) = conv{vertices[i], vertices[j]}
      @param faces      set of faces for sub-simplices s.th.
                        faces(i,j,k) = conv{vertices[i], vertices[j], vertices[k]}
   */
   SimplexFactory( const std::shared_ptr< Simplex3 > parent,
                   const std::vector< int >&         vertices ) //, const EdgeMap& edges, const FaceMap& faces)
   : _parent( parent )
   , _vertices( vertices )
   // , _edges(edges)
   // , _faces(faces)
   {
      // _vertices.fill(-1);
   }

   /* add existing face \\subset parent for creating sub-tetrahedra
      @param a    index < #{v \\in vertices | v \\in parent}
      @param b    index < #{v \\in vertices | v \\in parent}
      @param c    index < #{v \\in vertices | v \\in parent}
   */
   inline void add_face( int a, int b, int c, std::shared_ptr< Simplex2 > face ) { _faces[{ a, b, c }] = face; }
   /* add existing edge \\subset parent for creating sub-faces/-tetrahedra
      @param a    index < #{v \\in vertices | v \\in parent}
      @param b    index < #{v \\in vertices | v \\in parent}
   */
   inline void add_edge( int a, int b, std::shared_ptr< Simplex1 > edge ) { _edges[{ a, b }] = edge; }

   // void update_faces(const FaceMap& faces) { _faces = faces; }
   // void update_edges(const EdgeMap& edges) { _edges = edges; }

   /* create new edge \\subset interior(parent) and add it to list of sub-edges
      @param a    index < #{v \\in vertices | v \\in parent}
      @param b    index < #{v \\in vertices | v \\in parent}
      @return pointer to edge = conv{vertices[a], vertices[b]}
   */
   std::shared_ptr< Simplex1 > make_edge( int a, int b, Color color = RED );
   /* create new face \\subset interior(parent) and add it to list of sub-faces
      @param a    index < #{v \\in vertices | v \\in parent}
      @param b    index < #{v \\in vertices | v \\in parent}
      @param c    index < #{v \\in vertices | v \\in parent}
      @return pointer to face = conv{vertices[a], vertices[b], vertices[c]}
   */
   std::shared_ptr< Simplex2 > make_face( int a, int b, int c );
   /* create new cell \\subset parent
      @param a    index < #{v \\in vertices | v \\in parent}
      @param b    index < #{v \\in vertices | v \\in parent}
      @param c    index < #{v \\in vertices | v \\in parent}
      @param d    index < #{v \\in vertices | v \\in parent}
      @return pointer to cell = conv{vertices[a], vertices[b], vertices[c], vertices[d]}
   */
   std::shared_ptr< Simplex3 > make_cell( int a, int b, int c, int d );

 private:
   const std::shared_ptr< Simplex3 > _parent;
   std::vector< int >                _vertices;
   // std::array<int, 10> _vertices;
   EdgeMap _edges;
   FaceMap _faces;
};

} // namespace adaptiveRefinement
} // namespace hyteg