#pragma once

#include <set>

#include "simplex.hpp"

namespace adaptiveRefinement {

template < class K_Simplex >
class Mesh
{
 public:
   Mesh( const std::vector< Point3D >& vertices, const std::set< std::shared_ptr< K_Simplex > >& elements )
   : _vertices( vertices )
   , _T( elements )
   {}

   /* apply red-green refinement to this mesh
      @param elements_to_refine  subset of elements eligible for red refinement
   */
   void refineRG( const std::set< std::shared_ptr< K_Simplex > >& elements_to_refine );

   const std::vector< Point3D >& vertices() const { return _vertices; }

   const std::set< std::shared_ptr< K_Simplex > >& elements() const { return _T; }

   // get minimum and maximum angle of the elements in T
   std::pair< real_t, real_t > min_max_angle() const;

   // compute total volume of the triangulated domain
   real_t volume() const;

 private:
   /* remove green edges from _T and replace the corresponding faces in R with their parents
   */
   void remove_green_edges( std::set< std::shared_ptr< K_Simplex > >& R );

   /* find all elements in U which require a red refinement step
      @param U set of unprocessed elements
      @return set R of elements requiring red refinement
   */
   std::set< std::shared_ptr< K_Simplex > > find_elements_for_red_refinement( const std::set< std::shared_ptr< K_Simplex > >& U );

   /*
      apply red refinement to all elements in R and remove them from U
      @param R set of elements marked for refinement
      @param U set of elements which have not been subject to refinement yet
      @return R_refined
   */
   std::set< std::shared_ptr< K_Simplex > > refine_red( const std::set< std::shared_ptr< K_Simplex > >& R,
                                                        std::set< std::shared_ptr< K_Simplex > >&       U );

   /* apply green refinement to all elements in U which
      have a new vertex on one of their edges,
      remove these elements from U and
      add the new elements to _T
      @param U set of elements which have not been subject to refinement yet
      @return U_refined = set newly refined elements
   */
   std::set< std::shared_ptr< K_Simplex > > refine_green( std::set< std::shared_ptr< K_Simplex > >& U );

   /* apply red refinement to element
      @return sub-elements
   */
   std::set< std::shared_ptr< K_Simplex > > refine_element_red( std::shared_ptr< K_Simplex > element );

   std::vector< Point3D >                   _vertices;
   std::set< std::shared_ptr< K_Simplex > > _T; // set of elements of current refinement level
};

using Mesh2D = Mesh< Simplex2 >;
using Mesh3D = Mesh< Simplex3 >;

} // namespace adaptiveRefinement
