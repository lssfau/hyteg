/*
 * Copyright (c) 2021-2024 Benjamin Mann
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <set>

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "simplexData.hpp"

namespace hyteg {
namespace adaptiveRefinement {

enum Loadbalancing
{
   ROUND_ROBIN,  // cheap loadbalancer
   GREEDY,       // assign elements to processes using greedy algorithm
   SFC,          // use space filling curve for loadbalancing (NOT IMPLEMENTED YET)
   SFC_VISUALIZE // use space filling curve for loadbalancing. Also export the sfc to vtk (NOT IMPLEMENTED YET)
};

struct Strategy
{
   enum Type
   {
      /* e_mean(p) = (∑_i (n-i)^p e_i) / (∑_i (n-i)^p)
         e.g. e_mean(0) = (∑_i e_i)/n (arithmetic mean),
            e_mean(∞) = e_0 (greatest error),
            e_mean(-∞) = e_{n-1} (smallest error)
         Refinement:
            refine all elements j where e_j >= e_mean(p)
            fewer elements are refined as p gets larger
         Coarsening:
            coarsen all elements j where e_j < e_mean(p_c)
            more elements are coarsened as p gets larger
      */
      WEIGHTED_MEAN,
      /* Refinement:
            refine the p*n elements with largest error
         Coarsening:
            coarsen the p*n elements with smallest error
      */
      PERCENTILE,
      /* refine smallest subset R of T where ∑_t∈R e_t >= p ∑_t∈T e_t (p=0 ⇒ R=∅; p=1 ⇒ R=T)
         (when providing a list of squared errors, this is the classic Dörfler marking)
      */
      DOERFLER,
      /* coarsen a number of elements equal to p*|R|
         should be chosen to keep the number of elements approximately constant
      */
      MULTIPLE_OF_R
   };

   constexpr Strategy( Type type, real_t param )
   : t( type )
   , p( param )
   {}

   /* e_mean(p) = (∑_i (n-i)^p e_i) / (∑_i (n-i)^p)
      e.g. e_mean(0) = (∑_i e_i)/n (arithmetic mean),
         e_mean(∞) = e_0 (greatest error),
         e_mean(-∞) = e_{n-1} (smallest error)
      Refinement:
         refine all elements j where e_j >= e_mean(p)
         fewer elements are refined as p gets larger
      Coarsening:
         coarsen all elements j where e_j < e_mean(p_c)
         more elements are coarsened as p gets larger
      @param p weighting parameter; (-∞, ∞)
   */
   constexpr static Strategy mean( real_t p = real_t( 0.0 ) ) { return { WEIGHTED_MEAN, p }; }

   /* Refinement:
         refine the p*n elements with largest error
      Coarsening:
         coarsen the p*n elements with smallest error
      @param p proportion of elements to be marked; [0,1]
   */
   constexpr static Strategy percentile( real_t p ) { return { PERCENTILE, p }; }

   /* refine all elements j where e_j <= p*∑_j e_j
      (when providing a list of squared errors, this is the classic Dörfler marking)
      @param p proportion of total error; [0,1]
   */
   constexpr static Strategy doerfler( real_t p ) { return { DOERFLER, p }; }

   /* coarsen a number of elements p*|R| to enable keeping the number
      of elements approximately constant
   */
   constexpr static Strategy multiple_of_R( real_t p ) { return { MULTIPLE_OF_R, p }; }

   /* refine/coarsen all elements */
   constexpr static Strategy all() { return { PERCENTILE, real_t( 1.0 ) }; }

   /* refine/coarsen no elements */
   constexpr static Strategy none() { return { PERCENTILE, real_t( 0.0 ) }; }

   Type   t;
   real_t p;
};

// local error for each macro-cell
using ErrorVector = std::vector< std::pair< real_t, PrimitiveID > >;
// stores IDs of neighbor primitives
using Neighborhood = std::array< std::vector< PrimitiveID >, PrimitiveType::ALL >;
// stores neighborhood of all primitives
using NeighborhoodMap = std::array< std::map< PrimitiveID, Neighborhood >, PrimitiveType::ALL >;

// adaptively refinable mesh for K-dimensional domains
template < class K_Simplex >
class K_Mesh
{
 public:
   /* construct adaptable mesh from setupStorage
      @param setupStorage SetupPrimitiveStorage corresponding to initial coarse grid.
                           Note that Geometrymaps and boundary flags must be applied
                           to setupStorage before constructing the adaptive mesh.
                           Furthermore, afer calling this constructor, the original
                           setupStorage should not be used to construct a
                           PrimitiveStorage. Instead use this->make_storage().
   */
   K_Mesh( const SetupPrimitiveStorage& setupStorage );

   /* apply k steps of regular refinement
      @param k    number of refinement steps
   */
   void refine_uniform( uint_t k );

   /* apply red-green refinement to this mesh
      @param el_to_refine     subset of elements that shall be refined (red)
                              given by primitiveIDs w.r.t. K_Mesh::make_storage()
      @param el_to_coarsen    subset of elements that shall be coarsened. Note that
                              only previously refined elements can be coarsened this way
   */
   void refineRG( const std::vector< PrimitiveID >& el_to_refine, const std::vector< PrimitiveID >& el_to_coarsen = {} );

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param crit_r           criterion w.r.t. sorted (greatest first) global error list whether an element should be refined
      @param crit_c           criterion w.r.t. sorted (greatest first) global error list whether an element should be coarsened
   */
   void refineRG( const ErrorVector&                                         local_errors,
                  const std::function< bool( const ErrorVector&, uint_t ) >& crit_r,
                  const std::function< bool( const ErrorVector&, uint_t ) >& crit_c );

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param refinement       refinement strategy
      @param coarsening       coarsening strategy
      @param verbose          print information about refinement
   */
   void refineRG( const ErrorVector& local_errors, Strategy refinement, Strategy coarsening, bool verbose = false );

   // get minimum and maximum angle of the elements of T
   std::pair< real_t, real_t > min_max_angle() const;

   // get mean values of minimum and maximum angle over all elements of T
   std::pair< real_t, real_t > mean_min_max_angle() const;

   // compute total volume of the triangulated domain
   real_t volume() const;

   // compute minimum and maximum volume over all elements of T
   std::pair< real_t, real_t > min_max_volume() const;

   /* construct PrimitiveStorage corresponding to current refinement
   */
   std::shared_ptr< PrimitiveStorage > make_storage();

   /* apply loadbalancing scheme to current refinement
      @param lbScheme scheme used for load balancing
      @param migrationInfo_required if set to false, the return value will be empty
      @param allow_split_siblings if true, green siblings may be assigned to different processes (don't do this when interpolating between grids)
      @param verbose show information about distribution of volume elements over processes
      @return MigrationInfo to be used to migratePrimitives of the storage (in case loadbalancing is called after make_storage)
   */
   MigrationInfo loadbalancing( const Loadbalancing& lbScheme               = ROUND_ROBIN,
                                const bool           migrationInfo_required = false,
                                const bool           allow_split_siblings   = true,
                                const bool           verbose                = false );

   inline uint_t n_elements() const { return _n_elements; }
   inline uint_t n_vtx() const { return _n_vertices; }

   // reference to set of volume elements corresponding to current refinement
   const std::set< std::shared_ptr< K_Simplex > >& get_elements() const { return _T; }
   // reference to list of vertex coordinates corresponding to current refinement
   const EnumeratedList< Point3D >& get_vertices() const { return _coords; }
   // export mesh to gmsh file
   void exportMesh( const std::string& filename ) const;

 private:
   /* apply round robin loadbalancing to volume elements
   */
   std::vector< uint_t > loadbalancing_roundRobin( const bool allow_split_siblings );

   /* apply greedy loadbalancing to volume elements
   */
   std::vector< uint_t > loadbalancing_greedy( const NeighborhoodMap& nbrHood, const bool allow_split_siblings );

   /* remove green edges from _T and replace them with their parents.
      Also remove green children from all elements.
   */
   void remove_green_edges();

   /* refine all faces with more than one hanging node on their edges (3d only)
      @return true if any face was refined
   */
   bool refine_faces();

   /* find all elements which require a red refinement step
      @return set R of elements requiring red refinement
   */
   std::set< std::shared_ptr< K_Simplex > > find_elements_for_red_refinement() const;

   /*
      apply red refinement to all elements in R
      @param R set of elements marked for refinement
   */
   void refine_red( const std::set< std::shared_ptr< K_Simplex > >& R );

   /* apply green refinement to all elements in T with hanging nodes
   */
   void refine_green();

   /* unrefine all elements in P
      @param C set of elements marked for unrefinement, i.e.,
               each el in P should be a parent of an element in T
   */
   void unrefine( const std::set< std::shared_ptr< K_Simplex > >& P );

   /* generate sets R and Pc of elements marked for refinement and un-refinement, respectively
      @param id_r set of primitiveIDs of elements in T that shall be refined
      @param id_c set of primitiveIDs of elements in T that shall be coarsened
      @return [R,Pc] where R is a subset of T and Pc the set of parent elements of a subset C of T
   */
   std::pair< std::set< std::shared_ptr< K_Simplex > >, std::set< std::shared_ptr< K_Simplex > > >
       init_R_Pc( const std::vector< PrimitiveID >& id_r, const std::vector< PrimitiveID >& id_c ) const;

   /* check whether
      * all required vertices and only those are stored
      * elements in T have no children
      * all faces/edges of elements in T exist and have no children
      @param hanging_nodes_allowed  number of hanging nodes allowed on an edge
      @param unassigned_nodes_allowed  allow vertices which aren't part of any element
   */
   void check_integrity( uint_t hanging_nodes_allowed, bool unassigned_nodes_allowed = false ) const;

   /* extract geometryMap, boundaryFlags, etc. from all elements*/
   void extract_data( std::map< PrimitiveID, EdgeData >& edgeData,
                      std::map< PrimitiveID, FaceData >& faceData,
                      std::map< PrimitiveID, CellData >& cellData ) const;

   /* create PrimitiveStorage from SimplexData */
   std::shared_ptr< PrimitiveStorage > make_localPrimitives( std::map< PrimitiveID, EdgeData >& edges,
                                                             std::map< PrimitiveID, FaceData >& faces,
                                                             std::map< PrimitiveID, CellData >& cells );

   /// @brief create sorted global error vector from local error vectors
   /// @param err_loc         local error vectors
   /// @param err_glob_sorted container for output data
   void gatherGlobalError( const ErrorVector& err_loc, ErrorVector& err_glob_sorted ) const;

   static constexpr auto VOL = PrimitiveType( K_Simplex::TYPE );

   uint_t                                                  _n_vertices;
   uint_t                                                  _n_elements;
   uint_t                                                  _n_processes; // number of processes
   EnumeratedList< Point3D >                               _coords;      // vertex coordinates
   EnumeratedList< VertexData >                            _V;           // set of vertices of current refinement
   std::set< std::shared_ptr< K_Simplex > >                _T;           // set of elements of current refinement
   std::map< PrimitiveID, std::shared_ptr< GeometryMap > > _geometryMap; // geometryMaps of original mesh

   PrimitiveID _invalidID;
};

using Mesh2D = K_Mesh< Simplex2 >;
using Mesh3D = K_Mesh< Simplex3 >;

// wrapper for adaptively refinable mesh for domains with arbitrary dimension
class Mesh
{
 public:
   /* construct adaptable mesh from setupStorage
      @param setupStorage SetupPrimitiveStorage corresponding to initial coarse grid.
                           Note that Geometrymaps and boundary flags must be applied
                           to setupStorage before constructing the adaptive mesh.
                           Furthermore, afer calling this constructor, the original
                           setupStorage should not be used to construct a
                           PrimitiveStorage. Instead use this->make_storage().
   */
   Mesh( const SetupPrimitiveStorage& setupStorage )
   : _DIM( ( setupStorage.getNumberOfCells() > 0 ) ? 3 : 2 )
   {
      if ( _DIM == 3 )
      {
         _mesh2D = nullptr;
         _mesh3D = std::make_shared< Mesh3D >( setupStorage );
      }
      else
      {
         _mesh2D = std::make_shared< Mesh2D >( setupStorage );
         _mesh3D = nullptr;
      }
   }

   /* apply k steps of regular refinement
      @param k    number of refinement steps
   */
   void refine_uniform( uint_t k )
   {
      if ( _DIM == 3 )
      {
         _mesh3D->refine_uniform( k );
      }
      else
      {
         _mesh2D->refine_uniform( k );
      }
   }

   /* apply red-green refinement to this mesh
      @param el_to_refine     subset of elements that shall be refined (red)
                              given by primitiveIDs w.r.t. K_Mesh::make_storage()
      @param el_to_coarsen    subset of elements that shall be coarsened. Note that
                              only previously refined elements can be coarsened this way
   */
   void refineRG( const std::vector< PrimitiveID >& el_to_refine, const std::vector< PrimitiveID >& el_to_coarsen = {} )
   {
      if ( _DIM == 3 )
      {
         _mesh3D->refineRG( el_to_refine, el_to_coarsen );
      }
      else
      {
         _mesh2D->refineRG( el_to_refine, el_to_coarsen );
      }
   }

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param crit_r           criterion w.r.t. sorted (greatest first) global error list whether an element should be refined
      @param crit_c           criterion w.r.t. sorted (greatest first) global error list whether an element should be coarsened
   */
   void refineRG(
       const ErrorVector&                                         local_errors,
       const std::function< bool( const ErrorVector&, uint_t ) >& crit_r,
       const std::function< bool( const ErrorVector&, uint_t ) >& crit_c = []( const ErrorVector&, uint_t ) { return false; } )
   {
      if ( _DIM == 3 )
      {
         _mesh3D->refineRG( local_errors, crit_r, crit_c );
      }
      else
      {
         _mesh2D->refineRG( local_errors, crit_r, crit_c );
      }
   }

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param refinement       refinement strategy
      @param coarsening       coarsening strategy
      @param verbose          print information about refinement
   */
   void refineRG( const ErrorVector& local_errors,
                  Strategy           refinement,
                  Strategy           coarsening = Strategy::none(),
                  bool               verbose    = false )
   {
      if ( _DIM == 3 )
      {
         _mesh3D->refineRG( local_errors, refinement, coarsening, verbose );
      }
      else
      {
         _mesh2D->refineRG( local_errors, refinement, coarsening, verbose );
      }
   }

   // get minimum and maximum angle of the elements in T
   std::pair< real_t, real_t > min_max_angle() const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->min_max_angle();
      }
      else
      {
         return _mesh2D->min_max_angle();
      }
   }

   // compute total volume of the triangulated domain
   real_t volume() const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->volume();
      }
      else
      {
         return _mesh2D->volume();
      }
   }

   // get mean values of minimum and maximum angle over all elements of T
   std::pair< real_t, real_t > mean_min_max_angle() const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->mean_min_max_angle();
      }
      else
      {
         return _mesh2D->mean_min_max_angle();
      }
   }

   // compute minimum and maximum volume over all elements of T
   std::pair< real_t, real_t > min_max_volume() const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->min_max_volume();
      }
      else
      {
         return _mesh2D->min_max_volume();
      }
   }

   /* construct PrimitiveStorage corresponding to current refinement
   */
   std::shared_ptr< PrimitiveStorage > make_storage()
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->make_storage();
      }
      else
      {
         return _mesh2D->make_storage();
      }
   };

   /* apply loadbalancing scheme to current refinement
      @param lbScheme scheme used for load balancing
      @param migrationInfo_required if set to false, the return value will be empty
      @param allow_split_siblings if true, green siblings may be assigned to different processes (don't do this when interpolating between grids)
      @param verbose show information about distribution of volume elements over processes
      @return MigrationInfo to be used to migratePrimitives of the storage (in case loadbalancing is called after make_storage)
   */
   MigrationInfo loadbalancing( const Loadbalancing& lbScheme               = ROUND_ROBIN,
                                const bool           migrationInfo_required = false,
                                const bool           allow_split_siblings   = true,
                                const bool           verbose                = false )
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->loadbalancing( lbScheme, migrationInfo_required, allow_split_siblings, verbose );
      }
      else
      {
         return _mesh2D->loadbalancing( lbScheme, migrationInfo_required, allow_split_siblings, verbose );
      }
   }

   // get number of elements in current refinement
   inline uint_t n_elements() const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->n_elements();
      }
      else
      {
         return _mesh2D->n_elements();
      }
   }

   // get number of vertices in current refinement
   inline uint_t n_vtx() const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->n_vtx();
      }
      else
      {
         return _mesh2D->n_vtx();
      }
   }

   inline uint_t dim() const { return _DIM; }

   // export mesh to gmsh file
   void exportMesh( const std::string& filename ) const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->exportMesh( filename );
      }
      else
      {
         return _mesh2D->exportMesh( filename );
      }
   }

   // get all faces of the 2d mesh
   const std::set< std::shared_ptr< Simplex2 > >& get_elements2d() const { return _mesh2D->get_elements(); }

   // get all tets of the 3d mesh
   const std::set< std::shared_ptr< Simplex3 > >& get_elements3d() const { return _mesh3D->get_elements(); }

   // get vertex coordinates
   const EnumeratedList< Point3D >& get_vertices() const
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->get_vertices();
      }
      else
      {
         return _mesh2D->get_vertices();
      }
   }

 private:
   uint_t                    _DIM;    // spacial dimension
   std::shared_ptr< Mesh2D > _mesh2D; // internal mesh object for the case _DIM=2
   std::shared_ptr< Mesh3D > _mesh3D; // internal mesh object for the case _DIM=3
};

} // namespace adaptiveRefinement
} // namespace hyteg
