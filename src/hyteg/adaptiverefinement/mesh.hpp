/*
 * Copyright (c) 2021-2022 Benjamin Mann
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

// number of elements in the refined mesh by type of origin
struct RefinedElements
{
   uint_t n_U; // unrefined elements w.r.t. previous mesh
   uint_t n_R; // resulting from red refinement step
   uint_t n_G; // resulting from green refinement step (including green elements from the previous mesh)
   // return n_el = n_U + n_R + n_G
   inline uint_t n_el() const { return n_U + n_R + n_G; }
};

enum Loadbalancing
{
   ROUND_ROBIN,  // cheap loadbalancer
   GREEDY,       // assign elements to processes using greedy algorithm
   SFC,          // use space filling curve for loadbalancing (NOT IMPLEMENTED YET)
   SFC_VISUALIZE // use space filling curve for loadbalancing. Also export the sfc to vtk (NOT IMPLEMENTED YET)
};

enum RefinementStrategy
{
   /* refine all elements j where e_j >= (∑_i w_i e_i)/(∑_i w_i) with w_i = (n-i)^p,
      i.e., p=0: standard mean, p→∞: refine T_0, p→-∞: refine all elements */
   WEIGHTED_MEAN,
   PERCENTILE // refine the p*n elements where the error is largest
};

using ErrorVector = std::vector< std::pair< real_t, PrimitiveID > >;

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

   /* apply red-green refinement to this mesh
      @param el_to_refine     subset of elements that shall be refined (red)
                              given by primitiveIDs w.r.t. K_Mesh::make_storage()
      @param n_el_max         upper bound for number of elements in refined mesh
      @return |R\U|/|R| where R and U are the subsets of T which are marked for refinement and
                              remain unrefined, respectively. Note that the result will be 1
                              unless n_el_max has been reached.
   */
   real_t refineRG( const std::vector< PrimitiveID >& el_to_refine, uint_t n_el_max = std::numeric_limits< uint_t >::max() );

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param criterion        criterion w.r.t. sorted (greatest first) global error list whether an element should be refined
      @param n_el_max         upper bound for number of elements in refined mesh
      @return |R\U|/|R| where R and U are the subsets of T which are marked for refinement and
                              remain unrefined, respectively. Note that the result will be 1
                              unless n_el_max has been reached.
   */
   real_t refineRG( const ErrorVector&                                         local_errors,
                    const std::function< bool( const ErrorVector&, uint_t ) >& criterion,
                    uint_t                                                     n_el_max = std::numeric_limits< uint_t >::max() );

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param strategy         predefined refinement strategy
      @param p                parameter for refinement strategy
                              ratio*n_elements elements with the largest error will be refined
      @param n_el_max         upper bound for number of elements in refined mesh
      @param verbose          print information about refinement
      @return |R\U|/|R| where R and U are the subsets of T which are marked for refinement and
                              remain unrefined, respectively. Note that the result will be 1
                              unless n_el_max has been reached.
   */
   real_t refineRG( const ErrorVector& local_errors,
                    RefinementStrategy strategy,
                    real_t             p,
                    bool               verbose  = false,
                    uint_t             n_el_max = std::numeric_limits< uint_t >::max() );

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
   const std::vector< Point3D >& get_vertices() const { return _vertices; }
   // export mesh to gmsh file
   void exportMesh( const std::string& filename ) const;

 private:
   /* apply round robin loadbalancing to volume elements
   */
   std::vector< uint_t > loadbalancing_roundRobin( const bool allow_split_siblings );

   /* apply greedy loadbalancing to volume elements
   */
   std::vector< uint_t > loadbalancing_greedy( const std::map< PrimitiveID, Neighborhood >& nbrHood,
                                               const bool                                   allow_split_siblings );

   /* remove green edges from _T and replace them with their parents
   */
   void remove_green_edges();

   /* find all elements in U which require a red refinement step
      @param U set of unprocessed elements
      @param vtxs_added Flag that will be set to true if any new vertices are added during this step
      @return set R of elements requiring red refinement
   */
   std::set< std::shared_ptr< K_Simplex > > find_elements_for_red_refinement( const std::set< std::shared_ptr< K_Simplex > >& U,
                                                                              bool& vtxs_added );

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

   /* predict the number of additional elements that will be
      added during green refinement step
   */
   uint_t predict_n_el_green( const std::set< std::shared_ptr< K_Simplex > >& U ) const;

   /* find elements in _T corresponding to primitiveIDs
      @param primitiveIDs  set of primitiveIDs w.r.t. this->make_storage
      @return subset of _T for red refinement
   */
   std::vector< std::shared_ptr< K_Simplex > > init_R( const std::vector< PrimitiveID >& primitiveIDs ) const;

   /* extract connectivity, geometrymap, boundaryFlags, etc. from all elements*/
   void extract_data( std::map< PrimitiveID, VertexData >&   vtxData,
                      std::map< PrimitiveID, EdgeData >&     edgeData,
                      std::map< PrimitiveID, FaceData >&     faceData,
                      std::map< PrimitiveID, CellData >&     cellData,
                      std::map< PrimitiveID, Neighborhood >& nbrHood ) const;

   /* extract geometrymap, boundaryFlags, etc. from all elements*/
   void extract_data( std::map< PrimitiveID, VertexData >& vtxData,
                      std::map< PrimitiveID, EdgeData >&   edgeData,
                      std::map< PrimitiveID, FaceData >&   faceData,
                      std::map< PrimitiveID, CellData >&   cellData ) const;

   /* update target rank for all primitives */
   void update_targetRank( const std::map< PrimitiveID, VertexData >& vtxData,
                           const std::map< PrimitiveID, EdgeData >&   edgeData,
                           const std::map< PrimitiveID, FaceData >&   faceData,
                           const std::map< PrimitiveID, CellData >&   cellData );

   /* create PrimitiveStorage from SimplexData */
   std::shared_ptr< PrimitiveStorage > make_localPrimitives( std::map< PrimitiveID, VertexData >& vtxs,
                                                             std::map< PrimitiveID, EdgeData >&   edges,
                                                             std::map< PrimitiveID, FaceData >&   faces,
                                                             std::map< PrimitiveID, CellData >&   cells );

   /// @brief create sorted global error vector from local error vectors
   /// @param err_loc         local error vectors
   /// @param err_glob_sorted container for output data
   void gatherGlobalError( const ErrorVector& err_loc, ErrorVector& err_glob_sorted ) const;

   uint_t                                                  _n_vertices;
   uint_t                                                  _n_elements;
   uint_t                                                  _n_processes; // number of processes
   std::vector< Point3D >                                  _vertices;    // vertex coordinates
   std::vector< VertexData >                               _V;           // set of vertices of current refinement
   std::set< std::shared_ptr< K_Simplex > >                _T;           // set of elements of current refinement level
   std::map< PrimitiveID, std::shared_ptr< GeometryMap > > _geometryMap; // geometrymaps of original mesh

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

   /* apply red-green refinement to this mesh
      @param elements_to_refine  subset of elements that shall be refined (red)
                                 given by primitiveIDs w.r.t. K_Mesh::make_storage()
      @param n_el_max            upper bound for number of elements in refined mesh
      @return |R\U|/|R| where R and U are the subsets of T which are marked for refinement and
                              remain unrefined, respectively. Note that the result will be 1
                              unless n_el_max has been reached.
   */
   real_t refineRG( const std::vector< PrimitiveID >& elements_to_refine, uint_t n_el_max = std::numeric_limits< uint_t >::max() )
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->refineRG( elements_to_refine, n_el_max );
      }
      else
      {
         return _mesh2D->refineRG( elements_to_refine, n_el_max );
      }
   }

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param criterion        criterion whether an element should be refined
      @param n_el_max         upper bound for number of elements in refined mesh
      @return |R\U|/|R| where R and U are the subsets of T which are marked for refinement and
                              remain unrefined, respectively. Note that the result will be 1
                              unless n_el_max has been reached.
   */
   real_t refineRG( const ErrorVector&                                         local_errors,
                    const std::function< bool( const ErrorVector&, uint_t ) >& criterion,
                    uint_t                                                     n_el_max = std::numeric_limits< uint_t >::max() )
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->refineRG( local_errors, criterion, n_el_max );
      }
      else
      {
         return _mesh2D->refineRG( local_errors, criterion, n_el_max );
      }
   }

   /* apply red-green refinement to this mesh
      @param local_errors     list of elementwise errors for all local macro cells/faces
      @param strategy         predefined refinement strategy
      @param p                parameter for refinement strategy
                              ratio*n_elements elements with the largest error will be refined
      @param n_el_max         upper bound for number of elements in refined mesh
      @param verbose          print information about refinement
      @return |R\U|/|R| where R and U are the subsets of T which are marked for refinement and
                              remain unrefined, respectively. Note that the result will be 1
                              unless n_el_max has been reached.
   */
   real_t refineRG( const ErrorVector& local_errors,
                    RefinementStrategy strategy,
                    real_t             p,
                    bool               verbose  = false,
                    uint_t             n_el_max = std::numeric_limits< uint_t >::max() )
   {
      if ( _DIM == 3 )
      {
         return _mesh3D->refineRG( local_errors, strategy, p, verbose, n_el_max );
      }
      else
      {
         return _mesh2D->refineRG( local_errors, strategy, p, verbose, n_el_max );
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
   const std::vector< Point3D >& get_vertices() const
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
