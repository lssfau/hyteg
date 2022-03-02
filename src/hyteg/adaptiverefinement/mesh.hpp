/*
 * Copyright (c) 2021 Benjamin Mann
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
      @param elements_to_refine  subset of elements that shall be refined (red)
                                 given by primitiveIDs w.r.t. K_Mesh::make_storage()
   */
   // void refineRG( const std::vector< PrimitiveID >& elements_to_refine );

   /* apply red-green refinement to this mesh
      @param elements_to_refine  subset of elements that shall be refined (red)
                                 given by primitiveIDs w.r.t. K_Mesh::make_storage()
      @param n_el_max            upper bound for number of elements in refined mesh
   */
   void refineRG( const std::vector< PrimitiveID >& elements_to_refine, uint_t n_el_max = uint_t( -1 ) );

   // get minimum and maximum angle of the elements in T
   std::pair< real_t, real_t > min_max_angle() const;

   // compute total volume of the triangulated domain
   real_t volume() const;

   // construct PrimitiveStorage corresponding to current refinement
   std::shared_ptr< PrimitiveStorage > make_storage();

   inline uint_t n_elements() const { return _n_elements; }
   inline uint_t n_vtx() const { return _n_vertices; }

   // reference to set of volume elements corresponding to current refinement
   const std::set< std::shared_ptr< K_Simplex > >& get_elements() const { return _T; }
   // reference to list of vertex coordinates corresponding to current refinement
   const std::vector< Point3D >& get_vertices() const { return _vertices; }

 private:
   /* remove green edges from _T and replace them with their parents
   */
   void remove_green_edges();

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

   /* predict the number of additional elements that will be
      added during green refinement step
   */
   uint_t predict_n_el_green( const std::set< std::shared_ptr< K_Simplex > >& U ) const;

   /* find elements in _T corresponding to primitiveIDs
      @param primitiveIDs  set of primitiveIDs w.r.t. this->make_storage
      @return subset of _T for red refinement
   */
   std::vector< std::shared_ptr< K_Simplex > > init_R( const std::vector< PrimitiveID >& primitiveIDs ) const;

   /* extract connectivity, geometrymap and boundaryFlags from all elements and add PrimitiveIDs*/
   void extract_data( std::vector<EdgeData>& edgeData, std::vector<FaceData>& faceData, std::vector<CellData>& cellData ) const;

   /* create PrimitiveStorage and return id of first volume element */
   std::pair< uint_t, std::shared_ptr< PrimitiveStorage > >
       convert_to_storage( const std::vector<EdgeData>& edges, const std::vector<FaceData>& faces, const std::vector<CellData>& cells, const uint_t& n_processes );

   /* create SetupPrimitiveStorage*/
   std::shared_ptr< SetupPrimitiveStorage > make_setupStorage( const std::vector< uint_t >& vertexIDs,
                                                               const EdgeData&              edges,
                                                               const FaceData&              faces,
                                                               const CellData&              cells,
                                                               const uint_t&                n_processes );

   uint_t                                             _n_vertices;
   uint_t                                             _n_elements;
   uint_t                                             _n_processes;        // number of processes
   std::vector< Point3D >                             _vertices;           // vertex coordinates
   std::vector< uint_t >                              _vertexGeometryMap;  // geometrymap for vertices
   std::vector< uint_t >                              _vertexBoundaryFlag; // boundaryFlag for vertices
   std::set< std::shared_ptr< K_Simplex > >           _T;                  // set of elements of current refinement level
   std::map< uint_t, std::shared_ptr< GeometryMap > > _geometryMap;        // geometrymaps of original mesh
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

   // /* apply red-green refinement to this mesh
   //    @param elements_to_refine  subset of elements that shall be refined (red)
   //                               given by primitiveIDs w.r.t. K_Mesh::setupStorage()
   // */
   // void refineRG( const std::vector< PrimitiveID >& elements_to_refine )
   // {
   //    if ( _DIM == 3 )
   //    {
   //       return _mesh3D->refineRG( elements_to_refine );
   //    }
   //    else
   //    {
   //       return _mesh2D->refineRG( elements_to_refine );
   //    }
   // }

   /* apply red-green refinement to this mesh
      @param elements_to_refine  subset of elements that shall be refined (red)
                                 given by primitiveIDs w.r.t. K_Mesh::make_storage()
      @param n_el_max            upper bound for number of elements in refined mesh
   */
   void refineRG( const std::vector< PrimitiveID >& elements_to_refine, uint_t n_el_max = uint_t( -1 ) )
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

   // construct PrimitiveStorage corresponding to current refinement
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

 private:
   uint_t                    _DIM;    // spacial dimension
   std::shared_ptr< Mesh2D > _mesh2D; // internal mesh object for the case _DIM=2
   std::shared_ptr< Mesh3D > _mesh3D; // internal mesh object for the case _DIM=3
};

} // namespace adaptiveRefinement
} // namespace hyteg
