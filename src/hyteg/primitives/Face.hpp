/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <array>
#include <core/DataTypes.h>
#include <core/debug/CheckFunctions.h>
#include <hyteg/Math.hpp>
#include <hyteg/primitives/Primitive.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>
#include <hyteg/types/pointnd.hpp>
#include <hyteg/types/types.hpp>
#include <vector>

namespace hyteg {

class Vertex;
class Edge;

class Face : public Primitive
{
 public:
   friend class SetupPrimitiveStorage;
   friend class PrimitiveStorage;
   template < class K_Simplex >
   friend class adaptiveRefinement::K_Mesh;

   Face( const PrimitiveID&                  primitiveID,
         const std::array< PrimitiveID, 3 >& vertexIDs,
         const std::array< PrimitiveID, 3 >& edgeIDs,
         const std::array< int, 3 >&         edgeOrientation,
         const std::array< Point3D, 3 >&     coordinates )
   : Primitive( primitiveID )
   , edge_orientation( edgeOrientation )
   , coords( coordinates )
   {
      neighborVertices_.push_back( vertexIDs[0] );
      neighborVertices_.push_back( vertexIDs[1] );
      neighborVertices_.push_back( vertexIDs[2] );

      neighborEdges_.push_back( edgeIDs[0] );
      neighborEdges_.push_back( edgeIDs[1] );
      neighborEdges_.push_back( edgeIDs[2] );
      std::array< Point3D, 2 > B( { { coords[1] - coords[0], coords[2] - coords[0] } } );
      area = std::abs( 0.5 * math::det2( B ) );
   }

   Face( walberla::mpi::RecvBuffer& recvBuffer )
   : Primitive( recvBuffer )
   {
      deserializeSubclass( recvBuffer );
   }

   uint_t vertex_index( const PrimitiveID& vertex ) const;
   uint_t edge_index( const PrimitiveID& edge ) const;
   uint_t cell_index( const PrimitiveID& cell ) const;

   std::vector< PrimitiveID > adjacent_edges( const PrimitiveID& vertex ) const;
   PrimitiveID                get_vertex_opposite_to_edge( const PrimitiveID& edge ) const;
   PrimitiveID                getEdgeOppositeToVertex( const PrimitiveID& vertexID ) const;
   PrimitiveID                get_edge_between_vertices( const PrimitiveID& v0, const PrimitiveID& v1 ) const;

   friend std::ostream& operator<<( std::ostream& os, const Face& face );

   const std::array< int, 3 >&     getEdgeOrientation() const { return edge_orientation; }
   const std::array< Point3D, 3 >& getCoordinates() const { return coords; }
   real_t                          getArea() const { return area; }

   const PrimitiveID& getVertexID0() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 3 );
      return neighborVertices_[0];
   }
   const PrimitiveID& getVertexID1() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 3 );
      return neighborVertices_[1];
   }
   const PrimitiveID& getVertexID2() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 3 );
      return neighborVertices_[2];
   }

   const PrimitiveID& getEdgeID0() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborEdges(), 3 );
      return neighborEdges_[0];
   }
   const PrimitiveID& getEdgeID1() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborEdges(), 3 );
      return neighborEdges_[1];
   }
   const PrimitiveID& getEdgeID2() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborEdges(), 3 );
      return neighborEdges_[2];
   }

   const PrimitiveID& getCellID0() const
   {
      WALBERLA_CHECK_GREATER( getNumNeighborCells(), 0 );
      return neighborCells_[0];
   }
   const PrimitiveID& getCellID1() const
   {
      WALBERLA_CHECK_GREATER( getNumNeighborCells(), 1 );
      return neighborCells_[1];
   }

   /// Returns true if the data that belongs to the passed \ref PrimitiveDataID is allocated.
   /// \param index the \ref PrimitiveDataID of the data that shall be asked for
   template < typename DataType >
   bool hasData( const PrimitiveDataID< DataType, Face >& index ) const
   {
      return genericHasData< DataType >( index );
   }

   /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
   /// \param index the \ref PrimitiveDataID of the data that should be returned
   template < typename DataType >
   DataType* getData( const PrimitiveDataID< DataType, Face >& index ) const
   {
      return genericGetData< DataType >( index );
   }

   using Primitive::getData;

   virtual void getLowerDimNeighbors( std::vector< PrimitiveID >& lowerDimNeighbors ) const
   {
      getNeighborEdges( lowerDimNeighbors );
   }
   virtual void getHigherDimNeighbors( std::vector< PrimitiveID >& higherDimNeighbors ) const { higherDimNeighbors.clear(); }

   virtual const std::vector< PrimitiveID >& getLowerDimNeighbors() const { return neighborEdges(); }
   virtual const std::vector< PrimitiveID >& getHigherDimNeighbors() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborFaces(), 0 );
      return neighborFaces();
   }

   virtual uint_t getNumLowerDimNeighbors() const { return getNumNeighborEdges(); }
   virtual uint_t getNumHigherDimNeighbors() const { return 0; }

   /// Returns indirect neighbor macro-face IDs (not including self). An indirect neighbor face is a face
   /// that shares at least one macro-vertex with this face.
   ///
   /// The indirect neighbors are sorted by the corresponding local interface macro-edge IDs.
   const std::map< uint_t, PrimitiveID >& getIndirectNeighborFaceIDs() const { return indirectNeighborFaceIDs_; }

 protected:
   /// Not public in order to guarantee that data is only added through the governing structure.
   /// This ensures valid DataIDs.
   template < typename DataType, typename DataHandlingType >
   inline void addData( const PrimitiveDataID< DataType, Face >& index, const std::shared_ptr< DataHandlingType >& dataHandling )
   {
      genericAddData( index, dataHandling, this );
   }

   virtual void serializeSubclass( walberla::mpi::SendBuffer& sendBuffer ) const;
   virtual void deserializeSubclass( walberla::mpi::RecvBuffer& recvBuffer );

 private:
   real_t               area;
   std::array< int, 3 > edge_orientation;

   std::array< Point3D, 3 > coords;

   void                            addCell( const PrimitiveID& cellID ) { neighborCells_.push_back( cellID ); }
   std::map< uint_t, PrimitiveID > indirectNeighborFaceIDs_;
};

} // namespace hyteg
