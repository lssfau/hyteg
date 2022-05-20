/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "core/logging/Logging.h"

#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/types/pointnd.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

class Cell : public Primitive
{
 public:
   friend class SetupPrimitiveStorage;
   template < class K_Simplex >
   friend class adaptiveRefinement::K_Mesh;

   /// Creates a macro-cell instance
   ///
   /// \param primitiveID ID of this macro-cell
   /// \param vertexIDs neighbor macro-vertex IDs
   /// \param edgeIDs neighbor macro-edge IDs
   /// \param faceIDs neighbor macro-face IDs
   /// \param coordinates absolute coordinates of the four vertices of this macro-cell
   /// \param edgeLocalVertexToCellLocalVertexMaps Maps for each edge of the macro-cell that map the local vertex ID (one of 0, 1) of the edge
   ///                                             to the corresponding local vertex ID (one of 0, 1, 2, 3) of this macro-cell. \n
   ///                                             Refer to the documentation for detailed illustrations of that mapping.
   /// \param faceLocalVertexToCellLocalVertexMaps Maps for each face of the macro-cell that map the local vertex ID (one of 0, 1, 2) of the face
   ///                                             to the corresponding local vertex ID (one of 0, 1, 2, 3) of this macro-cell. \n
   ///                                             Refer to the documentation for detailed illustrations of that mapping.
   Cell( const PrimitiveID&                                 primitiveID,
         const std::vector< PrimitiveID >&                  vertexIDs,
         const std::vector< PrimitiveID >&                  edgeIDs,
         const std::vector< PrimitiveID >&                  faceIDs,
         const std::array< Point3D, 4 >&                    coordinates,
         const std::array< std::map< uint_t, uint_t >, 6 >& edgeLocalVertexToCellLocalVertexMaps,
         const std::array< std::map< uint_t, uint_t >, 4 >& faceLocalVertexToCellLocalVertexMaps );

   Cell( walberla::mpi::RecvBuffer& recvBuffer )
   : Primitive( recvBuffer )
   {
      deserializeSubclass( recvBuffer );
   }

   virtual void getLowerDimNeighbors( std::vector< PrimitiveID >& lowerDimNeighbors ) const
   {
      getNeighborFaces( lowerDimNeighbors );
   }
   virtual void getHigherDimNeighbors( std::vector< PrimitiveID >& higherDimNeighbors ) const { higherDimNeighbors.clear(); }

   virtual const std::vector< PrimitiveID >& getLowerDimNeighbors() const { return neighborFaces(); }
   virtual const std::vector< PrimitiveID >& getHigherDimNeighbors() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborCells(), 0 );
      return neighborCells();
   }

   virtual uint_t getNumLowerDimNeighbors() const { return getNumNeighborFaces(); }
   virtual uint_t getNumHigherDimNeighbors() const { return 0; }

   /// Returns true if the data that belongs to the passed \ref PrimitiveDataID is allocated.
   /// \param index the \ref PrimitiveDataID of the data that shall be asked for
   template < typename DataType >
   bool hasData( const PrimitiveDataID< DataType, Cell >& index ) const
   {
      return genericHasData< DataType >( index );
   }

   /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
   /// \param index the \ref PrimitiveDataID of the data that should be returned
   template < typename DataType >
   DataType* getData( const PrimitiveDataID< DataType, Cell >& index ) const
   {
      return genericGetData< DataType >( index );
   }

   const std::array< Point3D, 4 >&                    getCoordinates() const { return coordinates_; }
   const std::array< std::map< uint_t, uint_t >, 6 >& getEdgeLocalVertexToCellLocalVertexMaps() const
   {
      return edgeLocalVertexToCellLocalVertexMaps_;
   }
   const std::array< std::map< uint_t, uint_t >, 4 >& getFaceLocalVertexToCellLocalVertexMaps() const
   {
      return faceLocalVertexToCellLocalVertexMaps_;
   }
   const Point3D getFaceInwardNormal( uint_t oppositeLocalVertexID ) const
   {
      return faceInwardNormals_.at( oppositeLocalVertexID );
   }

   /// Returns all macro-cell IDs of macro-cell that share at least one macro-face with this cell.
   /// Neighbor cells are mapped from local face IDs.
   const std::map< uint_t, PrimitiveID >& getIndirectNeighborCellIDsOverFaces() const
   {
      return indirectNeighborCellIDsOverFaces_;
   }

   /// Returns all macro-cell IDs of macro-cells that share at least one macro-vertex with this cell.
   const std::vector< PrimitiveID >& getIndirectNeighborCellIDsOverVertices() const
   {
      return indirectNeighborCellIDsOverVertices_;
   }

   uint_t getLocalVertexID( const PrimitiveID& vertexID ) const;
   uint_t getLocalEdgeID( const PrimitiveID& edgeID ) const;
   uint_t getLocalFaceID( const PrimitiveID& faceID ) const;

   const PrimitiveID& getOppositeVertexID( const PrimitiveID& faceID ) const
   {
      auto   otherLocalVertexIDs   = indexing::cellLocalFaceIDsToSpanningVertexIDs.at( getLocalFaceID( faceID ) );
      uint_t localOppositeVertexID = 6;
      for ( auto id : otherLocalVertexIDs )
      {
         localOppositeVertexID -= id;
      }
      return neighborVertices().at( localOppositeVertexID );
   }

   const PrimitiveID& getOppositeEdgeID( const PrimitiveID& edgeID ) const
   {
      return neighborEdges().at( indexing::getCellLocalOppositeEdgeID( getLocalEdgeID( edgeID ) ) );
   }

 protected:
   /// Not public in order to guarantee that data is only added through the governing structure.
   /// This ensures valid DataIDs.
   template < typename DataType, typename DataHandlingType >
   inline void addData( const PrimitiveDataID< DataType, Cell >& index, const std::shared_ptr< DataHandlingType >& dataHandling )
   {
      genericAddData( index, dataHandling, this );
   }

   virtual void serializeSubclass( walberla::mpi::SendBuffer& sendBuffer ) const;
   virtual void deserializeSubclass( walberla::mpi::RecvBuffer& recvBuffer );

 private:
   std::array< Point3D, 4 >                    coordinates_;
   std::array< std::map< uint_t, uint_t >, 6 > edgeLocalVertexToCellLocalVertexMaps_;
   std::array< std::map< uint_t, uint_t >, 4 > faceLocalVertexToCellLocalVertexMaps_;
   // stores all 4 face inward normals, array idx == opposite vertex
   std::array< Point3D, 4 > faceInwardNormals_;

   std::vector< PrimitiveID >      indirectNeighborCellIDsOverVertices_;
   std::map< uint_t, PrimitiveID > indirectNeighborCellIDsOverFaces_;
};

} // namespace hyteg
