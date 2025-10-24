/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include <core/DataTypes.h>
#include <hyteg/mesh/MeshInfo.hpp>
#include <hyteg/primitives/Primitive.hpp>
#include <hyteg/types/types.hpp>
#include <vector>

#include "hyteg/types/PointND.hpp"

namespace hyteg {

class Edge;
class Face;

/// \brief  Macro-Vertex primitive
/// \author Daniel Drzisga (drzisga@ma.tum.de)
/// \date   March, 2017
///
/// The Vertex class represents a Macro-Vertex primitve. It saves geometrical and topological information
/// as well as pointers to memory reserved for \ref Function and \ref Operator.
class Vertex : public Primitive
{
 public:
   friend class SetupPrimitiveStorage;
   friend class PrimitiveStorage;
   template < class K_Simplex >
   friend class adaptiveRefinement::K_Mesh;

   /// Constructs a vertex with given id and coordinates
   /// \param primitiveID Id of vertex
   /// \param coordinates Spatial coordinates of vertex
   Vertex( const PrimitiveID& primitiveID, const Point3D& coordinates );

   Vertex( walberla::mpi::RecvBuffer& recvBuffer )
   : Primitive( recvBuffer )
   {
      deserializeSubclass( recvBuffer );
   }

   /// Returns VERTEX
   virtual PrimitiveTypeEnum getType() const final { return VERTEX; };

   /// Returns the Vertex-local index of \p edge
   /// \param edge Edge
   /// \returns Vertex-local index of \p edge
   uint_t edge_index( const PrimitiveID& edge ) const;

   /// Returns the Vertex-local index of \p face
   /// \param face Face
   /// \returns Vertex-local index of \p face
   uint_t face_index( const PrimitiveID& face ) const;

   /// Returns the Vertex-local index of \p cell
   /// \param cell Cell
   /// \returns Vertex-local index of \p cell
   uint_t cell_index( const PrimitiveID& cell ) const;

   /// Method overload for string formatting
   friend std::ostream& operator<<( std::ostream& os, const Vertex& vertex );

   /// Spatial coordinates of vertex
   const Point3D getCoordinates() const { return coordinates_; }

   /// Returns true if the data that belongs to the passed \ref PrimitiveDataID is allocated.
   /// \param index the \ref PrimitiveDataID of the data that shall be asked for
   template < typename DataType >
   bool hasData( const PrimitiveDataID< DataType, Vertex >& index ) const
   {
      return genericHasData< DataType >( index );
   }

   /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
   /// \param index the \ref PrimitiveDataID of the data that should be returned
   template < typename DataType >
   DataType* getData( const PrimitiveDataID< DataType, Vertex >& index ) const
   {
      return genericGetData< DataType >( index );
   }

   using Primitive::getData;

   virtual void getLowerDimNeighbors( std::vector< PrimitiveID >& lowerDimNeighbors ) const { lowerDimNeighbors.clear(); }
   virtual void getHigherDimNeighbors( std::vector< PrimitiveID >& higherDimNeighbors ) const
   {
      getNeighborEdges( higherDimNeighbors );
   }

   virtual const std::vector< PrimitiveID >& getLowerDimNeighbors() const
   {
      WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 0 );
      return neighborVertices();
   }
   virtual const std::vector< PrimitiveID >& getHigherDimNeighbors() const { return neighborEdges(); }

   virtual uint_t getNumLowerDimNeighbors() const { return 0; }
   virtual uint_t getNumHigherDimNeighbors() const { return getNumNeighborEdges(); }

 protected:
   /// Not public in order to guarantee that data is only added through the governing structure.
   /// This ensures valid DataIDs.
   template < typename DataType, typename DataHandlingType >
   inline void addData( const PrimitiveDataID< DataType, Vertex >& index,
                        const std::shared_ptr< DataHandlingType >& dataHandling )
   {
      genericAddData( index, dataHandling, this );
   }

   /// Deletes the data that belongs to the passed \ref PrimitiveDataID.
   /// Not public in order to guarantee that data is only deleted through the governing structure.
   /// \param index the \ref PrimitiveDataID of the data that shall be deleted
   template < typename DataType >
   void deleteData( const PrimitiveDataID< DataType, Vertex >& index )
   {
      return genericDeleteData< DataType >( index );
   }

   virtual void serializeSubclass( walberla::mpi::SendBuffer& sendBuffer ) const;
   virtual void deserializeSubclass( walberla::mpi::RecvBuffer& recvBuffer );

 private:
   void addEdge( const PrimitiveID& edgeID ) { neighborEdges_.push_back( edgeID ); }
   void addFace( const PrimitiveID& faceID ) { neighborFaces_.push_back( faceID ); }
   void addCell( const PrimitiveID& cellID ) { neighborCells_.push_back( cellID ); }

   Point3D coordinates_;
};

} // namespace hyteg
